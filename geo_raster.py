import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_bounds
import geopandas as gpd
import numpy as np
from pyproj import Proj
from matplotlib import pyplot as plt
from shapely.geometry import Polygon
from PIL import Image

ft_to_meters = 1.0/3.28084

class GeoRaster:
    def __init__(self, units_per_pixel,
                     proj, zone, north, datum, units):
        #Example: 
        #Proj("+proj=utm +zone=15, +north +datum=NAD83 +units=m +no_defs")
        self.geo_data = None
        self.geo_bounds = np.array([0,0,0,0])
        self.projection = proj
        self.zone = zone
        if north:
            self.pole = "north"
        else:
            self.pole = "south"        
        self.datum = datum
        self.units = units
        proj_str=f"+proj={self.projection} +zone={str(self.zone)}, +{self.pole} +datum={self.datum} +units={self.units} +no_defs"
        print(f"Building Projection with string '{proj_str}'")
        self.projector = Proj(proj_str)
        self.raster = np.array([0])
        self.units_per_pixel = units_per_pixel
    
    
def make_raster_from_geopackage(gr, geopackage_file, retain_geopackage = False):
    geo_data = gpd.read_file(geopackage_file)
    if(retain_geopackage):
        gr.geo_data = geo_data
    make_raster_from_geopandas(gr, geo_data)
    

def make_raster_from_geopandas(gr, geo_data):
    gr.geo_bounds = geo_data.total_bounds
    geo_width = gr.geo_bounds[2] - gr.geo_bounds[0]    
    geo_height = gr.geo_bounds[3] - gr.geo_bounds[1]
    
    meter_scale_x = gr.units_per_pixel
    meter_scale_y = gr.units_per_pixel
    
    if gr.pole == "north":
        meter_scale_y = -meter_scale_y
    
    transform = rasterio.Affine(meter_scale_x, 0.0, gr.geo_bounds[0], 
                                0.0, meter_scale_y, gr.geo_bounds[3])
    
    gr.raster = rasterize(
        [(shape, 1) for shape in geo_data['geometry']],
        out_shape=(int(geo_width/gr.units_per_pixel + 0.5),
                   int(geo_height/gr.units_per_pixel + 0.5)),
        transform=transform,
        fill=0,
        all_touched=True,
        dtype=rasterio.uint8)
    gr.raster = np.array(gr.raster,dtype=bool)   
    
def make_prefix_sum_table(gr,raster):
    dims = raster.shape
    d = np.zeros_like(raster,dtype=float)
    d[0][0] = raster[0][0]
    psum = raster[0][0]
    #first col
    for i in range(0,dims[0]):
        psum = psum + d[i][0]
        d[i][0] = psum;
    #first row    
    psum = raster[0][0]
    #first col
    for j in range(1,dims[1]):
        psum = psum + d[0][j]
        d[0][j] = psum;
    
    for i in range(1,dims[0]):
        for j in range(1,dims[0]):
            d[i][j] = d[i][j-1] + d[i][j-1] - d[i-1][j-1]
    
    gr.raster_prefix = d

def get_prefix_sum(gr,x1,x2, y1,y2):
    x1 = np.max(0,x1-1)
    y1 = np.max(y1-1,0)  
    y2 = np.min(y2,gr.raster_prefix.shape[1])
    x2 = np.min(x2,gr.raster_prefix.shape[0])
    region_sum = gr.raster_prefix[x2][y2] - gr.raster_prefix[x1][y2] - gr.raster_prefix[x2][y1] + gr.raster_prefix[x1][y1]
    return region_sum

def sample_raster(gr, i, j, box_radius):
    box_radius = int(box_radius)
    sample_array = gr.raster[(i-box_radius):(i+box_radius), 
                     (j-box_radius):(j+box_radius)]
    return sample_array
            
def sample_window(gr, lat, lon, radius):
    raster_coords = coord_to_pixel(gr, lat,lon)   
    #row, col --> y, x     
    return sample_raster(gr, raster_coords[1],raster_coords[0], radius/gr.units_per_pixel)

def sample_window_coverage(gr,lat,lon,radius):
    return np.mean(sample_window(gr, lat,lon,radius))


def plot_window(gr, lat, lon, radius):
    coverage_array = sample_window(gr, lat,lon,radius) * 255
    im = Image.fromarray(coverage_array)
    plt.imshow(im)
    print(np.mean(coverage_array))
    

def sample_window_coverage_prefix(gr,lat,lon,radius):
    im_coords = coord_to_pixel(gr, lat,lon) 
    half_box = int(radius/gr.units_per_pixel/2)
    x1 = im_coords[0]-half_box
    y1 = im_coords[1]-half_box
    x2 = im_coords[0] + half_box
    y2 = im_coords[1] + half_box
    
    return get_prefix_sum(x1,x2,y1,y2)


def sample_acres(gr, lat, lon, acres):            
    print(sample_sqft(gr,lat,lon,43560 * acres))

def sample_sqft(gr, lat, lon, sqft, variance_radius_ft=0):
    #box_size is the radius in <gr.units> of the box containing <sqft> square feet
    box_size = int(np.sqrt(sqft)/2.0)
    variance_margin = int(variance_radius_ft)
    
    if gr.units == "m":
        box_size = int(np.sqrt(sqft) / 3.28084 /2.0)        
        variance_margin = int(variance_radius_ft/3.28084)            
    
    coverage = 0
    variance = 0
    raster_coords = coord_to_pixel(gr, lat,lon) 
    samples = np.zeros(np.power(variance_margin * 2 + 1,2))
    ix = 0
    if variance_margin > 0:            
        for i in range(-variance_margin,variance_margin):
            for j in range(-variance_margin,variance_margin):        
                sample_array = sample_raster(gr, raster_coords[1]+i,raster_coords[0]+j,box_size)
                sample = np.mean(sample_array)
                
                samples[ix] = sample
                ix = ix+1
                if i == 0 and j == 0:
                    coverage = sample
            
        variance = np.var(samples)
    else:            
        coverage = sample_window_coverage(gr, lat, lon, int(box_size))        
    
    return coverage, variance
    
    
def coord_to_pixel(gr, lat,lon):
    easting, northing = gr.projector(lon,lat)    
    geo_width = gr.geo_bounds[2] - gr.geo_bounds[0]
    geo_height = gr.geo_bounds[3] - gr.geo_bounds[1]
    scale = gr.units_per_pixel
    pixel_x = int((easting - gr.geo_bounds[0])/scale + 0.5)
    # for image, 0 starts at top so need to flip by height
    pixel_y = int(geo_height/scale + 0.5) - int((northing - gr.geo_bounds[1])/scale + 0.5)
    return (pixel_x, pixel_y)


def save_raster(gr,file_path):    
    np.save(file_path,gr.raster)
    
def load_raster(gr, file_path):
    gr.raster = np.load(file_path)


def plot_geometry(gr, window_size):
    transform = rasterio.transform.from_bounds(*gr.geo_bounds, 
                                               *window_size)
    rasterized = rasterio.features.rasterize(
        [(shape, 1) for shape in gr.geo_data['geometry']],
        out_shape=window_size,
        transform=transform,
        fill=0,
        all_touched=True,
        dtype=rasterio.uint8)
    plt.imshow(rasterized)
