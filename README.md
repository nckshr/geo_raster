# geo_raster
This repository leverages the geopandas and rasterio python libraries to facilitate rendering binary rasters of geological shape data for quick sampling of area coverage based on lat/lon coordinates and a given radius.

Here's a small example using the [Minnesota National Wetland Inventory Data](https://gisdata.mn.gov/dataset/water-nat-wetlands-inv-2009-2014),
you can download the geopackage [here](https://resources.gisdata.mn.gov/pub/gdrs/data/pub/us_mn_state_dnr/water_nat_wetlands_inv_2009_2014/gpkg_water_nat_wetlands_inv_2009_2014.zip).

To Build the raster data (which is just a 2D binary numpy array):
```python
# Builds a raster of wetland from the Minnesota National Wetland Inventory Data
mn_graster = geo_raster.GeoRaster(units_per_pixel=10,
                 proj="utm", zone="15", north=True, datum="NAD83", units="m")

# resulting array is stored in mn_graster.raster
geo_raster.make_raster_from_geopackage(mn_graster, "water_nat_wetlands_inv_2009-2014.gpkg", retain_geopackage=False)
```

There are convenience methods for sampling around latitude / longitude coordiantes. For example,
here is a lake with an island:
<img src="https://i.imgur.com/Xw1xRQD.jpg" width=400)></img>

with the red marker having lat/lon coordinates of ```(46.792370, -93.219519)```.
We can view the region about this coordinate from the raster as follows:
```
# Radius is in the units specified when the GeoRaster object was constructed,
# in our case this is meters.
geo_raster.plot_window(mn_graster,46.792381, -93.219508,radius=500)
```
<img src="https://i.imgur.com/fT5szGw.png" width=400></img>
In this case, the white regions are 10m squares that are covered by a wetland region in the data.
We can see the island rendered as the black spot in the center of the window.

Other functions such as ```geo_raster.sample_window``` return the region as a numpy array
for general purpose use.
