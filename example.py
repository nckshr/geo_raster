
import geo_raster
import json 
import glob
import numpy as np

#Example usage that annotates real estate listings with estimated wetland coverage from the 2009-2014 MN wetland geological survey data

#sets up values for coordinate mapping
mn_graster = geo_raster.GeoRaster(units_per_pixel=10,
                 proj="utm", zone="15", north=True, datum="NAD83", units="m")

# as the geopackage can be quite large, you might not want it to hang around in memory
# this will consider all shapes in the geopackage as "covered" areas.
# to filter shapes, load the geopackage first with geopandas, filter and use 
# geo_raster.make_raster_from_geopanas(geopandas_data_frame)

geo_raster.make_raster_from_geopackage(mn_graster, "water_nat_wetlands_inv_2009-2014.gpkg", retain_geopackage=False)

#saving the raster can save time & space later (wraps numpy.save)
geo_raster.save_raster(mn_graster, "mn_wetland_10m")
# to load:
# geo_raster.load_raster(mn_graster, "mn_wetland_10m.npy")

#annotate real estate listings with wetland coverage (values in 0-1)
listing_file = "mn_listings.json"
with  open(listing_file) as listing_file_obj:
listing_data = json.load(listing_file_obj)

for lx in range(0,len(listing_data)):
    listing = listing_data[lx]
    sqft = listing['description']['lot_sqft']
    wetness = 0
    if sqft is not None and sqft > 100:
        if 'location' in listing and 'address' in listing['location'] and 'coordinate' in listing['location']['address']:
            # returns proportion of sqrt(sqft) x sqrt(sqft) window centered on listing coordinate
            # that is covered by wetland
            wetness = geo_raster.sample_sqft(mn_graster,
                                        float(listing['location']['address']['coordinate']['lat']),
                                        float(listing['location']['address']['coordinate']['lon']),
                                        int(sqft))
            if not np.isnan(wetness[0]):
                listing['description']['wetness'] = wetness[0]
            listing_data[lx] = listing
print(f"Writing annotated listings for {listing_file}")
with open(f"{listing_file}_annotated.json", "w") as file_out:
    json.dump(listing_data,file_out)  
