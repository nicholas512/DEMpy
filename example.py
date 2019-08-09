from DEMpy.DEM import *


''' Find URL of DEM tiles '''
# USGS NED (30m, north america)
print(get_tile_path_NED(lon=-114, lat=52))

# CDED
print(get_tile_path_CDED('075C'))
print(get_tile_path_CDED('075C02'))


''' download single DEM tiles by ID ''' 
DEM_directory = "C:/DEM/"

# CDED
download_single_DEM('075C02', DEM_directory, product="CDED")

# NED
NED_id = NED_tile_name(lon=-114, lat=52)
download_single_DEM(NED_id, DEM_directory, product="NED")

# Won't re-download them once they exist
download_single_DEM('075C02', DEM_directory, product="CDED")


''' download many DEM tiles by extent ''' 
ext = {'ymin':53, 'ymax': 54, 'xmin' : -120, 'xmax' : -117}

# use scale parameter to choose either 250k or 50k
print(NTS_tiles_from_extent(ext, scale=1)) # 250k
print(NTS_tiles_from_extent(ext, scale=2)) # 50k

create_DEM_mosaic_from_extent(ext=ext, dstfile="C:/DEM/DEM_mosaic.tif", DEM_dir=DEM_directory, product='CDED', scale=1)