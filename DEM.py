import re
import tarfile
import zipfile
import glob
import gdal
import ogr
import osr
import itertools

import numpy as np
import urllib.request

from os import path, remove, listdir, makedirs
from .NTS import nts


def get_tile_path_CDED(NTS):
    '''Get FTP path for a CDED NTS tile '''
    NTS = NTS.lower()
    
    # test resolution from NTS specification
    if len(NTS) == 6:
        resolution = "50k_dem" 
    elif len(NTS) == 4:
        resolution = "250k_dem" 
    else:
        raise Exception("Invalid NTS sheet!")
    
    # base path for all files
    basepath =  "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec"
    
    # build ftp path
    tile = NTS[0:3]
    ftp_path = "{}/{}/{}/{}".format(basepath, resolution, tile, NTS)
    ftp_path = ftp_path + ".zip"
    
    return(ftp_path)
    
def get_tile_path_NED(lon=None, lat=None,  name=None, test=True):
    ''' Get http path to download USGS NED elevation tile '''
    if not (lon or name):
        raise Exception("provide lat/lon or tile name")
        
    basepath = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/GridFloat/"
    
    if name is None:
        name = NED_tile_name(lon, lat, fext=".zip")
    else:
        name = name + ".zip"
        
    fullpath = basepath + name
    

    if test:
        try:
            code = urllib.request.urlopen(fullpath).getcode()
        
        except urllib.error.HTTPError:
            if lon is not None:
                name = NED_tile_name(lon, lat, fext=".zip", v="2017")
                fullpath = basepath + name
            else:
                fullpath = re.sub(name, "USGS_NED_1_" + name, fullpath)
                fullpath = re.sub("\\.zip", "_GridFloat.zip", fullpath)
        
    return(fullpath)
        
def NED_tile_name(lon, lat, fext="", v="2013"):
    ''' '''
    if v == "2013":
        name = "n{:02d}w{:03d}{}".format(np.abs(lat), np.abs(lon), fext)
    elif v == "2017":
        name = "usgs_ned_1_n{:02d}w{:03d}_gridfloat{}".format(np.abs(lat), np.abs(lon), fext)
 
    return(name)

def get_tile_path_CDEM():
    raise NotImplementedError()
    
def download_single_DEM(DEM_id, DEM_dir, replace=False, product="NED"):
    '''
    
    '''
    output = True
    
    if product.upper() == "NED":
        ftp_path = get_tile_path_NED(name = DEM_id)
    elif product.upper() == "CDED":
        ftp_path = get_tile_path_CDED(NTS = DEM_id)
    else:
        raise NotImplementedError("DEM product not implemented")
    
    # create appropriate file name / directory structure based on ftp path
    file_paths = ftp_path.split('/')[-3:]
    
    save_dir = path.join(DEM_dir, *file_paths[0:2])
    if not path.isdir(save_dir):
        makedirs(save_dir)
    
    destfile = path.join(DEM_dir, *file_paths)
    
    if destfile.endswith("tar.gz"):
        dest_dir = re.sub("\\.tar\\.gz", "", destfile)
    else:
        dest_dir = re.sub("\\.zip", "", destfile)

    # Check to see if file already exists
    if not replace and path.isdir(dest_dir):
        print("{} exists locally and was not downloaded\n".format(dest_dir))
    else:
        output = download_and_unzip(url = ftp_path, destfile = destfile, exdir = dest_dir)
        
    # If an appropriate file was downloaded, return the corresponding file paths
    if output:
        pattern_dict = {"CDED" : "dem[ew_].*[td][ie][fm]$",
                        "CDEM" : "dem[ew_].*[td][ie][fm]$",
                        "CDSM" : ".*_cdsm_final_[ew]\\.tif",
                        "NED"  : "flt$",
                        "ADEM" : "reg_dem\\.tif$"}
    
        pattern = pattern_dict[product]

        dem = [f for f in listdir(dest_dir) if re.search(pattern, f)]
        dem = [path.join(dest_dir, x) for x in dem]
        
        return(dem)
    
def download_and_unzip(url, destfile, exdir, rmzip=True):
    ''' 
    Downloads and unzips a file

    @args
        url:       character url path
        destfile:  character, filepath of output zipfile
        exdir:     character,  the directory to which files are extracted
        rmzip:     logical, whether or not to remove zipfile after extraction. 
        
    @return
        path(s) to target tiles
    '''
    # try:
    # download file
    try:
        print("Downloading file from {}".format(url))
        urllib.request.urlretrieve(url, destfile)
   
    # if the url doesn't exist    
    except Exception as e:
        print(url)
        print(e)
        return(False)
    
    if destfile.endswith("tar.gz"):
        with tarfile.open(destfile, "r:gz") as tarf:
            tarf.extractall(exdir)
        
    elif destfile.endswith("zip"):
        with zipfile.ZipFile(destfile, "r") as zipf:
            zipf.extractall(exdir)
            
    else:
        rmzip = False
        
    if rmzip:
        remove(destfile)
    
    return(True)
        


def download_multiple_DEM(DEM, DEM_dir, product="NED"):
    ''' Download a list of DEM URLs. If they exist already, they are not downloaded
    @args
        DEM:      list of DEM urls (NED) or NTS tiles (CDED)
        DEM_dir:  path to which files are downloaded
        product:  which DEM tile series should be downloaded: ('NED', 'CDED')
        
    @return 
        a list of file paths for target DEMs
    '''
    #sanity check - make sure NTS names are well-formed
    if product.upper() in ["CDED", "CDEM", "CDSM"]:
        if not all([re.search("^\\d{3}\\w(\\d{2})?$", x) for x in DEM]):
            raise Exception("Bad format for one or more NTS strings")
            
    
    # download each DEM file using the map function
    get_single = lambda x: download_single_DEM(x, DEM_dir = DEM_dir, product=product)        
    files = map(get_single, DEM)
    
    # return list of files
    files = [f for f in files if f is not None]
    files = [dem for sublist in files for dem in sublist]
    return(files)
        
def create_DEM_mosaic(DEM, DEM_dir, dstfile, product="NED", ellipsoidal=False,
                        vrt_only=False, format="GTiff"):
    ''' Create a Mosaic from a list of DEM urls or NTS tiles. Missing tiles will
    be downloaded'''
    
    
    files = download_multiple_DEM(DEM, DEM_dir, product)
    
    # build VRT
    VRT_path = path.join(path.dirname(dstfile), "tmp.VRT")
    VRT = gdal.BuildVRT(VRT_path, files)
    VRT.FlushCache()
    VRT = None
    
    # return VRT-only if desired
    if vrt_only:
        if ellipsoidal:
            raise Exception("Cannot convert VRT to ellipsoidal heights")
        return(VRT_path)
    
    # set warp parameters
    if ellipsoidal:
        wo = gdal.WarpOptions(srcSRS=DEMproj4(product),  dstSRS='epsg: 4623',
                                format=format)
        ds = gdal.Warp(dstfile,  VRT_path, options=wo)
    else:
        ds = gdal.Translate(dstfile, VRT_path, format=format)

    ds.FlushCache()
    ds = None
    
    remove(VRT_path)
    return(dstfile)
    

def NED_tiles_from_extent(ext):
    ''' get list of NED tiles required to cover a spatial extent'''
    # unpack extent dictionary
    xmin = ext['xmin']
    xmax = ext['xmax']
    ymin = ext['ymin']
    ymax = ext['ymax']
    
    # get all corner coordinates to cover extent
            # +1 beacuse of 0 indexing 
    xrange = range(int(np.floor(xmin)), int(np.floor(xmax)) +1) 
    yrange = range(int(np.ceil(ymin)), int(np.ceil(ymax)) +1)
    pts = itertools.product(xrange, yrange)
  
    # get tile index for all 
    f = lambda x: NED_tile_name(lon = x[0], lat = x[1])
    tiles =  [pth for pth in map(f, pts)]
    
    return(tiles)
    
def NTS_tiles_from_extent(ext, scale=1):
    ''' Determine which NTS tiles are required to cover a target spatial extent 
    
    @usage:
        ext = {'xmin': 52, 'xmax': 53, 'ymin' : -114, 'ymax' : -112}
        NTS_tiles_from_extent(ext)
    '''
    # unpack extent dictionary
    w = ext['xmin']
    e = ext['xmax']
    s = ext['ymin']
    n = ext['ymax']
    
    # find NTS tiles
    bbox = nts.makebbox(n=n, e=e, s=s, w=w)
    tiles = nts.bybbox(bbox, scale)
    
    # convert to list of strings
    tile_list = [''.join(tile) for tile in tiles]
    
    return(tile_list)

def get_spatial_extent(raster_path, target_EPSG = 4326, tol=0.1):
    ''' Get spatial extent for raster file. If not georeferenced, attempts
    to use the GCPs in the image (e.g. for radarsat 2). Howver, this doesn't
    produce exact results, so it is advisable to use an extra buffer tolerance
    in your spatial extent (maybe ~0.1 decimal degrees)'''

    # open file
    src = gdal.Open(raster_path)
    prj = src.GetProjection()
    
    # check if georeferencing information is available
    if (sum(src.GetGeoTransform()) == 2) and src.GetGCPCount() > 3:
        gt = gdal.GCPsToGeoTransform(src.GetGCPs())
        prj = src.GetGCPProjection()
    else:
        gt = src.GetGeoTransform()
    
    # get untransformed upper left and lower right corner coordinates
    ulx, xres, xskew, uly, yskew, yres  = gt
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    
    
    # Setup the source projection - you can also import from epsg, proj4...
    source = osr.SpatialReference()
    source.ImportFromWkt(prj)

    # The target projection
    target = osr.SpatialReference()
    target.ImportFromEPSG(target_EPSG)
    
    # set up transform
    transform = osr.CoordinateTransformation(source, target)

    # transform coordinates
    upper = transform.TransformPoint(ulx, uly)
    lower = transform.TransformPoint(lrx, lry)
    
    ext = {'xmin' : min(upper[0], lower[0]) - tol,
           'xmax' : max(upper[0], lower[0]) + tol,
           'ymin' : min(upper[1], lower[1]) - tol,
           'ymax' : max(upper[1], lower[1]) + tol}
    

    return(ext)
    
def create_DEM_mosaic_from_extent(ext, dstfile, DEM_dir, product="CDED", 
                                    vrt_only=False, ellipsoidal=False):
    ''' 
    generate DEM mosaic covering extent
    
    @args
        ext : dictionary with keys xmin, xmax, ymin, ymax  and values in decimal degrees
        product : one of "NED", "CDED".  DEM source
    
    @example
        from os import path
        home = path.expanduser('~')
        ext = {'ymin': 52, 'ymax': 53, 'xmin' : -114, 'xmax' : -112}
        create_DEM_mosaic_from_extent(ext, 
                                      dstfile = path.join(home, 'mosaic.tif'),
                                      DEM_dir = path.join(home, 'DEM'),
                                      product = "CDED",
                                      vrt_only = False,
                                      ellipsoidal = True)
    ''' 
    
    # 
    if product.upper() == "NED":
        tiles = NED_tiles_from_extent(ext)
    elif product.upper() == "CDED":
        tiles = NTS_tiles_from_extent(ext)
    else:
        raise NotImplementedError
    
    # create mosaic
    mosaic = create_DEM_mosaic(DEM = tiles, DEM_dir=DEM_dir, dstfile=dstfile, 
                                product=product, vrt_only=vrt_only, 
                                ellipsoidal=ellipsoidal)
    
    return(mosaic)
  
def DEMproj4(product):
    ''' return proj4 string or equivalent for DEM product'''
    P = {"CDED": "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +geoidgrids=HT2_0.gtx +no_defs",
         "NED": "epsg: 4269 + 5703",
         "CDEM": "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +geoidgrids=HT2_0.gtx +no_defs"
         }
    return(P[product.upper()])
    

    
def gdalslope(DEM, dst, latlon = True):
    if latlon:
        scale = 111120
    else:
        scale  = 1
    gdal.DEMProcessing(dst, DEM, "slope", scale=scale) 
    
def gdalTPI(DEM, dst, latlon = True):
    gdal.DEMProcessing(dst, DEM, "TPI") 

# gdalslope(r"C:\NB\DEM_m.tif",  r"C:\NB\DEM_m_gdalslope.tif")
    
def TWI(DEM):
    pass
    
    
