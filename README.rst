=====
DEMpy
=====

DEMPY makes it easy to download and mosaic DEM tiles from multiple providers based on areal extents or coordinates.


Installation
============

From the package
^^^^^^^^^^^^^^^^

To install DEMpy, download the package and navigate to the installation directory then install using setuptools: 

.. code-block:: bash

    cd DEMpy
    python setup.py install

If you are interested in the DEM mosaicking functionality, you'll need the GDAL libraries and python wrappers. These can be difficult to set up on Windows, but I suggest using the terrific resources from `<http://www.gisinternals.com/archive.php>`_ and  `<https://www.lfd.uci.edu/~gohlke/pythonlibs/>`_ 

Using docker
^^^^^^^^^^^^
Alternately, you can use docker to run the scripts and avoid having to compile GDAL. This is only important if you want to run the mosaicking functions. To do so, run the following from within the main DEMpy directory.


.. code-block:: bash

    docker build . -f docker/Dockerfile -t dempy
    docker run -it dempy


Examples
========

.. code-block:: python

    from dempy.DEM import *
    from os import path

    DEM_directory = path.join(path.expanduser("~"), "dempy")

    ''' Find URL of DEM tiles '''
    # USGS NED (30m, north america)
    print(get_tile_path_NED(lon=-114, lat=52))

    # CDED
    print(get_tile_path_CDED('075C'))
    print(get_tile_path_CDED('075C02'))


    ''' download single DEM tiles by ID ''' 


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
    
