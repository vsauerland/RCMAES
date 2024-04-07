import numpy as np
import netCDF4 as nc
from andacc import load
from tmm_mops_andacc_utils import tmmmops2xy, readFromPETScRestart
from pytmmutils import matrixToGrid

# NOTE by IK: This create will create a netCDF-4 file named tracers.nc, that uses HDF as underlying storage layer.
# It may create problems when reading into older ferret version. If you don't want that, convert the file with 
# nccopy -3 tracers.nc <outfile> 

def petsc2netcdf( ncFile ):
  """ Read all single tracer PETSCs files and write to nc file """

  data = load( "data_for_tmm_mops" )
  data.trFileNamesEnd = [ "po4end.petsc", "dopend.petsc", "oxyend.petsc", "phyend.petsc", "zooend.petsc", "detend.petsc", "no3end.petsc", "dicend.petsc", "alkend.petsc", "runoffend.bin" ]
  tmparr = readFromPETScRestart( data.trNames[ 0 : 10 ], data.trFileNamesEnd[ 0 : 10 ] )
  x, y, Cb = tmmmops2xy( tmparr, data )
  nb = Cb.shape[ 0 ] # number of wet boxes
  ntr = Cb.shape[ 1 ] # number of tracers

  # According to the Matlab code of Samar Khatiwala and assuming that
  # "matrixToGrid.py" is the python version of "grid_boxes3d.m",
  # Ib has to be calculated as follows
  # (Ib tells us how to rearrange the tracer vectors stored in Cb)
  b = load( "boxes" )
  Ibs = np.where( b.izBox == 1 )[ 0 ]  # indices of all ocean surface grid boxes
  nbs = len( Ibs )  # number of ocean surface grid boxes (water columns)
  Ips = {}  # dictionary to store grid box indices of each water column
  for i in range( nbs ):
    ibs = Ibs[ i ]  # current ocean surface grid box index
    Ips[ i ] = np.where( ( b.Xboxnom == b.Xboxnom[ ibs ] ) & ( b.Yboxnom == b.Yboxnom[ ibs ] ) )[ 0 ]
    # Ips{ i } contains all grid box indices of the current water column 
    # note, that the list may not be ordered w.r.t. the depth levels of the boxes
    izps = np.argsort( b.Zboxnom[ Ips[ i ] ] )
    # izps tells how to rearrange the entries of Ips[ i ] to be orderd by depth
    Ips[ i ] = Ips[ i ][ izps ]
  Ibi = np.concatenate( list( Ips.values() ) )
  # Ibi is the index vector that stacks all water columns on top of each other
  Ib = np.argsort( Ibi )  # Ib tells us how to rearrange any tracer vector
  # Now, sort all columns of Cb (tracer vectors) acc. to Ib
  for i in range( ntr ):
    Cb[ :, i ] = Cb[ Ib, i ]

  # Apply matrixToGrid which does the same as Matlab "grid_boxes3d.m"
  Cg, xg, yg, zg = matrixToGrid( Cb, np.arange( nb ), 'boxes', 'grid', keepEmptyLayers = 0, emptyLayerFillValue = 1 )

  # Write all 3D tracer variables as given by Cg to ncFile
  with nc.Dataset( ncFile, 'w' ) as ncid:
    # Define dimensions
    lon_dimid = ncid.createDimension( 'Longitude', len( xg ) )
    lat_dimid = ncid.createDimension( 'Latitude', len( yg ) )
    dep_dimid = ncid.createDimension( 'Depth', len( zg ) )

    # Define coordinate variables
    lon_varid = ncid.createVariable( 'Longitude', 'f8', ( 'Longitude', ) )
    lat_varid = ncid.createVariable( 'Latitude', 'f8', ( 'Latitude', ) )
    dep_varid = ncid.createVariable( 'Depth', 'f8', ( 'Depth', ) )

    # Assign values to coordinate variables
    lon_varid[ : ] = xg
    lat_varid[ : ] = yg
    dep_varid[ : ] = zg

    for i in range( 9 ):
      # Define i-th tracer variable and add missing_value attribute NaN
      data_varid = ncid.createVariable( data.trNames[ i ], 'f8', ( 'Depth', 'Latitude', 'Longitude' ) )
      data_varid.setncattr( 'missing_value', np.nan ) 

      # Assign i-th tracer variable
      # (transposition of arrays before writing to NetCDF
      # appears to be required by "NetCDF order of dimension convention")
      data_varid[ :, :, : ] = np.transpose( Cg[ :, :, :, i ], (2, 1, 0) )

  # The return fields are for comparison with Matlab version:
  # (Cb corresponds to the vi's in n7tracers.m, Ib to Irr, and Cg to the Vi's)
  print( "Written 3d tracers to NetCDF file", ncFile ) 
  return( Cb, Ib, Cg )
