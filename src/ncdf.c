//#include <netcdf.h>
#include <string>
#include <cstring>
#include "ncdf.h"

int total_written = 1;

void write_to_netcdf(const double *X, const double *Y, const double *Z, int N,
                     const std::string &name)
{
  const int NDIMS = 1;
  int ncid, dimids[NDIMS], varid_x, varid_y, varid_z;

  std::string fname = name + std::to_string(total_written) + ".nc";
  nc_create(fname.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid);

  /* Define the dimensions. NetCDF will hand back an ID for each. */
  nc_def_dim(ncid, "Atom", N, &dimids[0]);
  nc_def_var(ncid, "X", NC_DOUBLE, NDIMS, dimids, &varid_x);
  nc_def_var(ncid, "Y", NC_DOUBLE, NDIMS, dimids, &varid_y);
  nc_def_var(ncid, "Z", NC_DOUBLE, NDIMS, dimids, &varid_z);

  /* set COMPRESSION!!!! This works better for non-contiguous data*/
  int shuffle = 0, deflate = 1, deflate_level = 1;
  nc_def_var_deflate(ncid, varid_x, shuffle, deflate, deflate_level);
  nc_def_var_deflate(ncid, varid_y, shuffle, deflate, deflate_level);
  nc_def_var_deflate(ncid, varid_z, shuffle, deflate, deflate_level);
  nc_enddef(ncid);

  nc_put_var_double(ncid, varid_x, &X[0]); // write all data
  nc_put_var_double(ncid, varid_y, &Y[0]);
  nc_put_var_double(ncid, varid_z, &Z[0]);
  nc_close(ncid);
  total_written++;
}
