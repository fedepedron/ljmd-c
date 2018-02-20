#include <netcdf.h>
#include <string>
void write_to_netcdf(const double *X, const double *Y, const double *Z, int N,
                     const std::string &name);
