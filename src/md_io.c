//#include <netcdf.h>
#include <string>
#include <cstring>
#include <netcdf.h>
#include <stdio.h>
#include <mpi.h>
#include "md_utils.h"

#define BLEN 200

int total_written = 1;

// Read a line and then return the first string with whitespace stripped off.
int get_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;

        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}

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


// Handles output printing.
void output(mdsys_t *sys, FILE *erg, FILE *traj)
{
  printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
         sys->ekin, sys->epot, sys->ekin + sys->epot);
  fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
          sys->ekin, sys->epot, sys->ekin + sys->epot);

  fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
  for (int i=0; i < sys->all_atoms; ++i) {
      fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
  }

  // Print trajectory using NetCDF format.
  //write_to_netcdf(&sys->coordinates[0], &sys->coordinates[sys->all_atoms],
  //                &sys->coordinates[2*sys->all_atoms], sys->all_atoms, trajfile);
}

void write_output(mdsys_t *sys, FILE *erg, FILE *traj) {
  if (sys->my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, (void*) &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, (void*) &sys->ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    output(sys, erg, traj);
  } else {
    MPI_Reduce((void*) &sys->epot, (void*) &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce((void*) &sys->ekin, (void*) &sys->ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }
}

void write_initial_output(mdsys_t *sys, FILE *erg, FILE *traj) {
  if (sys->my_rank == 0) {

    printf("Starting simulation with %d atoms for %d steps.\n",sys->natoms, sys->nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    MPI_Reduce(MPI_IN_PLACE, (void*) &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, (void*) &sys->ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    output(sys, erg, traj);
  } else {
    MPI_Reduce((void*) &sys->epot, (void*) &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce((void*) &sys->ekin, (void*) &sys->ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }
}
