#include <stdlib.h>
#include "md_utils.h"

// These zero arrays of int and double.
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

void azzero_i(int *a, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        a[i]=0.0;
    }
}

// Apply minimum image convention
double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

// Sets up the arrays for MPI_AllgatherV.
void ALLGV_setup(mdsys_t *sys, allgv_t *allgv_data) {

  allgv_data->sendbuf   = (double *)malloc(sys->my_atoms *sizeof(double)*3);
  allgv_data->recvbuf   = (double *)malloc(sys->all_atoms*sizeof(double)*3);
  allgv_data->recvcounts=    (int *)malloc(sys->n_tasks  *sizeof(int));
  allgv_data->displs    =    (int *)malloc(sys->n_tasks  *sizeof(int));

  allgv_data->displs[0] = 0;
  if (sys->task_rest > 0) {
    int accum_count = 0;
    for (int t = 0; t < sys->task_rest; t++) {
      allgv_data->recvcounts[t] = 3*sys->atom_count_L;
      accum_count              += allgv_data->recvcounts[t];
      allgv_data->displs[t+1]   = accum_count;
    }
    for (int t = sys->task_rest; t < sys->n_tasks -1; t++) {
      allgv_data->recvcounts[t] = 3*sys->atom_count_S;
      accum_count              += allgv_data->recvcounts[t];
      allgv_data->displs[t+1]   = accum_count;
    }
  } else {
    int accum_count = 0;
    for (int t = 0; t < sys->n_tasks -1; t++) {
      allgv_data->recvcounts[t] = 3*sys->atom_count_S;
      accum_count              += allgv_data->recvcounts[t];
      allgv_data->displs[t+1]   = accum_count;
    }
  }
  allgv_data->recvcounts[sys->n_tasks - 1] = 3*sys->atom_count_S;
}
