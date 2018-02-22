/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "timers.h"
#include "md_utils.h"
#include "md_io.h"

/* generic file- or pathname buffer length */
#define BLEN 200
#define NGB_MAX 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */


//##############################################################################
// HELPER FUNCTIONS


//##############################################################################
// MAIN SUBS
/* compute kinetic energy */
static void create_atom_list(mdsys_t *sys)
{
  int offset;
  // Sets parameters for array usage.
  sys->atom_count_S  = (int) (sys->natoms / sys->n_tasks );
  sys->atom_count_L  = sys->atom_count_S + 1;
  sys->task_rest     = (int) (sys->natoms % sys->n_tasks );

  sys->my_atoms = sys->atom_count_S;
  if (sys->my_rank < sys->task_rest) sys->my_atoms++;

  // Initializes local atom variables.
  sys->my_atom_list = (int*)malloc((sys->my_atoms)*sizeof(int));
  if (sys->task_rest > 0) {
    if (sys->my_rank < sys->task_rest) {
      offset = sys->atom_count_L * sys->my_rank;
    } else {
      offset = sys->atom_count_L * sys->task_rest +
              (sys->my_rank - sys->task_rest) * sys->atom_count_S;
    }
  } else {
    offset = sys->my_rank * sys->my_atoms;
  }
  for (int i = 0; i < sys->my_atoms ; i++) {
    sys->my_atom_list[i] = i + offset;
  }

  // Initializes local ghost atom variables.
  sys->ghost_atoms  = sys->natoms - sys->my_atoms;
  sys->ghost_atom_list = (int*) malloc((sys->ghost_atoms)*sizeof(int));
  for (int i = 0; i < offset ; i++) {
    sys->ghost_atom_list[i] = i;
  }
  for (int i = offset; i < sys->ghost_atoms ; i++) {
    sys->ghost_atom_list[i] = i + sys->my_atoms;
  }

  sys->all_atoms = sys->ghost_atoms + sys->my_atoms;
}

static void ekin(mdsys_t *sys)
{
  sys->ekin=0.0;
  for (int id=0; id<sys->my_atoms; ++id) {
    int i = sys->my_atom_list[id];
    sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] +
                 sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
  }
  sys->temp = 2.0*sys->ekin/(3.0*sys->my_atoms-3.0)/kboltz;
}

/* compute forces */
static void force_ngb(mdsys_t *sys)
{
  sys->epot=0.0;
  azzero(sys->fx,sys->my_atoms);
  azzero(sys->fy,sys->my_atoms);
  azzero(sys->fz,sys->my_atoms);

  double epot = 0.0;

  for(int iloc=0; iloc < (sys->my_atoms); iloc++) {
    int i = sys->my_atom_list[iloc];
    // Iterates through all local atoms.
    for(int ingb=0; ingb < (sys->ngb_loc[iloc]); ingb++) {
      int jloc = sys->ngb_list[ingb + iloc*NGB_MAX];
      int j    = sys->my_atom_list[jloc];

      if (iloc > jloc) continue;
      double rx    = pbc(sys->coordinates[i] -
                         sys->coordinates[j], 0.5*sys->box);
      double ry    = pbc(sys->coordinates[i + sys->all_atoms] -
                         sys->coordinates[j + sys->all_atoms],
                         0.5*sys->box);
      double rz    = pbc(sys->coordinates[i + 2*(sys->all_atoms)] -
                         sys->coordinates[j + 2*(sys->all_atoms)],
                         0.5*sys->box);
      double r  = sqrt(rx*rx + ry*ry + rz*rz);

      double ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r +
                                           6*pow(sys->sigma/r,6.0)/r);
      epot += 4.0*sys->epsilon*(pow(sys->sigma/r,12.0) -
                                pow(sys->sigma/r,6.0));
      sys->fx[iloc] += rx/r*ffac;
      sys->fy[iloc] += ry/r*ffac;
      sys->fz[iloc] += rz/r*ffac;
      sys->fx[jloc] -= rx/r*ffac;
      sys->fy[jloc] -= ry/r*ffac;
      sys->fz[jloc] -= rz/r*ffac;
    }
    // Iterates through all ghost atoms.
    for(int ingb = sys->ngb_loc[iloc]; ingb < sys->n_ngb[iloc]; ingb++) {
      int jloc = sys->ngb_list[ingb + iloc*NGB_MAX];
      int j    = sys->ghost_atom_list[jloc];
      double rx    = pbc(sys->coordinates[i] -
                         sys->coordinates[j], 0.5*sys->box);
      double ry    = pbc(sys->coordinates[i + sys->all_atoms] -
                         sys->coordinates[j + sys->all_atoms],
                         0.5*sys->box);
      double rz    = pbc(sys->coordinates[i + 2*(sys->all_atoms)] -
                         sys->coordinates[j + 2*(sys->all_atoms)],
                         0.5*sys->box);
      double r = sqrt(rx*rx + ry*ry + rz*rz);

      double ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r +
                                           6*pow(sys->sigma/r,6.0)/r);
      epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0) -
                                    pow(sys->sigma/r,6.0));
      sys->fx[iloc] += rx/r*ffac;
      sys->fy[iloc] += ry/r*ffac;
      sys->fz[iloc] += rz/r*ffac;
    }
  }

  sys->epot = epot;
  return;
}

/* Create Neighbour list using cutoff */
void make_ngb_list(mdsys_t *sys)
{
  azzero_i(sys->ngb_loc, sys->my_atoms);
  azzero_i(sys->n_ngb  , sys->my_atoms);
  for(int i=0; i < sys->my_atoms; i++) {
    int my_id = sys->my_atom_list[i];
    for(int j = i+1; j < sys->my_atoms; j++) {

      int atom_id  = sys->my_atom_list[j];
      double rx    = pbc(sys->coordinates[my_id] -
                         sys->coordinates[atom_id], 0.5*sys->box);
      double ry    = pbc(sys->coordinates[my_id + sys->all_atoms] -
                         sys->coordinates[atom_id + sys->all_atoms],
                         0.5*sys->box);
      double rz    = pbc(sys->coordinates[my_id + 2*(sys->all_atoms)] -
                         sys->coordinates[atom_id + 2*(sys->all_atoms)],
                         0.5*sys->box);
      double r     = rx*rx + ry*ry + rz*rz;
      if (r < (sys->rcut * sys->rcut)) {
        sys->ngb_list[i*NGB_MAX + sys->n_ngb[i]] = j;
        sys->ngb_list[j*NGB_MAX + sys->n_ngb[j]] = i;
        sys->n_ngb[i]++;
        sys->n_ngb[j]++;
        sys->ngb_loc[i]++;
        sys->ngb_loc[j]++;
      }
    }

    for(int j=0; j < sys->ghost_atoms; j++) {
      int atom_id  = sys->ghost_atom_list[j];
      double rx    = pbc(sys->coordinates[my_id] -
                         sys->coordinates[atom_id], 0.5*sys->box);
      double ry    = pbc(sys->coordinates[my_id + sys->all_atoms] -
                         sys->coordinates[atom_id + sys->all_atoms],
                         0.5*sys->box);
      double rz    = pbc(sys->coordinates[my_id + 2*(sys->all_atoms)] -
                         sys->coordinates[atom_id + 2*(sys->all_atoms)],
                         0.5*sys->box);
      double r     = rx*rx + ry*ry + rz*rz;

      if (r < (sys->rcut * sys->rcut)) {
        sys->ngb_list[i*NGB_MAX + sys->n_ngb[i]] = j;
        sys->n_ngb[i]++;
      }
    }
  }
  return;
}

static void update_coords(mdsys_t *sys, allgv_t *allgv_data)
{
  int i;
  for (int id = 0; id < sys->my_atoms; id++) {
    i = sys->my_atom_list[id];
    allgv_data->sendbuf[id] = sys->coordinates[i];
  }
  for (int id = 0; id < sys->my_atoms; id++) {
    i = sys->my_atom_list[id];
    allgv_data->sendbuf[id + sys->my_atoms] =
                                      sys->coordinates[i + sys->all_atoms];
  }
  for (int id = 0; id < sys->my_atoms; id++) {
    i = sys->my_atom_list[id];
    allgv_data->sendbuf[id + (sys->my_atoms)*2] =
                                      sys->coordinates[i + (sys->all_atoms)*2];
  }

  MPI_Allgatherv((void*)allgv_data->sendbuf, sys->my_atoms*3, MPI_DOUBLE,
                 (void*)allgv_data->recvbuf, allgv_data->recvcounts,
                 allgv_data->displs, MPI_DOUBLE, MPI_COMM_WORLD);
  // Copies coordinates from coordinate buffer.
  if (sys->task_rest > 0) {
    int offset = 0;
    int index_a = 0, index_b = 0;
    for (int t = 0; t < sys->task_rest; t++) {
      for (int i = 0; i < sys->atom_count_L; i++) {
        index_a = i + t*sys->atom_count_L;
        index_b = i + offset;
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      for (int i = 0; i < sys->atom_count_L; i++) {
        index_a = i + t*sys->atom_count_L + sys->all_atoms;
        index_b = i + offset + sys->atom_count_L;
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      for (int i = 0; i < sys->atom_count_L; i++) {
        index_a = i + t*sys->atom_count_L + 2*sys->all_atoms;
        index_b = i + offset + 2*(sys->atom_count_L);
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      offset += 3*sys->atom_count_L;
    }
    for (int t = sys->task_rest; t < sys->n_tasks; t++) {
      for (int i = 0; i < sys->atom_count_S; i++) {
        index_a = i + sys->task_rest + t*sys->atom_count_S;
        index_b = i + offset;
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      for (int i = 0; i < sys->atom_count_S; i++) {
        index_a = i + sys->task_rest + t*sys->atom_count_S + sys->all_atoms;
        index_b = i + offset + sys->atom_count_S;
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      for (int i = 0; i < sys->atom_count_S; i++) {
        index_a = i + sys->task_rest + t*sys->atom_count_S + 2*sys->all_atoms;
        index_b = i + offset + 2*(sys->atom_count_S);
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      offset += 3*sys->atom_count_S;
    }
  } else {
    int offset = 0;
    int index_a = 0, index_b = 0;
    for (int t = 0; t < sys->n_tasks; t++) {
      for (int i = 0; i < sys->atom_count_S; i++) {
        index_a = i + t*sys->atom_count_S;
        index_b = i + offset;
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      for (int i = 0; i < sys->atom_count_S; i++) {
        index_a = i + t*sys->atom_count_S + sys->all_atoms;
        index_b = i + offset + sys->atom_count_S;
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      for (int i = 0; i < sys->atom_count_S; i++) {
        index_a = i + t*sys->atom_count_S + 2*sys->all_atoms;
        index_b = i + offset + 2*(sys->atom_count_S);
        sys->coordinates[index_a] = allgv_data->recvbuf[index_b];
      }
      offset += 3*sys->atom_count_S;
    }
  }
}

/* velocity verlet */
static void velverlet(mdsys_t *sys, allgv_t *allgv_data)
{
  /* first part: propagate velocities by half and positions by full step */
  for (int id=0; id<sys->my_atoms; ++id) {
    int i = sys->my_atom_list[id];
    sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[id] / sys->mass;
    sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[id] / sys->mass;
    sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[id] / sys->mass;
    sys->coordinates[i]                      += sys->dt*sys->vx[i];
    sys->coordinates[i + sys->all_atoms]     += sys->dt*sys->vy[i];
    sys->coordinates[i + 2*(sys->all_atoms)] += sys->dt*sys->vz[i];
  }

  if (sys->my_rank == 0)  timer_start("Coordinate Messaging");
  update_coords(sys, allgv_data);

  if (sys->my_rank == 0)  {
    timer_pause("Coordinate Messaging");
    timer_start("Neighbour List");
  }
  if ((sys->nfi % sys->ngb_freq) == 0 ) make_ngb_list(sys);
  if (sys->my_rank == 0) {
    timer_pause("Neighbour List");
    timer_start("Force Calculation");
  }
  /* compute forces and potential energy */
  force_ngb(sys);
  if (sys->my_rank == 0) timer_pause("Force Calculation");

  /* second part: propagate velocities by another half step */
  for (int id=0; id<sys->my_atoms; ++id) {
    int i = sys->my_atom_list[id];
    sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[id] / sys->mass;
    sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[id] / sys->mass;
    sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[id] / sys->mass;
  }
}

//##############################################################################
int main(int argc, char **argv)
{
  MPI_Init(NULL,NULL);

  char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
  FILE *fp,*traj,*erg;
  mdsys_t sys;
  mdsys_inp_t input_buffer;

  MPI_Comm_rank(MPI_COMM_WORLD, &sys.my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &sys.n_tasks);

  // Rank 0 reads the input.
  if (sys.my_rank == 0) {
    timer_start("Total");
    timer_start("Input Read");
    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    input_buffer.natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    input_buffer.mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    input_buffer.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    input_buffer.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    input_buffer.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    input_buffer.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    input_buffer.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    input_buffer.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    input_buffer.nprint=atoi(line);
  }

  // Create MPI types for messaging inputs and broadcasts data.
  int array_count = 2;
  int array_blk[2]  = {3, 6};
  MPI_Aint array_disp[2] = {0, 3*sizeof(int)};
  MPI_Datatype array_typ[2] = {MPI_INT, MPI_DOUBLE};
  MPI_Datatype MD_INP;
  MPI_Type_struct(array_count, &array_blk[0], &array_disp[0], &array_typ[0],
                  &MD_INP);
  MPI_Type_commit(&MD_INP);
  MPI_Bcast((void*) &input_buffer, 1, MD_INP, 0, MPI_COMM_WORLD);

  // Copies variables from input buffer.
  sys.natoms  = input_buffer.natoms;
  sys.mass    = input_buffer.mass;
  sys.epsilon = input_buffer.epsilon;
  sys.sigma   = input_buffer.sigma;
  sys.rcut    = input_buffer.rcut;
  sys.box     = input_buffer.box;
  sys.nsteps  = input_buffer.nsteps;
  sys.nprint  = input_buffer.nprint;
  sys.dt      = input_buffer.dt;

  // Creates local atom list and arrays.
  create_atom_list(&sys);

  // Creates buffer for restart reading.
  sys.coordinates =
           (double *)malloc(sys.all_atoms*sizeof(double)*3);
  sys.vx = (double *)malloc(sys.all_atoms*sizeof(double));
  sys.vy = (double *)malloc(sys.all_atoms*sizeof(double));
  sys.vz = (double *)malloc(sys.all_atoms*sizeof(double));

  // Rank 0 reads restart and broadcasts velocities and coordinates to the rest.
  if (sys.my_rank == 0) {
    sys.rx = (double *)malloc(sys.all_atoms*sizeof(double));
    sys.ry = (double *)malloc(sys.all_atoms*sizeof(double));
    sys.rz = (double *)malloc(sys.all_atoms*sizeof(double));

    fp=fopen(restfile,"r");
    if(fp) {
      for (int i = 0; i < sys.all_atoms; ++i) {
        fscanf(fp, "%lf%lf%lf", sys.rx + i, sys.ry + i,
               sys.rz + i);
      }
      for (int i = 0; i < sys.all_atoms; ++i) {
        fscanf(fp, "%lf%lf%lf", sys.vx + i, sys.vy + i,
               sys.vz+i);
      }
      fclose(fp);
    } else {
      perror("cannot read restart file");
      return 3;
    }

    for (int i = 0; i < sys.all_atoms; i++) {
      sys.coordinates[i] = sys.rx[i];
    }
    for (int i = 0; i < sys.all_atoms; i++) {
      sys.coordinates[i + sys.all_atoms] = sys.ry[i];
    }
    for (int i = 0; i < sys.all_atoms; i++) {
      sys.coordinates[i + 2*(sys.all_atoms)] = sys.rz[i];
    }

    timer_stop("Input Read");
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
  }

  MPI_Datatype MD_CRDVEL;
  MPI_Type_contiguous(sys.all_atoms, MPI_DOUBLE, &MD_CRDVEL);
  MPI_Type_commit(&MD_CRDVEL);
  MPI_Bcast((void*) sys.coordinates, 3, MD_CRDVEL, 0, MPI_COMM_WORLD);
  MPI_Bcast((void*) sys.vx, 1, MD_CRDVEL, 0, MPI_COMM_WORLD);
  MPI_Bcast((void*) sys.vy, 1, MD_CRDVEL, 0, MPI_COMM_WORLD);
  MPI_Bcast((void*) sys.vz, 1, MD_CRDVEL, 0, MPI_COMM_WORLD);

  // Initializes other arrays
  sys.fx      = (double *)malloc(sys.my_atoms*sizeof(double));
  sys.fy      = (double *)malloc(sys.my_atoms*sizeof(double));
  sys.fz      = (double *)malloc(sys.my_atoms*sizeof(double));
  sys.n_ngb   =     (int*)malloc(sys.my_atoms*sizeof(int));
  sys.ngb_loc =     (int*)malloc(sys.my_atoms*sizeof(int));
  sys.ngb_list=     (int*)malloc(sys.my_atoms*NGB_MAX*sizeof(int));

  azzero(sys.fx, sys.my_atoms);
  azzero(sys.fy, sys.my_atoms);
  azzero(sys.fz, sys.my_atoms);
  azzero_i(sys.ngb_list, sys.my_atoms*NGB_MAX);

  // Creates datatypes for MPI_AllgatherV
  allgv_t allgv_data;
  ALLGV_setup(&sys, &allgv_data);

  // Initializes forces, energies and outputs.
  sys.nfi=0;
  make_ngb_list(&sys);
  force_ngb(&sys);
  ekin(&sys);

  erg=fopen(ergfile,"w");
  traj=fopen(trajfile,"w");
  write_initial_output(&sys, erg, traj);

  // MAIN LOOP
  for(sys.nfi=1; sys.nfi <= sys.nsteps; sys.nfi++) {
    if ((sys.nfi % sys.nprint) == 0) {
      if (sys.my_rank == 0) timer_start("Output Writing");
      write_output(&sys, erg, traj);
      if (sys.my_rank == 0) timer_pause("Output Writing");
    }
    if (sys.my_rank == 0) timer_start("Verlet Propagation");
    velverlet(&sys, &allgv_data);

    if (sys.my_rank == 0) {
     timer_pause("Verlet Propagation");
     timer_start("Energy Calculation");
    }
    ekin(&sys);
    if (sys.my_rank == 0)  timer_pause("Energy Calculation");
  }

  /* clean up: close files, free memory */
  if (sys.my_rank == 0) printf("Simulation Done.\n");

  free(sys.my_atom_list);
  free(sys.ghost_atom_list);
  free(sys.coordinates);
  free(sys.vx);
  free(sys.vy);
  free(sys.vz);
  free(sys.fx);
  free(sys.fy);
  free(sys.fz);
  free(sys.n_ngb);
  free(sys.ngb_loc);
  free(sys.ngb_list);
  free(allgv_data.recvbuf);
  free(allgv_data.sendbuf);
  free(allgv_data.recvcounts);
  free(allgv_data.displs);

  if (sys.my_rank == 0) {
    timer_stop("Total");
    print_timers();
    fclose(erg);
    //fclose(traj);
  }
  MPI_Finalize();
}
