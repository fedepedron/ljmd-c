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
#include "timers.h"
#include "ncdf.h"
#include <mpi.h>

/* generic file- or pathname buffer length */
#define BLEN 200
#define NGB_MAX 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
  int natoms, nfi, nsteps;
  int nprint;
  double dt, mass, epsilon, sigma, box, rcut;
  double ekin, epot, temp;
  double *rx, *ry, *rz; // These are only used in input read.
  double *vx, *vy, *vz;
  double *fx, *fy, *fz;
  double *coordinates;  // XYZ coordinates for each atom.

  int *n_ngb;           // Number of total neighbours
  int *ngb_loc;         // Number of neighbours that are not ghosts.
  int *ngb_list;        // List of original IDs of neighbours.

  int my_atoms;         // N° of atoms in current task.
  int ghost_atoms;      // N° of ghost atoms in current task.
  int all_atoms;        // my_atoms + ghost_atoms;
  int *my_atom_list;    // List of original IDs of atoms in current task.
  int *ghost_atom_list; // List of original IDs of ghost atoms in current task.
  int atom_count_L;     // N° of atoms in tasks with id > natoms%ntasks.
  int atom_count_S;     // N° of atoms in tasks with id < natoms%ntasks.

  int my_rank;
  int n_tasks;
  int task_rest;        // Equals to natoms%ntasks.
};
typedef struct _mdsys mdsys_t;

struct __attribute__((packed)) _mdsysb {
  int    natoms, nsteps, nprint;
  double dt, mass, epsilon, sigma, box, rcut;
};
typedef struct _mdsysb mdsys_inp_t;

//##############################################################################
// HELPER FUNCTIONS
// Read a line and then return the first string with whitespace stripped off.
static int get_a_line(FILE *fp, char *buf)
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

// These zero arrays of int and double.
static void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

static void azzero_i(int *a, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        a[i]=0.0;
    }
}

// Apply minimum image convention
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

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

      if (iloc < jloc) continue;
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
      sys->fx[jloc] -= sys->fx[iloc];
      sys->fy[jloc] -= sys->fy[iloc];
      sys->fz[jloc] -= sys->fz[iloc];
    }
    if ((sys->my_rank == 0) && iloc < 1) printf("Energy %10.6f\n", epot);
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
    if ((sys->my_rank == 0) && iloc < 1) printf("Energy %10.6f\n", epot);
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
      double r     = sqrt(rx*rx + ry*ry + rz*rz);
      if (r < sys->rcut) {
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
      double r     = sqrt(rx*rx + ry*ry + rz*rz);

      if (r < sys->rcut) {
        sys->ngb_list[i*NGB_MAX + sys->n_ngb[i]] = j;
        sys->n_ngb[i]++;
      }
    }
  }
  return;
}

static void task_update_coords(mdsys_t *sys)
{

}


/* velocity verlet */
static void velverlet(mdsys_t *sys)
{
  /* first part: propagate velocities by half and positions by full step */
  for (int id=0; id<sys->my_atoms; ++id) {
    int i = sys->my_atom_list[id];
    sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
    sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
    sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    sys->coordinates[i] += sys->dt*sys->vx[i];
    sys->coordinates[i + sys->all_atoms] += sys->dt*sys->vy[i];
    sys->coordinates[i + 2*(sys->all_atoms)] += sys->dt*sys->vz[i];
  }

  /* compute forces and potential energy */
  if (sys->my_rank == 0) timer_start("Neighbour List");
  make_ngb_list(sys);

  if (sys->my_rank == 0) {
    timer_pause("Neighbour List");
    timer_start("Force Calculation");
  }

  force_ngb(sys);
  if (sys->my_rank == 0) timer_pause("Force Calculation");

  /* second part: propagate velocities by another half step */
  for (int id=0; id<sys->my_atoms; ++id) {
    int i = sys->my_atom_list[id];
    sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
    sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
    sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
  }
}

/* append data to output. */
static void output(mdsys_t *sys, FILE *erg, FILE *traj, char* trajfile)
{
  printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
         sys->ekin, sys->epot, sys->ekin+sys->epot);
  fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp,
          sys->ekin, sys->epot, sys->ekin+sys->epot);
  // Print trajectory using NetCDF format.
  write_to_netcdf(sys->rx, sys->ry, sys->rz, sys->natoms, trajfile);
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

  if (sys.my_rank == 1) {
    printf("HOLI %d %d %d  %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f ",
      sys.natoms, sys.nprint, sys.nsteps, sys.mass, sys.epsilon, sys.sigma,
     sys.rcut, sys.box, sys.dt);
  }
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
  azzero_i(sys.ngb_list, sys.my_atoms*sys.all_atoms);

  // Initializes forces, energies and outputs.
  sys.nfi=0;
  make_ngb_list(&sys);
  force_ngb(&sys);
  ekin(&sys);
  printf("My energy %10.6f \n",sys.epot );
  if (sys.my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, (void*) &sys.epot, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, (void*) &sys.ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj, trajfile);
  } else {
    MPI_Reduce((void*) &sys.epot, (void*) &sys.epot, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce((void*) &sys.ekin, (void*) &sys.ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }
  return 0;
  // MAIN LOOP
  //timer_start("Main Loop");
  for(sys.nfi=1; sys.nfi <= sys.nsteps; sys.nfi++) {
    if ((sys.nfi % sys.nprint) == 0) {
      if (sys.my_rank == 0) {
        timer_start("Output Writing");
        MPI_Reduce(MPI_IN_PLACE, (void*) &sys.epot, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, (void*) &sys.ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        output(&sys, erg, traj, trajfile);
        timer_pause("Output Writing");
      } else {
        MPI_Reduce((void*) &sys.epot, (void*) &sys.epot, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce((void*) &sys.ekin, (void*) &sys.ekin, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
      }
    }
    if (sys.my_rank == 0)  timer_start("Coordinate Messaging");
    task_update_coords(&sys);

    if (sys.my_rank == 0)  {
      timer_pause("Coordinate Messaging");
      timer_start("Verlet Propagation");
    }

    velverlet(&sys);

    if (sys.my_rank == 0) {
     timer_pause("Verlet Propagation");
     timer_start("Energy Calculation");
    }

    ekin(&sys);

    if (sys.my_rank == 0)  timer_pause("Energy Calculation");

  }
  //timer_stop("Main Loop");

  /* clean up: close files, free memory */
  printf("Simulation Done. Rank: %d.\n", sys.my_rank);

  free(sys.my_atom_list);
  free(sys.ghost_atom_list);
  free(sys.coordinates);
  free(sys.rx);
  free(sys.ry);
  free(sys.rz);
  free(sys.vx);
  free(sys.vy);
  free(sys.vz);
  free(sys.fx);
  free(sys.fy);
  free(sys.fz);
  free(sys.n_ngb);
  free(sys.ngb_loc);
  free(sys.ngb_list);

  if (sys.my_rank == 0) {
    timer_stop("Total");
    print_timers();
    fclose(erg);
    fclose(traj);
  }
  MPI_Finalize();
}
