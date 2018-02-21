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
  int natoms,nfi,nsteps;
  double dt, mass, epsilon, sigma, box, rcut;
  double ekin, epot, temp;
  double *rx, *ry, *rz;
  double *vx, *vy, *vz;
  double *fx, *fy, *fz;
  int *n_ngb;
  int *ngb_list;
  double *d_ngb;
};
typedef struct _mdsys mdsys_t;

struct _mdsysb {
  int natoms, nsteps;
  double dt, mass, epsilon, sigma, box, rcut;
};
typedef struct _mdsysb mdsys_inp_t;

struct _crdvel_buf {
  double *rx, *ry, *rz;
  double *vx, *vy, *vz;
};
typedef struct _crdvel_buf crdvel_buf_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
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

/* helper function: zero out an array */
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

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute kinetic energy */
static void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i]
                     + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

/* compute forces */
static void force_ngb(mdsys_t *sys)
{
    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    double epot = 0.0;

    #pragma omp parallel for reduction(+:epot)
    for(int i=0; i < (sys->natoms); i++) {
        double fx   = sys->fx[i];
        double fy   = sys->fy[i];
        double fz   = sys->fz[i];
        double loc_epot = 0.0;

        #pragma omp parallel for reduction(+:fx, fy, fz, loc_epot)
        for(int ingb=0; ingb < (sys->n_ngb[i]); ingb++) {
            int j = sys->ngb_list[ingb + i*NGB_MAX];
            double rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            double ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            double rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            double r = sqrt(rx*rx + ry*ry + rz*rz);

            double ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);
            loc_epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                          -pow(sys->sigma/r,6.0));
            fx += rx/r*ffac;
            fy += ry/r*ffac;
            fz += rz/r*ffac;
        }

        sys->fx[i] = fx;
        sys->fy[i] = fy;
        sys->fz[i] = fz;
        epot += loc_epot;
    }

    sys->epot = epot;
    return;
}

/* Create Neighbour list using cutoff */
void make_ngb_list(mdsys_t *sys)
{
  if (sys->natoms < 15000){
    for(int i=0; i < (sys->natoms); i++) {
      int tmp_ngb = 0;

      #pragma omp parallel for
      for(int j=0; j < (sys->natoms); j++) {
        if (i==j) continue;
        double rx  = pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
        double ry  = pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
        double rz  = pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
        sys->d_ngb[i*(sys->natoms) + j] = sqrt(rx*rx + ry*ry + rz*rz);
      }

      for(int j=0; j < (sys->natoms); j++) {
        if (i==j) continue;
        if (sys->d_ngb[i*(sys->natoms) + j] < sys->rcut) {
          sys->ngb_list[i*NGB_MAX + tmp_ngb] = j;
          tmp_ngb++;
        }
      }
      sys->n_ngb[i] = tmp_ngb;
    }
  } else {
    azzero_i(sys->n_ngb, sys->natoms);
    for(int i=0; i < (sys->natoms); i++) {
      for(int j=0; j < i; j++) {
        double rx  = pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
        double ry  = pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
        double rz  = pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
        double r   = sqrt(rx*rx + ry*ry + rz*rz);

        if (r < sys->rcut) {
          sys->ngb_list[i*NGB_MAX + sys->n_ngb[i]] = j;
          sys->ngb_list[j*NGB_MAX + sys->n_ngb[j]] = j;
          sys->n_ngb[i]++;
          sys->n_ngb[j]++;
        }
      }
    }
  }
  return;
}

/* velocity verlet */
static void velverlet(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }

    /* compute forces and potential energy */
    make_ngb_list(sys);
    force_ngb(sys);

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
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
    //fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi,
    // sys->ekin+sys->epot);
    //for (int i=0; i<sys->natoms; ++i) {
    //    fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i],
    // sys->rz[i]);
    //}
}


/* main */
int main(int argc, char **argv)
{
  MPI_Init(NULL,NULL);

  char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
  int nprint;
  FILE *fp,*traj,*erg;
  int my_rank=1, n_tasks;
  mdsys_t sys;
  mdsys_inp_t input_buffer;
  crdvel_buf_t crdvel_buffer;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);

  // Rank 0 reads the input.
  if (my_rank == 0) {
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
    nprint=atoi(line);
  }

  // Create MPI types for messaging inputs and broadcasts data.
  int array_count = 2;
  int array_blk[2]  = {2, 6};
  MPI_Aint array_disp[2] = {0, 0};
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
  sys.dt      = input_buffer.dt;

  // Creates buffer for restart reading.
  crdvel_buffer.rx = (double *)malloc(sys.natoms*sizeof(double));
  crdvel_buffer.ry = (double *)malloc(sys.natoms*sizeof(double));
  crdvel_buffer.rz = (double *)malloc(sys.natoms*sizeof(double));
  crdvel_buffer.vx = (double *)malloc(sys.natoms*sizeof(double));
  crdvel_buffer.vy = (double *)malloc(sys.natoms*sizeof(double));
  crdvel_buffer.vz = (double *)malloc(sys.natoms*sizeof(double));

  // Rank 0 reads restart.
  if (my_rank == 0) {
    fp=fopen(restfile,"r");
    if(fp) {
      for (int i=0; i<sys.natoms; ++i) {
        fscanf(fp, "%lf%lf%lf", crdvel_buffer.rx + i, crdvel_buffer.ry + i,
               crdvel_buffer.rz + i);
      }
      for (int i=0; i<sys.natoms; ++i) {
        fscanf(fp, "%lf%lf%lf", crdvel_buffer.vx + i, crdvel_buffer.vy + i,
               crdvel_buffer.vz+i);
      }
      fclose(fp);
    } else {
      perror("cannot read restart file");
      return 3;
    }
    timer_stop("Input Read");
  }

  // Rank 0 broadcasts velocities and coordinates to the rest.
  MPI_Datatype MD_CRDVEL;
  //
  int array_countb = 1;
  int array_blkb[2]  = {6*sys.natoms};
  MPI_Aint array_dispb[1] = {0};
  MPI_Datatype array_typb[1] = {MPI_DOUBLE};
  MPI_Type_struct(array_countb, &array_blkb[0], &array_dispb[0], &array_typb[0],
                  &MD_CRDVEL);
  //MPI_Type_contiguous(sys.natoms*6, MPI_DOUBLE, &MD_CRDVEL);
  MPI_Type_commit(&MD_CRDVEL);
  MPI_Bcast((void*) &crdvel_buffer, 1, MD_CRDVEL, 0, MPI_COMM_WORLD);

  /* allocate memory */
  sys.rx      = (double *)malloc(sys.natoms*sizeof(double));
  sys.ry      = (double *)malloc(sys.natoms*sizeof(double));
  sys.rz      = (double *)malloc(sys.natoms*sizeof(double));
  sys.vx      = (double *)malloc(sys.natoms*sizeof(double));
  sys.vy      = (double *)malloc(sys.natoms*sizeof(double));
  sys.vz      = (double *)malloc(sys.natoms*sizeof(double));
  sys.fx      = (double *)malloc(sys.natoms*sizeof(double));
  sys.fy      = (double *)malloc(sys.natoms*sizeof(double));
  sys.fz      = (double *)malloc(sys.natoms*sizeof(double));
  sys.n_ngb   =     (int*)malloc(sys.natoms*sizeof(int));
  sys.ngb_list=     (int*)malloc(sys.natoms*NGB_MAX*sizeof(int));

  if (sys.natoms < 15000)
    sys.d_ngb = (double *)malloc(sys.natoms*sys.natoms*sizeof(double));

  for (int i=0; i < sys.natoms; i++) {
    if (my_rank == 1) printf("I am at %d \n", i);
    sys.rx[i] = crdvel_buffer.rx[i];
    sys.ry[i] = crdvel_buffer.ry[i];
    sys.rz[i] = crdvel_buffer.rz[i];
    sys.vx[i] = crdvel_buffer.vx[i];
    sys.vy[i] = crdvel_buffer.vy[i];
    sys.vz[i] = crdvel_buffer.vz[i];
  }
  printf("equal %d \n", my_rank);

  MPI_Barrier(MPI_COMM_WORLD);

  if(my_rank == 0) {
    free(crdvel_buffer.rx);
    free(crdvel_buffer.ry);
    free(crdvel_buffer.rz);
    free(crdvel_buffer.vx);
    free(crdvel_buffer.vy);
    free(crdvel_buffer.vz);
  }

  printf("Freed %d \n", my_rank);

  azzero(sys.fx, sys.natoms);
  azzero(sys.fy, sys.natoms);
  azzero(sys.fz, sys.natoms);
  azzero_i(sys.n_ngb, sys.natoms);
  azzero_i(sys.ngb_list, sys.natoms*sys.natoms);

  /* initialize forces and energies.*/
  sys.nfi=0;
  make_ngb_list(&sys);
  force_ngb(&sys);
  ekin(&sys);

  if (my_rank == 0) {
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj, trajfile);
  }

  printf("Main loop %d \n", my_rank);

  // MAIN LOOP
  //timer_start("Main Loop");
  for(sys.nfi=1; sys.nfi <= sys.nsteps; sys.nfi++) {
    if (my_rank == 0) {
      timer_start("Output Writing");
      if ((sys.nfi % nprint) == 0) output(&sys, erg, traj, trajfile);
      timer_pause("Output Writing");
      timer_start("Verlet Propagation");
    }

    velverlet(&sys);

    if (my_rank == 0) {
     timer_pause("Verlet Propagation");
     timer_start("Energy Calculation");
    }

    ekin(&sys);

    if (my_rank == 0)  timer_pause("Energy Calculation");

    }
    //timer_stop("Main Loop");
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);

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
    free(sys.ngb_list);
    if (sys.natoms < 15000) free(sys.d_ngb);

    timer_stop("Total");
    print_timers();
  MPI_Finalize();
}
