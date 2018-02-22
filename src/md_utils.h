#ifndef __MDUH__
#define __MDUH__
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
  int ngb_freq = 20;   // Frequency for neighbour list update.      

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

struct _ALLGV_struct {
  double *sendbuf;
  double *recvbuf;
  int    *recvcounts;
  int    *displs;
};
typedef struct _ALLGV_struct allgv_t;

extern void azzero(double *d, const int n);
extern void azzero_i(int *d, const int n);
extern double pbc(double x, const double boxby2);
extern void ALLGV_setup(mdsys_t *sys, allgv_t *allgv_data);
#endif
