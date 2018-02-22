#ifndef __IOHEAD__
#define __IOHEAD__
#include <stdio.h>
#include "md_utils.h"

void write_output(mdsys_t *sys, FILE *erg, FILE *traj);
void write_initial_output(mdsys_t *sys, FILE *erg, FILE *traj);
int get_a_line(FILE *fp, char *buf);
#endif
