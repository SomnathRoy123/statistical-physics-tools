#ifndef GTACCUM_H
#define GTACCUM_H

#include "utils.h"
#include <stdbool.h>

typedef struct GtAccum GtAccum;

/* Create/destroy */
GtAccum *gtaccum_create(double dr);
void gtaccum_free(GtAccum *A);

/* Accumulate translational correlation for one frame. */
void gtaccum_accumulate(GtAccum *A,
                        const Vec2Array *coms,
                        double theta_G,
                        double a_lattice,
                        bool use_pbc,
                        double box_x,
                        double box_y);

/* Write averaged g_T(r): r_center  gT_avg  pair_count */
int gtaccum_write(GtAccum *A,
                  const char *outpath,
                  int t0, int t1,
                  double a_lattice,
                  bool use_pbc,
                  double box_x,
                  double box_y);

#endif /* GTACCUM_H */
