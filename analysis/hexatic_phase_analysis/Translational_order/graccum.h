#ifndef GRACCUM_H
#define GRACCUM_H

#include "utils.h"   /* Vec2Array, mic_delta */
#include <stdbool.h>

/* Opaque accumulator */
typedef struct GrAccum GrAccum;

/* Create/destroy */
GrAccum *graccum_create(double dr);
void graccum_free(GrAccum *A);

/* Accumulate one snapshot's pair-distance statistics for g(r). */
void graccum_accumulate(GrAccum *A,
                        const Vec2Array *coms,
                        bool use_pbc,
                        double box_x,
                        double box_y);

/* Write averaged g(r):
 * columns:
 *   r_center  pair_count  shell_area  pair_density  ideal_pairs  g_r
 *
 * pair_density = pair_count / shell_area
 * ideal_pairs  = accumulated ideal-gas pair expectation in that shell
 * g_r          = pair_count / ideal_pairs (if ideal_pairs > 0)
 */
int graccum_write(GrAccum *A,
                  const char *outpath,
                  int t0, int t1,
                  bool use_pbc,
                  double box_x,
                  double box_y);

#endif /* GRACCUM_H */
