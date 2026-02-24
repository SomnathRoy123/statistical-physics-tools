#include "graccum.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

typedef struct {
    double r_center;
    double shell_area_sum;   /* sum over frames of annulus area for this bin */
    double ideal_pairs_sum;  /* sum over frames of ideal unordered pair count in shell */
    long   pair_count;       /* observed unordered pairs */
} GrBin;

struct GrAccum {
    GrBin *bins;
    int    nbins;
    double dr;
    long   nframes;
};

static double annulus_area(double r_in, double r_out){
    return M_PI * (r_out * r_out - r_in * r_in);
}

GrAccum *graccum_create(double dr){
    if(dr <= 0.0){
        fprintf(stderr, "graccum_create: dr must be > 0\n");
        return NULL;
    }
    GrAccum *A = (GrAccum*)malloc(sizeof(GrAccum));
    if(!A){
        fprintf(stderr, "graccum_create: OOM\n");
        return NULL;
    }
    A->bins = NULL;
    A->nbins = 0;
    A->dr = dr;
    A->nframes = 0;
    return A;
}

void graccum_free(GrAccum *A){
    if(!A) return;
    free(A->bins);
    A->bins = NULL;
    A->nbins = 0;
    free(A);
}

static void graccum_ensure_bins(GrAccum *A, int bmax){
    if(!A || bmax < 0) return;
    if(bmax < A->nbins) return;

    int new_n = bmax + 1;
    GrBin *nb = (GrBin*)realloc(A->bins, (size_t)new_n * sizeof(GrBin));
    if(!nb){
        fprintf(stderr, "graccum_ensure_bins: OOM\n");
        exit(1);
    }

    for(int b = A->nbins; b < new_n; ++b){
        nb[b].r_center = (b + 0.5) * A->dr;
        nb[b].shell_area_sum = 0.0;
        nb[b].ideal_pairs_sum = 0.0;
        nb[b].pair_count = 0;
    }

    A->bins = nb;
    A->nbins = new_n;
}

void graccum_accumulate(GrAccum *A,
                        const Vec2Array *coms,
                        bool use_pbc,
                        double box_x,
                        double box_y)
{
    if(!A || !coms) return;
    int M = (int)coms->n;
    if(M < 2) return;

    if(use_pbc && (box_x <= 0.0 || box_y <= 0.0)){
        fprintf(stderr, "graccum_accumulate: invalid box dimensions for PBC\n");
        return;
    }

    /* First pass: frame rmax for dynamic bin growth */
    double rmax2 = 0.0;
    for(int i = 0; i < M - 1; ++i){
        for(int j = i + 1; j < M; ++j){
            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double r2 = dx*dx + dy*dy;
            if(r2 > rmax2) rmax2 = r2;
        }
    }

    double rmax = sqrt(rmax2);
    int bmax = (int)floor(rmax / A->dr);
    graccum_ensure_bins(A, bmax);

    /* Pair counting */
    for(int i = 0; i < M - 1; ++i){
        for(int j = i + 1; j < M; ++j){
            double dx = coms->data[j].x - coms->data[i].x;
            double dy = coms->data[j].y - coms->data[i].y;
            if(use_pbc){
                dx = mic_delta(dx, box_x);
                dy = mic_delta(dy, box_y);
            }
            double r = sqrt(dx*dx + dy*dy);
            int b = (int)floor(r / A->dr);
            if(b < 0 || b >= A->nbins) continue;
            A->bins[b].pair_count += 1;
        }
    }

    /* Add shell geometry + ideal-gas expectation for this frame */
    double area_box = use_pbc ? (box_x * box_y) : 0.0;
    double rho = (area_box > 0.0) ? ((double)M / area_box) : 0.0;

    for(int b = 0; b < A->nbins; ++b){
        double r_in = b * A->dr;
        double r_out = (b + 1) * A->dr;
        double shell_area = annulus_area(r_in, r_out);

        A->bins[b].shell_area_sum += shell_area;

        if(rho > 0.0){
            /* unordered ideal pairs in shell: 0.5 * N * rho * shell_area */
            double ideal_pairs = 0.5 * (double)M * rho * shell_area;
            A->bins[b].ideal_pairs_sum += ideal_pairs;
        }
    }

    A->nframes += 1;
}

int graccum_write(GrAccum *A,
                  const char *outpath,
                  int t0, int t1,
                  bool use_pbc,
                  double box_x,
                  double box_y)
{
    if(!A || !outpath){
        fprintf(stderr, "graccum_write: invalid args\n");
        return 1;
    }

    FILE *f = fopen(outpath, "w");
    if(!f){
        fprintf(stderr, "graccum_write: cannot open %s: %s\n", outpath, strerror(errno));
        return 2;
    }

    fprintf(f, "# Averaged g(r) over snapshots time_%d .. time_%d\n", t0, t1);
    fprintf(f, "# Columns: r_center pair_count shell_area pair_density ideal_pairs g_r\n");
    fprintf(f, "# Params: dr=%.8g USE_PBC=%s frames=%ld\n",
            A->dr, use_pbc ? "true" : "false", A->nframes);
    if(use_pbc){
        fprintf(f, "# Box dims: %.8g %.8g\n", box_x, box_y);
    }

    for(int b = 0; b < A->nbins; ++b){
        double shell_area = A->bins[b].shell_area_sum;
        double pairs = (double)A->bins[b].pair_count;
        double pair_density = (shell_area > 0.0) ? (pairs / shell_area) : 0.0;
        double ideal_pairs = A->bins[b].ideal_pairs_sum;
        double g_r = (ideal_pairs > 0.0) ? (pairs / ideal_pairs) : 0.0;

        fprintf(f, "%.10g %.10g %.10g %.10g %.10g %.10g\n",
                A->bins[b].r_center,
                pairs,
                shell_area,
                pair_density,
                ideal_pairs,
                g_r);
    }

    fclose(f);
    return 0;
}
