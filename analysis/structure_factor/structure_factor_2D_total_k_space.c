#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>


#define MAX_PARTICLES 100001
#define PI 3.14159265359
#define q_max 600 // Number of k-points in each dimension (q_max x q_max grid)
#define MAX_FILENAME 2500
// The required size is q_max * q_max for the 2D k-space grid.
#define MAX_WAVEVECTORS (q_max * q_max) 

typedef struct {
    double x;
    double y;
} VECTOR;

// Function to extract time number from filename
int get_time_number(const char* filename) {
    int time_num;
    if (sscanf(filename, "time_%d.dat", &time_num) == 1) {
        return time_num;
    }
    return -1;
}

int main() {
    // Allocate memory
    VECTOR *particles = malloc(MAX_PARTICLES * sizeof(VECTOR));
    double *S_qc = calloc(MAX_WAVEVECTORS, sizeof(double));
    double *S_qs = calloc(MAX_WAVEVECTORS, sizeof(double));
    
    if (!particles || !S_qc || !S_qs) {
        fprintf(stderr, "Memory allocation failed\n");
        free(particles); free(S_qc); free(S_qs);
        return 1;
    }

    int num_particles = 0, snapshot_count = 0;
    double BoxSize_x = 180.0, BoxSize_y = 180.0;
    double dq = 2 * PI / BoxSize_x;

    // --- Directory and File Processing ---
    char directory_path[] = "/home/somnath2/Codes/Trial_102/TA160_R3.3/evo/";
    struct dirent *entry;
    DIR *dp = opendir(directory_path);

    if (!dp) {
        perror("Error opening directory");
        free(particles); free(S_qc); free(S_qs);
        return 1;
    }

    while ((entry = readdir(dp)) != NULL) {
        int time_num = get_time_number(entry->d_name);
        if (time_num < 1500 || time_num > 1700) {
            continue;
        }

        char file_path[MAX_FILENAME];
        snprintf(file_path, MAX_FILENAME, "%s%s", directory_path, entry->d_name);

        FILE *fp_input = fopen(file_path, "r");
        if (!fp_input) continue;

        num_particles = 0;
        double dummy;
        while (fscanf(fp_input, "%lf %lf %lf", &particles[num_particles].x, 
                     &particles[num_particles].y, &dummy) == 3) {
            num_particles++;
            if (num_particles >= MAX_PARTICLES) break;
        }
        fclose(fp_input);

        if (num_particles == 0) continue;

        // --- Structure Factor Calculation for the current snapshot ---
        // Loop over a 2D grid of wavevectors, including negative values
        for (int qy_idx = 0; qy_idx < q_max; qy_idx++) {
            for (int qx_idx = 0; qx_idx < q_max; qx_idx++) {
                
                // Map indices to a symmetric range, e.g., -100 to 99 for q_max=200
                int qx = qx_idx - q_max / 2;
                int qy = qy_idx - q_max / 2;

                if (qx == 0 && qy == 0) {
                    continue;
                }
                
                // Get the 1D index for storage
                int k = qy_idx * q_max + qx_idx;

                double sqc_temp = 0.0;
                double sqs_temp = 0.0;

                for (int i = 0; i < num_particles; i++) {
                    double dot_product = dq * (qx * particles[i].x + qy * particles[i].y);
                    sqc_temp += cos(dot_product);
                    sqs_temp += sin(dot_product);
                }

                S_qc[k] += sqc_temp * sqc_temp;
                S_qs[k] += sqs_temp * sqs_temp;
            }
        }

        snapshot_count++;
        printf("Processed snapshot %d: %s (%d particles)\n", 
               snapshot_count, entry->d_name, num_particles);
    }
    closedir(dp);

    if (snapshot_count == 0) {
        fprintf(stderr, "No snapshots processed. Exiting.\n");
        free(particles); free(S_qc); free(S_qs);
        return 1;
    }

    // --- Averaging and Writing Output ---
    for (int k = 0; k < MAX_WAVEVECTORS; k++) {
        S_qc[k] /= snapshot_count;
        S_qs[k] /= snapshot_count;
    }

    FILE *sk_file = fopen("/home/somnath2/Codes/Radius_st_fac_4/structure_factor_2D_full_kspace_trail65_160_3.3_t_1500_1700.dat", "w");
    if (!sk_file) {
        perror("Error opening output file");
        free(particles); free(S_qc); free(S_qs);
        return 1;
    }

    fprintf(sk_file, "# kx\tky\tS(k)\n");

    // Loop through wavevectors again to write the final data
    for (int qy_idx = 0; qy_idx < q_max; qy_idx++) {
        for (int qx_idx = 0; qx_idx < q_max; qx_idx++) {
            
            // Map indices to the same symmetric range
            int qx = qx_idx - q_max / 2;
            int qy = qy_idx - q_max / 2;

            int k = qy_idx * q_max + qx_idx;
            // Skip writing the k=0 (origin) point to the file
            if (qx == 0 && qy == 0) {
                continue;
            }


            double kx = dq * qx;
            double ky = dq * qy;
            
            double Sk = (S_qc[k] + S_qs[k]) / ((double)num_particles);

            fprintf(sk_file, "%lf %lf %lf\n", kx, ky, Sk);
        }
    }
    fclose(sk_file);

    printf("\nFull 2D structure factor data averaged over %d snapshots written.\n", snapshot_count);
    printf("Output now includes negative kx and ky values.\n");

    free(particles);
    free(S_qc);
    free(S_qs);
    return 0;
}