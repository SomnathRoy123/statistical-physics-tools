#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>

#define MAX_WAVEVECTORS 40000
#define MAX_PARTICLES 100001
#define PI 3.14159265359
#define q_max 150
#define MAX_FILENAME 2000

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
    return -1;  // Return -1 if filename doesn't match expected format
}

int main() {
    VECTOR *particles = malloc(MAX_PARTICLES * sizeof(VECTOR));
    double *S_qc = calloc(MAX_WAVEVECTORS, sizeof(double));
    double *S_qs = calloc(MAX_WAVEVECTORS, sizeof(double));
    
    if (!particles || !S_qc || !S_qs) {
        fprintf(stderr, "Memory allocation failed\n");
        free(particles);
        free(S_qc);
        free(S_qs);
        return 1;
    }

    int num_particles = 0, snapshot_count = 0;
    double BoxSize_x =180.0, BoxSize_y = 180.0;
    double dq = 2 * PI / BoxSize_x;

    char directory_path[] = "/home/somnath2/Codes/Trial_115/TA400_R1.7/evo/";
    struct dirent *entry;
    DIR *dp = opendir(directory_path);

    if (!dp) {
        perror("Error opening directory_1");
        free(particles);
        free(S_qc);
        free(S_qs);
        return 1;
    }

    while ((entry = readdir(dp)) != NULL) {
        // Extract time number and check if it's in range 25-30
        int time_num = get_time_number(entry->d_name);
        if (time_num < 100|| time_num > 400) {
            continue;  // Skip files outside our time range
        }

        char file_path[MAX_FILENAME];
        if (snprintf(file_path, MAX_FILENAME, "%s%s", directory_path, entry->d_name) >= MAX_FILENAME) {
            fprintf(stderr, "File path too long\n");
            continue;
        }

        FILE *fp_input = fopen(file_path, "r");
        if (!fp_input) {
            perror("Error opening input file");
            continue;
        }

        num_particles = 0;
        double dummy;
        while (fscanf(fp_input, "%lf %lf %lf", &particles[num_particles].x, 
                     &particles[num_particles].y,&dummy) == 3) {
            num_particles++;
            if (num_particles >= MAX_PARTICLES) {
                fprintf(stderr, "Warning: Reached maximum particle limit (%d)\n", MAX_PARTICLES);
                break;
            }
        }
        fclose(fp_input);

        if (num_particles == 0) {
            fprintf(stderr, "Warning: No particles read from %s\n", file_path);
            continue;
        }

        // Compute Static Structure Factor S(k) for this snapshot
        for (int k = 0; k < MAX_WAVEVECTORS; k++) {
            int qx = k % q_max;
            int qy = k / q_max;
            
            double sqc = 0.0;
            double sqs = 0.0;

            for (int i = 0; i < num_particles; i++) {
                double dot_product = dq * ((qx + 1.0) * particles[i].x + 
                                        (qy + 1.0) * particles[i].y);
                sqc += cos(dot_product);
                sqs += sin(dot_product);
            }

            S_qc[k] += sqc * sqc;
            S_qs[k] += sqs * sqs;
        }

        snapshot_count++;
        printf("Processed snapshot %d: %s (%d particles)\n", 
               snapshot_count, entry->d_name, num_particles);
    }
    closedir(dp);

    if (snapshot_count == 0) {
        fprintf(stderr, "No snapshots in time range 25-30 processed. Exiting.\n");
        free(particles);
        free(S_qc);
        free(S_qs);
        return 1;
    }

    // Normalize structure factor by the number of snapshots
    for (int k = 0; k < MAX_WAVEVECTORS; k++) {
        S_qc[k] /= snapshot_count;
        S_qs[k] /= snapshot_count;
    }

    FILE *sk_file = fopen("/home/somnath2/Codes/Radius_st_fac_4/structure_factor_trial_115_400_1.7_t=100_400.dat", "w");
    if (!sk_file) {
        perror("Error opening output file");
        free(particles);
        free(S_qc);
        free(S_qs);
        return 1;
    }

    for (int k = 0; k < MAX_WAVEVECTORS; k++) {
        int qx = k % q_max;
        int qy = k / q_max;

        double k_magnitude = dq * sqrt((qx + 1.0) * (qx + 1.0) + 
                                     (qy + 1.0) * (qy + 1.0));
        fprintf(sk_file, "%lf %lf\n", k_magnitude, 
                (S_qc[k] + S_qs[k]) / ((double)num_particles));
    }
    fclose(sk_file);

    printf("Structure factor averaged over %d snapshots (time 25-30) written to file.\n", 
           snapshot_count);

    free(particles);
    free(S_qc);
    free(S_qs);
    return 0;
}