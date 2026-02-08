/*
This is the code ( among all codes in this file Code_rdf ) to calcualte the Radial Distribution function)
Things to modify and check before running the code 

1. Box_size depends on the file you are running 
2. time num from what file num to what you want to run your code 
3. Directory path 
4. Output path
 */ 



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>

#define MAX_PARTICLES 100001
#define MAX_BINS 1000
#define PI 3.14159265359
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
    if (!particles) {
        fprintf(stderr, "Memory allocation failed for particles\n");
        return 1;
    }

    int num_particles = 0, snapshot_count = 0;
    double BoxSize_x = 180.0, BoxSize_y = 180.0;
    
    // RDF parameters
    double r_max = 20;  // Maximum radius (half of box diagonal)   
    int num_bins = 400;  // Number of bins for the histogram
    double dr = r_max / num_bins;  // Bin width
    
    // Allocate memory for RDF histogram
    double *g_r = calloc(num_bins, sizeof(double));
    if (!g_r) {
        fprintf(stderr, "Memory allocation failed for g_r\n");
        free(particles);
        return 1;
    }

    char directory_path[] = "/home/somnath2/Codes/Trial_88/TA160_R3.5/evo/";
    struct dirent *entry;
    DIR *dp = opendir(directory_path);

    if (!dp) {
        perror("Error opening directory");
        free(particles);
        free(g_r);
        return 1;
    }

    // Process each snapshot in the directory
    while ((entry = readdir(dp)) != NULL) {
        // Extract time number and check if it's in range 120-170
        int time_num = get_time_number(entry->d_name);
        if (time_num < 1000|| time_num > 2000) {
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
                     &particles[num_particles].y, &dummy) == 3) {

        // while (fscanf(fp_input, "%lf %lf", &particles[num_particles].x, 
        //     &particles[num_particles].y) == 2) {
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

        // Calculate RD/F for this snapshot
        for (int i = 0; i < num_particles - 1; i++) {
            for (int j = i + 1; j < num_particles; j++) {
                // Calculate distance between particles with periodic boundary conditions
                double dx = particles[i].x - particles[j].x;
                dx = dx - BoxSize_x * round(dx / BoxSize_x);  // Apply periodic boundary
                
                double dy = particles[i].y - particles[j].y;
                dy = dy - BoxSize_y * round(dy / BoxSize_y);  // Apply periodic boundary
                
                double r = sqrt(dx*dx + dy*dy);
                
                // Add to histogram if within range
                if (r < r_max) {
                    int bin = (int)(r / dr);
                    g_r[bin] += 2.0;  // Count each pair twice (i->j and j->i)
                }
            }
        }

        snapshot_count++;
        printf("Processed snapshot %d: %s (%d particles)\n", 
               snapshot_count, entry->d_name, num_particles);
    }
    closedir(dp);

    if (snapshot_count == 0) {
        fprintf(stderr, "No snapshots in time range 120-170 processed. Exiting.\n");
        free(particles);
        free(g_r);
        return 1;
    }

    // Calculate the 2D density (particles per unit area)
    double area = BoxSize_x * BoxSize_y;
    double number_density = num_particles / area;
    
    // Normalize RDF
    FILE *rdf_file = fopen("/home/somnath2/Codes/rdf_data_5/rdf_trial_88_160_3.5_t1000_2000.dat", "w");
    if (!rdf_file) {
        perror("Error opening output file");
        free(particles);
        free(g_r);
        return 1;
    }

    // Normalize and write the RDF
    for (int i = 0; i < num_bins; i++) {
        double r = (i + 0.5) * dr;  // Radius at the bin center
        
        // Normalize by:
        // 1. Number of snapshots
        // 2. Number of particles
        // 3. Expected number of particles in an ideal gas (2πr dr × number_density)
        double normalization = snapshot_count * num_particles * (2.0 * PI * r * dr * number_density);
        
        if (normalization > 0) {
            double g_r_normalized = g_r[i] / normalization;
            fprintf(rdf_file, "%lf %lf\n", r, g_r_normalized);
        }
    }
    
    fclose(rdf_file);

    printf("Radial distribution function averaged over %d snapshots (time 120-170) written to file.\n", 
           snapshot_count);

    free(particles);
    free(g_r);
    return 0;
}