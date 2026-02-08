#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <string.h>

#define MAX_PARTICLES 100001
#define PI 3.14159265359
#define MAX_FILENAME 2000

#define IDX(r, p, t, p_bins, t_bins) ((r) * (p_bins) * (t_bins) + (p) * (t_bins) + (t))

typedef struct {
    double x;
    double y;
    double orientation;
} PARTICLE;

// safe_acos is no longer needed for theta, but still good practice to keep
double safe_acos(double value) {
    if (value <= -1.0) return PI;
    if (value >= 1.0) return 0.0;
    return acos(value);
}

int get_time_number(const char* filename) {
    int time_num;
    if (sscanf(filename, "time_%d.dat", &time_num) == 1) {
        return time_num;
    }
    return -1;
}

int main() {
    // ========== USER-DEFINED PARAMETERS ==========
    double BoxSize_x = 182.0;
    double BoxSize_y = 182.0;
    int start_time = 1500;
    int end_time = 2000;
    char directory_path[] = "/home/somnath2/Codes/Trial_111/TA160_R3.5/evo/";
    char output_path[] = "/home/somnath2/Codes/rdf_data_3D_1/RDF_data_3D_Trial_111_160_3.5_t2.dat";

    // Binning parameters
    double r_max = 5.0;
    int num_r_bins = 50;

    double phi_max = 2.0 * PI; // Full circle for phi
    int num_phi_bins = 72;

    // MODIFIED PARAMETERS FOR THETA
    double theta_max = 2.0 * PI; // Full circle for theta
    int num_theta_bins = 4;
    // =============================================

    // ... (Memory allocation and file reading setup remains the same) ...
    PARTICLE *particles = malloc(MAX_PARTICLES * sizeof(PARTICLE));
    double *Histogram = calloc((long)num_r_bins * num_phi_bins * num_theta_bins, sizeof(double));
    if (!particles || !Histogram) { /* handle error */ return 1; }
    double dr = r_max / num_r_bins;
    double dphi = phi_max / num_phi_bins;
    double dtheta = theta_max / num_theta_bins;
    int snapshot_count = 0;
    int num_particles = 0;
    DIR *dp = opendir(directory_path);
    if (!dp) { /* handle error */ return 1; }
    struct dirent *entry;
    while ((entry = readdir(dp)) != NULL) {
        int time_num = get_time_number(entry->d_name);
        if (time_num < start_time || time_num > end_time) continue;
        char file_path[MAX_FILENAME];
        snprintf(file_path, MAX_FILENAME, "%s%s", directory_path, entry->d_name);
        FILE *fp_input = fopen(file_path, "r");
        if (!fp_input) continue;
        num_particles = 0;
        double raw_orientation;
        while (fscanf(fp_input, "%lf %lf %lf", &particles[num_particles].x, 
                     &particles[num_particles].y, &raw_orientation) == 3) 
        {
            double normalized_orientation = fmod(raw_orientation, 2.0 * PI);
            if (normalized_orientation < 0) 
            {
                normalized_orientation += 2.0 * PI;
            }
            particles[num_particles].orientation = normalized_orientation;
            num_particles++;
            if (num_particles >= MAX_PARTICLES) break;
        }
        fclose(fp_input);
        if (num_particles == 0) continue;

   
        for (int i = 0; i < num_particles; i++) {
            double ni_x = cos(particles[i].orientation);
            double ni_y = sin(particles[i].orientation);

            for (int j = 0; j < num_particles; j++) {
                if (i == j) continue;

                double dx = particles[j].x - particles[i].x;
                dx -= BoxSize_x * round(dx / BoxSize_x);
                double dy = particles[j].y - particles[i].y;
                dy -= BoxSize_y * round(dy / BoxSize_y);
                double r = sqrt(dx*dx + dy*dy);
                
                if (r > 1e-9 && r < r_max) {
                    // --- Calculate phi over the full [0, 2*PI] range ---
                    double dx_rotated = dx * ni_x + dy * ni_y;
                    double dy_rotated = -dx * ni_y + dy * ni_x;
                    double phi = atan2(dy_rotated, dx_rotated);
                    if (phi < 0) phi += 2.0 * PI;

                    // --- MODIFIED: Calculate theta over the full [0, 2*PI] range ---
                    double nj_x = cos(particles[j].orientation);
                    double nj_y = sin(particles[j].orientation);

                    // Rotate particle j's orientation vector into i's frame
                    double nj_x_rotated = nj_x * ni_x + nj_y * ni_y;
                    double nj_y_rotated = -nj_x * ni_y + nj_y * ni_x;
                    
                    double theta = atan2(nj_y_rotated, nj_x_rotated);
                    if (theta < 0) theta += 2.0 * PI;
                    // -----------------------------------------------------------

                    // Find the correct bin for this (r, phi, theta) combination
                    int r_idx = (int)(r / dr);
                    int phi_idx = (int)(phi / dphi);
                    int theta_idx = (int)(theta / dtheta);

                    if (r_idx < num_r_bins && phi_idx < num_phi_bins && theta_idx < num_theta_bins) {
                       Histogram[IDX(r_idx, phi_idx, theta_idx, num_phi_bins, num_theta_bins)] += 1.0;
                    }
                }
            }
        }
        snapshot_count++;
        printf("Processed snapshot %d: %s (%d particles)\n", snapshot_count, entry->d_name, num_particles);
    }
    closedir(dp);

    // ... (Normalization and File Writing section remains IDENTICAL) ...
    if (snapshot_count == 0) { /* handle error */ return 1; }
    FILE *output_file = fopen(output_path, "w");
    if (!output_file) { /* handle error */ return 1; }
    fprintf(output_file, "# r_center  phi_center  theta_center  g(r,phi,theta)\n");
    double area = BoxSize_x * BoxSize_y;
    double number_density = num_particles / area;
    for (int i_r = 0; i_r < num_r_bins; i_r++) {
        double r_center = (i_r + 0.5) * dr;
        double annular_area = r_center * dr * dphi; 
        for (int i_phi = 0; i_phi < num_phi_bins; i_phi++) {
            double phi_center = (i_phi + 0.5) * dphi;
            for (int i_theta = 0; i_theta < num_theta_bins; i_theta++) {
                double theta_center = (i_theta + 0.5) * dtheta;
                long current_idx = IDX(i_r, i_phi, i_theta, num_phi_bins, num_theta_bins);
                double N_l = Histogram[current_idx];
                double denominator = annular_area * num_particles * snapshot_count * number_density * dtheta;
                double g_val = 0.0;
                if (denominator > 1e-9) {
                    g_val = (2.0 * PI * N_l) / denominator;
                }
                fprintf(output_file, "%lf %lf %lf %lf\n", r_center, phi_center, theta_center, g_val);
            }
        }
    }
    fclose(output_file);
    printf("Full angle g(r, phi, theta) written to %s.\n", output_path);
    free(particles);
    free(Histogram);
    return 0;
}