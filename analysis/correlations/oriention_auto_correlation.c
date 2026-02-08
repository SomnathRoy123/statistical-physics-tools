/**
 * @file orientation_autocorrelation.c
 * @brief Calculates the orientational autocorrelation function from particle data files.
 *
 * This program scans a directory for time-stamped data files, reads the orientation
 * (3rd column) for each particle, and computes the orientational autocorrelation
 * function.
 *
 * The formula used is C(t) = <cos(theta(t) - theta(0))>, where the average is
 * taken over all particles. This method is suitable for unwrapped angular data.
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <dirent.h>
 #include <math.h>
 
 // If M_PI is not defined in math.h, define it.
 #ifndef M_PI
 #define M_PI 3.14159265358979323846
 #endif
 
 // --- Configuration ---
 const char* DIRECTORY_PATH = "/home/somnath2/Codes/Trial_70/TA160_R3.5/evo/";
 const char* OUTPUT_FILE_PATH = "/home/somnath2/Codes/polar_order/orient_time_corr/Trial_102/orientation_autocorrelation_160_3.5.dat";
 const int TIME_MIN = 100;
 const int TIME_MAX = 3000;
 
 // --- Data Structures ---
 
 /**
  * @struct FileInfo
  * @brief Stores the path and extracted time number for a data file.
  */
 typedef struct {
     int time_num;
     char path[1024];
 } FileInfo;
 
 // --- Function Prototypes ---
 
 /**
  * @brief Comparison function for qsort to sort files by their time number.
  * @param a Pointer to first FileInfo element
  * @param b Pointer to second FileInfo element
  * @return Integer indicating relative order
  */
 int compare_files(const void* a, const void* b);
 
 /**
  * @brief Counts the number of lines (particles) in a given data file.
  * @param filepath Path to the data file
  * @return Number of particles, or -1 on error
  */
 int count_particles_in_file(const char* filepath);
 
 /**
  * @brief Reads the 3rd column (orientation) from a file into an array.
  * @param filepath Path to the data file
  * @param orientations Array to store orientation values
  * @param num_particles Expected number of particles
  * @return 0 on success, -1 on error
  */
 int read_orientations(const char* filepath, double* orientations, int num_particles);
 
 // --- Main Execution ---
 
 int main() {
     printf("Scanning directory for data files: %s\n", DIRECTORY_PATH);
 
     DIR* directory = opendir(DIRECTORY_PATH);
     if (!directory) {
         fprintf(stderr, "Error opening directory '%s'\n", DIRECTORY_PATH);
         return 1;
     }
 
     // --- 1. Find and sort all relevant data files ---
     struct dirent* dir_entry;
     FileInfo* file_list = NULL;
     int file_count = 0;
     int file_capacity = 0;
 
     while ((dir_entry = readdir(directory)) != NULL) {
         int time_num;
         
         // Check if filename matches pattern "time_*.dat"
         if (sscanf(dir_entry->d_name, "time_%d.dat", &time_num) == 1) {
             
             // Check if time is within specified range
            //  if (time_num >= TIME_MIN && time_num <= TIME_MAX && time_num%2==0) {     // making the change here
             if (time_num >= TIME_MIN && time_num <= TIME_MAX ) {     // making the change here

                 // Expand array if needed
                 if (file_count >= file_capacity) {
                     file_capacity = (file_capacity == 0) ? 16 : file_capacity * 2;
                     FileInfo* new_list = realloc(file_list, file_capacity * sizeof(FileInfo));
                     
                     if (!new_list) {
                         fprintf(stderr, "Memory allocation failed for file list.\n");
                         free(file_list);
                         closedir(directory);
                         return 1;
                     }
                     file_list = new_list;
                 }
                 
                 // Store file information
                //  file_list[file_count].time_num = time_num/2;                 // changes made here
                 file_list[file_count].time_num = time_num;                 // changes made here

                 snprintf(file_list[file_count].path, sizeof(file_list[file_count].path), 
                         "%s%s", DIRECTORY_PATH, dir_entry->d_name);
                 file_count++;
             }
         }
     }
     closedir(directory);
 
     if (file_count < 1) {
         fprintf(stderr, "No data files found in the specified time range.\n");
         free(file_list);
         return 1;
     }
 
     // Sort files by time number
     qsort(file_list, file_count, sizeof(FileInfo), compare_files);
     printf("Found and sorted %d files.\n", file_count);
 
     // --- 2. Setup memory based on particle count ---
     int num_particles = count_particles_in_file(file_list[0].path);
     if (num_particles <= 0) {
         fprintf(stderr, "Could not determine particle count from first file: %s\n", 
                 file_list[0].path);
         free(file_list);
         return 1;
     }
     printf("Detected %d particles per snapshot.\n", num_particles);
 
     // Allocate memory for orientation data
     double* initial_orientations = malloc(num_particles * sizeof(double));
     double* current_orientations = malloc(num_particles * sizeof(double));
     double* autocorrelation_results = malloc(file_count * sizeof(double));
 
     if (!initial_orientations || !current_orientations || !autocorrelation_results) {
         fprintf(stderr, "Memory allocation failed for orientation arrays.\n");
         free(file_list);
         free(initial_orientations);
         free(current_orientations);
         free(autocorrelation_results);
         return 1;
     }
 
     // --- 3. Read initial orientations (t=0) ---
     printf("Reading initial orientations from: %s\n", file_list[0].path);
     if (read_orientations(file_list[0].path, initial_orientations, num_particles) != 0) {
         fprintf(stderr, "Error reading initial orientation file. Aborting.\n");
         // Cleanup and exit
         free(file_list);
         free(initial_orientations);
         free(current_orientations);
         free(autocorrelation_results);
         return 1;
     }
 
     // --- 4. Calculate Autocorrelation for each time step ---
     printf("Calculating autocorrelation...\n");
     for (int t = 0; t < file_count; ++t) {
         
         // Read the orientations for the current time step
         if (read_orientations(file_list[t].path, current_orientations, num_particles) != 0) {
             fprintf(stderr, "Error reading file '%s', stopping calculation.\n", 
                     file_list[t].path);
             break; // Stop if a file is unreadable
         }
 
         double sum_cos_diff = 0.0;
         
         // Average over all particles for this time step
         for (int p = 0; p < num_particles; ++p) {
             double delta_theta = current_orientations[p] - initial_orientations[p];
             sum_cos_diff += cos(delta_theta);
         }
         
         autocorrelation_results[t] = sum_cos_diff / (double)num_particles;
     }
 
     // --- 5. Save results to output file ---
     FILE* output_file = fopen(OUTPUT_FILE_PATH, "w");
     if (!output_file) {
         fprintf(stderr, "Could not open output file: %s\n", OUTPUT_FILE_PATH);
     } else {
         fprintf(output_file, "# Time_Step\tOrientational_Autocorrelation\n");
         
         for (int t = 0; t < file_count; ++t) {
             // Use the time number from the file for robustness
             fprintf(output_file, "%d\t%.8f\n", 
                     file_list[t].time_num, autocorrelation_results[t]);
         }
         
         fclose(output_file);
         printf("Processing complete. Results saved to %s\n", OUTPUT_FILE_PATH);
     }
 
     // --- 6. Cleanup ---
     free(initial_orientations);
     free(current_orientations);
     free(autocorrelation_results);
     free(file_list);
 
     return 0;
 }
 
 // --- Function Implementations ---
 
 /**
  * @brief Comparison function for qsort to sort files by their time number.
  */
 int compare_files(const void* a, const void* b) {
     const FileInfo* file_a = (const FileInfo*)a;
     const FileInfo* file_b = (const FileInfo*)b;
     return (file_a->time_num - file_b->time_num);
 }
 
 /**
  * @brief Counts the number of lines (particles) in a given data file.
  */
 int count_particles_in_file(const char* filepath) {
     FILE* file = fopen(filepath, "r");
     if (!file) {
         perror("Error opening file for counting particles");
         return -1;
     }
     
     int count = 0;
     char line[256];
     
     while (fgets(line, sizeof(line), file)) {
         // Basic check to not count empty lines
         if (strspn(line, " \t\n\r") != strlen(line)) {
            count++;
         }
     }
     
     fclose(file);
     return count;
 }
 
 /**
  * @brief Reads the 3rd column (orientation) from a file into an array.
  */
 int read_orientations(const char* filepath, double* orientations, int num_particles) {
     FILE* file = fopen(filepath, "r");
     if (!file) {
         fprintf(stderr, "Error opening file '%s'\n", filepath);
         return -1;
     }
     
     for (int i = 0; i < num_particles; ++i) {
         double col1, col2, col3;
         
         // Try to read three columns
         if (fscanf(file, "%lf %lf %lf", &col1, &col2, &col3) != 3) {
             // Check if we've reached end of file unexpectedly
             if (feof(file)) {
                 fprintf(stderr, "Error: Unexpected end of file in '%s'. "
                                "Expected %d particles, found %d.\n", 
                         filepath, num_particles, i);
                 fclose(file);
                 return -1;
             }
             
             // Other read error
             fprintf(stderr, "Error reading particle %d from file '%s'\n", i, filepath);
             fclose(file);
             return -1;
         }
         
         orientations[i] = col3;  // Store the 3rd column (orientation)
     }
     
     fclose(file);
     return 0;
 }