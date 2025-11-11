#include "global_parameters.h"

double g = 9.81;
int nom = 50;
double m = 1.0;
int n_phi = 50;
int n = 2;
double l = 1.0;
struct mass* masses_array;

void row_to_col_major(double* src, double* dst, int rows, int cols);
void save_energy(double* energy, int n);
void buildB(double* ret);
void buildJ(double* arr);
void combine_matrixA(double* A, double* M, double* J);
void buildM(double* arr);
void visualize_double_pendulum(const double* results, int steps, int dt_ms, double playback_time_s, int n_masses);
void rk4(int steps, double dt, gsl_vector* ddr_temp1, gsl_vector* ddr_temp2, gsl_vector* ddr_temp3, gsl_vector* ddr_temp4, double* temp, double* results, double* M, double* J, double* A, double* B, gsl_permutation* p, double* energy);
void write_results_txt(double* results, int steps, double dt);
double* read_results_txt(int* nom, int* steps, double* dt);

int main(void) {
    char mode[30];
    char txt_file_path[200];
    printf("Enter mode: \"Solve\" or \"Vis\": ");
    fgets(mode, sizeof(mode), stdin);
    mode[strcspn(mode, "\n")] = '\0';

    if (strcmp(mode, "Solve")==0 || strcmp(mode, "solve") == 0) {
        clock_t start, end;
        double cpu_time_used;
        start = clock();

        int steps = 500000;
        double dt = 0.0002;
        masses_array = malloc(nom * sizeof(struct mass)); if (!masses_array) { return 1; }
        for (int i = 0;i < nom;i++) {
            struct mass lala;
            lala.m = 1;
            lala.x = i + 1;
            lala.y = 0;
            lala.dx = 0;
            lala.dy = 0;
            masses_array[i] = lala;
        }
        double* results = malloc(steps * nom * 2 * sizeof(double));if (!results) { return 1; }
        double* energy = malloc(250 * sizeof(double)); if (!energy) { return 1; }
        double* temp = malloc(4 * nom * sizeof(double)); if (!temp) { return 1; }
        gsl_vector* ddr_temp1 = gsl_vector_alloc(n * nom + n_phi);
        gsl_vector* ddr_temp2 = gsl_vector_alloc(n * nom + n_phi);
        gsl_vector* ddr_temp3 = gsl_vector_alloc(n * nom + n_phi);
        gsl_vector* ddr_temp4 = gsl_vector_alloc(n * nom + n_phi);
        double* M = calloc((size_t)(nom * n) * (size_t)(nom * n), sizeof(double)); if (!M) { return 1; }
        double* B = calloc(nom * n + n_phi, sizeof(double)); if (!B) { return 1; }
        double* J = (double*)calloc(n_phi * n * nom * 2, sizeof(double)); if (!J) { return 1; }
        int n_dof = nom * n;
        int total = n_dof + n_phi;
        double* A = calloc((nom * n + n_phi) * (nom * n + n_phi), sizeof(double)); if (!A) { return 1; }
        gsl_permutation* p = gsl_permutation_alloc((n * nom + n_phi));

        rk4(steps, dt, ddr_temp1, ddr_temp2, ddr_temp3, ddr_temp4, temp, results, M, J, A, B, p, energy);

        save_energy(energy, 250);
        free(temp);
        free(B);
        free(M);
        free(J);
        free(A);
        free(energy);                                                              
        free(masses_array);
        gsl_vector_free(ddr_temp1);
        gsl_vector_free(ddr_temp2);
        gsl_vector_free(ddr_temp3);
        gsl_vector_free(ddr_temp4);
        gsl_permutation_free(p);
        end = clock();

        cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Time taken: %f seconds\n", cpu_time_used);
        
        char io_bool[20];
        printf("Want to write the results? (Y/N): ");
        fgets(io_bool, sizeof(io_bool), stdin);
        io_bool[strcspn(io_bool, "\n")] = '\0';
        if (strcmp(io_bool, "Y") == 0 || strcmp(io_bool, "y") == 0) {
            write_results_txt(results, steps, dt);
        }
        printf("Want to visualize the results? (Y/N): ");
        fgets(io_bool, sizeof(io_bool), stdin);
        io_bool[strcspn(io_bool, "\n")] = '\0';
        if (strcmp( io_bool, "Y") == 0 || strcmp(io_bool, "y") == 0) {
            double play_back_time;
            printf("Normal playback time: %lf seconds. Input: ", steps * dt);
            fscanf_s(stdin, "%lf", &play_back_time);
            visualize_double_pendulum(results, steps, dt, play_back_time, nom);
        }
        free(results);
    }
    else if (strcmp(mode, "Vis") == 0 || strcmp(mode, "vis") == 0) {   
        int steps;
        double dt;
        double* results= read_results_txt(&nom, &steps, &dt);
        double play_back_time;
        printf("Normal playback time: %lf seconds. Input: ", steps * dt);
        fscanf_s(stdin, "%lf", &play_back_time);
        visualize_double_pendulum(results, steps, dt, play_back_time, nom);

    }



    
    return 0;
}