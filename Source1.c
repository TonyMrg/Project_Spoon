#include "global_parameters.h"

double g = 9.81;
int nom = 50;
double m = 1.0;
int n_phi = 50;
int n = 3;
double l = 1.0;
struct mass* masses_array;
double matrix_building_time;

void row_to_col_major(double* src, double* dst, int rows, int cols);
void save_energy(double* energy, int n);
void buildB(double* ret);
void buildJ(double* arr);
void combine_matrixA(double* A, double* M, double* J, double* JT);
void buildM(double* arr);
void visualize_double_pendulum(const double* results, int steps, int dt_ms, double playback_time_s, int n_masses);
double* rk4(int steps, double dt, 
    double* ddr_temp1, double* ddr_temp2, 
    double* ddr_temp3, double* ddr_temp4, 
    double* temp, double* results, 
    double* M, double* J, double* A, 
    double* B,double* JT, double* energy);
void write_results_txt(double* results, int steps, double dt);
double* read_results_txt(int* nom, int* steps, double* dt);
void solve(double* x, double* M, double* J, double* A, double* B, double *JT);

int main(void) {
    int n_dof = nom * n;
    int total = n_dof + n_phi;
    char mode[30];
    char txt_file_path[200];
    printf("Enter mode: \"Solve\" or \"Vis\": ");
    fgets(mode, sizeof(mode), stdin);
    mode[strcspn(mode, "\n")] = '\0';

    if (strcmp(mode, "Solve") == 0 || strcmp(mode, "solve") == 0) {

        int steps = 5000;
        double dt = 0.01;
        masses_array = malloc(nom * sizeof(struct mass)); if (!masses_array) { return 1; }
        for (int i = 0;i < nom;i++) {
            struct mass lala;
            lala.m = 1;
            lala.x = i * sqrt(1) + sqrt(1);                                     //Initial conditions
            lala.y = 0;
            lala.z = i * sqrt(1) + sqrt(1);
            lala.dx = 0;
            lala.dy = 0;
            lala.dz = 0;
            masses_array[i] = lala;
        }
        double* results = malloc(steps * nom * n * sizeof(double));if (!results) { return 1; }             // [x,y]
        double* energy = malloc(250 * sizeof(double)); if (!energy) { return 1; }
        double* temp = malloc(n*n * nom * sizeof(double)); if (!temp) { return 1; }
        double* ddr_temp1 = calloc(n * nom + n_phi, sizeof(double));                                       // ddr_temp -> Each step at Runge Kutta 4
        double* ddr_temp2 = calloc(n * nom + n_phi, sizeof(double));                                       //
        double* ddr_temp3 = calloc(n * nom + n_phi, sizeof(double));                                       //
        double* ddr_temp4 = calloc(n * nom + n_phi, sizeof(double));                                       //
        double* M = calloc((size_t)(nom * n) * (size_t)(nom * n), sizeof(double)); if (!M) { return 1; }   // Mass matrix
        double* B = calloc(nom * n + n_phi, sizeof(double)); if (!B) { return 1; }                         // RHS of Ax=B
        double* J = (double*)calloc(n_phi * n * nom , sizeof(double)); if (!J) { return 1; }               // Jacobian
        double* JT = malloc(n_dof * n_phi * sizeof(double)); if (!JT) { return; }                          // transposed
        double* A = calloc((nom * n + n_phi) * (nom * n + n_phi), sizeof(double)); if (!A) { return 1; }   // LHS of Ax=B
   
        double* timer_ = rk4(steps, dt, ddr_temp1, ddr_temp2, ddr_temp3, ddr_temp4, temp, results, M, J, A, B,JT, energy);

        //save_energy(energy, 250);
        free(temp);
        free(B);
        free(M);
        free(J);
        free(A);
        free(energy);
        free(masses_array);
        free(ddr_temp1);
        free(ddr_temp2);
        free(ddr_temp3);
        free(ddr_temp4);
        free(JT);

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
        if (strcmp(io_bool, "Y") == 0 || strcmp(io_bool, "y") == 0) {
            double play_back_time;
            printf("Normal playback time: %lf seconds. Input: ", steps * dt);
            fscanf_s(stdin, "%lf", &play_back_time);
            visualize_double_pendulum(results, steps, dt, play_back_time, nom);
        }
        free(results);
        printf("CombineA() took in total: %lf seconds", timer_[0]);
        printf("Solve() took in total: %lf seconds", timer_[1]);
    }
    else if (strcmp(mode, "Vis") == 0 || strcmp(mode, "vis") == 0) {
        int steps;
        double dt;
        double* results = read_results_txt(&nom, &steps, &dt);
        double play_back_time;
        printf("Normal playback time: %lf seconds. Input: ", steps * dt);
        fscanf_s(stdin, "%lf", &play_back_time);
        visualize_double_pendulum(results, steps, dt, play_back_time, nom);

    } 
    return 0;
}