#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <SDL.h>
#include <stdarg.h>


// Simulation parameters
const double g = 9.81;
const int nom = 20;
const double m = 1.0;
const int n_phi = 20;  //number of constraints
const int n = 2;
const double l = 1.0;

double* L;

struct mass* masses_array;

double* buildB();
double* buildJ();
double* combine_matrixA(double* M, double* J);
double* buildM();
void row_to_col_major(double* src, double* dst, int rows, int cols);

struct mass {
    double x;
    double y;
    double dx;
    double dy;
    double ddx;
    double ddy;
    double m;
};

void visualize_double_pendulum(const double* results, int steps, int dt_ms, int n_masses)
{
    const int WIDTH = 1920;
    const int HEIGHT = 1400;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL init failed: %s\n", SDL_GetError());
        return;
    }

    SDL_Window* window = SDL_CreateWindow(
        "Pendulum Visualization",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIDTH, HEIGHT, SDL_WINDOW_SHOWN);

    if (!window) {
        fprintf(stderr, "Window creation failed: %s\n", SDL_GetError());
        SDL_Quit();
        return;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        fprintf(stderr, "Renderer creation failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return;
    }

    const double SCALE = 50;    // pixels per meter
    const int CX = WIDTH / 2;
    const int CY = HEIGHT / 4;

    int running = 1;
    SDL_Event e;

    // optional colors for masses (cycled if more than 6)
    SDL_Color mass_colors[] = {
        {255, 80, 80, 255}, {80, 160, 255, 255}, {80, 255, 120, 255},
        {255, 200, 50, 255}, {180, 100, 255, 255}, {255, 100, 180, 255}
    };
    int n_colors = sizeof(mass_colors) / sizeof(SDL_Color);

    for (int step = 0; step < steps && running; ++step) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
        }

        SDL_SetRenderDrawColor(renderer, 20, 20, 30, 255);
        SDL_RenderClear(renderer);

        // draw rods between consecutive masses
        int prev_px = CX, prev_py = CY;
        for (int i = 0; i < n_masses; ++i) {
            double x = results[step * n_masses * 2 + 2 * i];
            double y = results[step * n_masses * 2 + 2 * i + 1];

            int px = CX + (int)(SCALE * x);
            int py = CY - (int)(SCALE * y); // flip y

            // draw rod from previous mass (or pivot) to current mass
            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
            SDL_RenderDrawLine(renderer, prev_px, prev_py, px, py);

            // draw current mass
            SDL_Rect mrect = { px - 5, py - 5, 10, 10 };
            SDL_Color c = mass_colors[i % n_colors];
            SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, c.a);
            SDL_RenderFillRect(renderer, &mrect);

            prev_px = px;
            prev_py = py;
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(dt_ms);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void row_to_col_major(double* src, double* dst, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            dst[j * rows + i] = src[i * cols + j];
        }
    }
}

double* buildM() {
    double* arr = calloc((size_t)(nom * n) * (size_t)(nom * n) , sizeof(double));
    if (!arr) {
        perror("malloc failed in BuildM()");
        return NULL;
    }
   
    int j = -1;
    for (int i = 0;i < n * nom;i++) {
        int index = i * n * nom;
        if (i % 2 == 0) {
            j++;
        }
        arr[index + i] = masses_array[j].m;
    }

    
    return arr;
}

double* buildJ() {

    double* arr = (double*)calloc(n_phi * n * nom*2, sizeof(double));
    if (!arr) {
        fprintf(stderr, "Error: calloc failed in BuildJ\n");
        return NULL;
    }
    int index = 0;
    for (int i = 0; i < n_phi;i++) {
        if (i == 0) {
            arr[index] = 2 * masses_array[i].x;
            arr[index + 1] = 2 * masses_array[i].y;
        }
        else {
            arr[i * nom * n + index] = -2 * (masses_array[i].x - masses_array[i - 1].x);
            arr[i * nom * n + index + 1] = -2 * (masses_array[i].y - masses_array[i - 1].y);
            arr[i * nom * n + index + 2] = 2 * (masses_array[i].x - masses_array[i - 1].x);
            arr[i * nom * n + index + 3] = 2 * (masses_array[i].y - masses_array[i - 1].y);
            index += 2;
            
        }
        
    }
    //row_to_col_major(arr, arr_clmn, 2, 4);
    
    return arr;
}

double* combine_matrixA(double* M, double* J) {
    int n_dof = nom * n;                       // number of dynamic DOFs
    int total = n_dof + n_phi;                 // total matrix size
    double* A = calloc(total * total, sizeof(double));
    if (!A) {
        fprintf(stderr, "Error: calloc failed in combine_matrixA()\n");
        return NULL;
    }
    if (!A) return NULL;

    // Upper-left block: M (size n_dof x n_dof)
    for (int i = 0; i < n_dof; ++i)
        for (int j = 0; j < n_dof; ++j)
            A[i * total + j] = M[i * n_dof + j];

    // Upper-right block: J^T (size n_dof x n_phi)
    double* JT = malloc(n_dof * n_phi * sizeof(double));
    if (!JT) {
        fprintf(stderr, "Error: JT calloc failed in combine_matrixA()\n");
        return NULL;
    }
    row_to_col_major(J, JT, n_phi, n_dof);
    for (int i = 0; i < n_dof; ++i)
        for (int j = 0; j < n_phi; ++j)
            A[i * total + (n_dof + j)] = JT[i * n_phi + j];
    free(JT);

    // Lower-left block: J (size n_phi x n_dof)
    for (int i = 0; i < n_phi; ++i)
        for (int j = 0; j < n_dof; ++j)
            A[(n_dof + i) * total + j] = J[i * n_dof + j];

    // Lower-right block: zeros (already zeroed by calloc)
    return A;
}

double* buildB(void) {
    double* ret = calloc(nom * n + n_phi, sizeof(double));
    if (!ret) {
        fprintf(stderr, "Error: calloc failed in buildB()\n");
        return NULL;
    }
    for (int i = 1;i < nom * n;i += 2) {
        ret[i] = -m * g;
    }
    
    for (int i = nom * n,j=0;i < nom * n + n_phi;i++,j++) {
        if (i == nom * n) {
            ret[i] = -2 * (pow(masses_array[0].dx, 2)+ pow(masses_array[0].dy, 2));
        }
        else {
            ret[i] = -2 * (pow(masses_array[j].dx - masses_array[j-1].dx, 2) + pow(masses_array[j].dy - masses_array[j-1].dy, 2));
        }
    }
    return ret;
}

void save_energy(double* energy, int n) {
    FILE* f = fopen("C://Users//kosta//Downloads//energy.csv", "w");
    if (!f) return;
    for (int i = 0; i < n; i++) {
        fprintf(f, "%d,%f\n", i, energy[i]);
    }
    fclose(f);
}


int main(void) {
    masses_array = malloc(nom * sizeof(struct mass));
    if (!masses_array) {
        fprintf(stderr, "Error: malloc failed in masses array\n");
        return 1;
    }

    for (int i = 0;i < nom;i++) {
        struct mass lala;
        lala.m = 1;
        lala.x = i + 1;
        lala.y = 0;
        lala.dx = 0;
        lala.dy = 0;
        masses_array[i] = lala;
    }
    
    int failed_sum = 0;
    int steps = 10000;
    double* results = malloc(steps * nom * 2 * sizeof(double));
    double* energy = malloc(250*sizeof(double));
    if (!results) {
        fprintf(stderr, "Error: malloc failed in results array\n");
        return 1;
    }
    if (!energy) {
        fprintf(stderr, "Error: malloc failed in energy array\n");
        return 1;
    }
    int energy_index = 0;
    for (int step = 0; step < steps; ++step) {
        double* J = buildJ();
        double* M = buildM();
        double* A = combine_matrixA(M, J);
        double* B = buildB();

        gsl_matrix_view A_view = gsl_matrix_view_array(A, n * nom + n_phi, n * nom + n_phi);
        gsl_vector_view B_view = gsl_vector_view_array(B, n * nom + n_phi);
        gsl_vector* ddr = gsl_vector_alloc(n * nom + n_phi);
        int signum;
        gsl_permutation* p = gsl_permutation_alloc((n * nom + n_phi));
        gsl_linalg_LU_decomp(&A_view.matrix, p, &signum);
        gsl_linalg_LU_solve(&A_view.matrix, p, &B_view.vector, ddr);
        
        for (int k = 0; k < n * nom + n_phi; ++k) {
            double val = gsl_vector_get(ddr, k);
            if (!isfinite(val)) {
                failed_sum += 1;
                //fprintf(stderr, "NaN or Inf detected in ddr[%d]        %d!\n", k,step);            
            }
        }
        double* temp = malloc(4 *nom* sizeof(double));
        if (!temp) {
            fprintf(stderr, "Error: malloc failed for temp array\n");
            return 1;
        }

        double kinetic=0;
        double potential=0;
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            kinetic += 0.5 * (pow(masses_array[j].dx, 2)+pow(masses_array[j].dy,2));
            potential += g * masses_array[j].y;
            results[step * nom * 2 + 2*j] = masses_array[j].x;
            results[step * nom * 2 + 2*j + 1] = masses_array[j].y;
            masses_array[j].ddx = gsl_vector_get(ddr, 2*j);
            masses_array[j].ddy = gsl_vector_get(ddr, 2*j + 1);
            temp[index] = masses_array[j].dx;
            temp[index + 1] = masses_array[j].dy;
            temp[index + 2] = masses_array[j].x;
            temp[index + 3] = masses_array[j].y;
            masses_array[j].dx = masses_array[j].dx + masses_array[j].ddx * 0.001;
            masses_array[j].dy = masses_array[j].dy + masses_array[j].ddy * 0.001;
            masses_array[j].x = masses_array[j].x + masses_array[j].dx * 0.001;
            masses_array[j].y = masses_array[j].y + masses_array[j].dy * 0.001;
        }
        
        double* J_temp = buildJ();
        double* M_temp= buildM();
        double* A_temp = combine_matrixA(M_temp, J_temp);
        double* B_temp = buildB();
        gsl_matrix_view A_temp_view1 = gsl_matrix_view_array(A_temp, n * nom + n_phi, n * nom + n_phi);
        gsl_vector_view B_temp_view1 = gsl_vector_view_array(B_temp, n * nom + n_phi);
        gsl_vector* ddr_temp = gsl_vector_alloc(n * nom + n_phi);
        int signum_temp;
        gsl_permutation* p_temp = gsl_permutation_alloc((n * nom + n_phi));
        gsl_linalg_LU_decomp(&A_temp_view1.matrix, p_temp, &signum_temp);
        gsl_linalg_LU_solve(&A_temp_view1.matrix, p_temp, &B_temp_view1.vector, ddr_temp);

        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            masses_array[j].ddx = (gsl_vector_get(ddr_temp, 2 * j) + masses_array[j].ddx)/2;
            masses_array[j].ddy = (gsl_vector_get(ddr_temp, 2 * j + 1)+ masses_array[j].ddy)/2;
            masses_array[j].dx = temp[index] + masses_array[j].ddx * 0.001;
            masses_array[j].dy = temp[index+1] + masses_array[j].ddy * 0.001;
            masses_array[j].x = temp[index+2] + masses_array[j].dx * 0.001;
            masses_array[j].y = temp[index+3] + masses_array[j].dy * 0.001;
        }

        if (step % 40 == 0) {
            energy[energy_index] = kinetic + potential;
            energy_index++;

        }
        free(temp);
        free(J_temp);
        free(M_temp);
        free(A_temp);
        free(B_temp);
        gsl_vector_free(ddr_temp);
        gsl_permutation_free(p_temp);
        free(J);
        free(M);
        free(A);
        free(B);
        gsl_vector_free(ddr);
        gsl_permutation_free(p);
    }

    save_energy(energy, 250);
    printf("NaN or Inf: %d", failed_sum);
    visualize_double_pendulum(results, steps, 1,nom);
    return 0;
}
