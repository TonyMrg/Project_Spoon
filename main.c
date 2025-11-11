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
const int nom = 2;
const double m = 1.0;
const int n_phi = 2;  //number of constraints
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

void visualize_double_pendulum(const double* results, int steps, int dt_ms)
{
    // Larger window
    const int WIDTH = 1200;
    const int HEIGHT = 1000;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL init failed: %s\n", SDL_GetError());
        return;
    }

    SDL_Window* window = SDL_CreateWindow(
        "Double Pendulum Visualization",
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

    // scale factor (pixels per meter)
    const double SCALE = 250.0;     // slightly larger to keep rods visible
    const int CX = WIDTH / 2;
    const int CY = HEIGHT / 4;      // pivot higher to fit long swings

    int running = 1;
    SDL_Event e;

    for (int step = 0; step < steps && running; ++step) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
        }

        double x1 = results[step * 4 + 0];
        double y1 = results[step * 4 + 1];
        double x2 = results[step * 4 + 2];
        double y2 = results[step * 4 + 3];

        int px0 = CX;
        int py0 = CY;
        int px1 = CX + (int)(SCALE * x1);
        int py1 = CY - (int)(SCALE * y1);
        int px2 = CX + (int)(SCALE * x2);
        int py2 = CY - (int)(SCALE * y2);

        // Clear screen
        SDL_SetRenderDrawColor(renderer, 20, 20, 30, 255);
        SDL_RenderClear(renderer);

        // Draw rods
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderDrawLine(renderer, px0, py0, px1, py1);
        SDL_RenderDrawLine(renderer, px1, py1, px2, py2);

        // Draw masses
        SDL_Rect m1 = { px1 - 5, py1 - 5, 10, 10 };
        SDL_Rect m2 = { px2 - 5, py2 - 5, 10, 10 };
        SDL_SetRenderDrawColor(renderer, 255, 80, 80, 255);
        SDL_RenderFillRect(renderer, &m1);
        SDL_SetRenderDrawColor(renderer, 80, 160, 255, 255);
        SDL_RenderFillRect(renderer, &m2);

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

    double* arr = (double*)calloc(n_phi * n * nom, sizeof(double));
    if (!arr) {
        fprintf(stderr, "Error: calloc failed in BuildJ\n");
        return NULL;
    }

    for (int i = 0; i < n_phi;i++) {
        int index = i*nom*n;
        if (i == 0) {
            arr[index] = 2 * masses_array[i].x;
            arr[index + 1] = 2 * masses_array[i].y;
        }
        else {
            arr[index] = -2 * (masses_array[i].x - masses_array[i-1].x);
            arr[index + 1] = -2 * (masses_array[i].y - masses_array[i - 1].y);
            arr[index + 2] = 2 * (masses_array[i].x - masses_array[i - 1].x);
            arr[index + 3] = 2 * (masses_array[i].y - masses_array[i - 1].y);
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

double* truncateResults(double* results, int steps, int nom, int* newLength) {
    int total = steps * nom * 2;
    double* truncated = malloc(total * sizeof(double));
    if (!truncated) {
        perror("malloc failed");
        return NULL;
    }

    int j = 0; // index in truncated array

    for (int i = 0; i < total; i++) {
        double val = results[i];
        if (val > 2.0 || val < -2.0) {
            // If y-coordinate, remove also previous x
            if (i % 2 == 1) {
                if (j > 0) j--; // remove previous x
            }
            break; // stop copying
        }
        truncated[j++] = val;
    }

    *newLength = j; // return new length
    return truncated;
}

int main(void) {
    masses_array = malloc(nom * sizeof(struct mass));
    if (!masses_array) {
        fprintf(stderr, "Error: malloc failed in masses array\n");
        return 1;
    }

    for (int i = 0;i < 2;i++) {
        struct mass lala;
        lala.m = 1;
        lala.x = i + 1;
        lala.y = 0;
        lala.dx = 0;
        lala.dy = 0;
        masses_array[i] = lala;
    }
    
    double d1 = sqrt(masses_array[0].x * masses_array[0].x + masses_array[0].y * masses_array[0].y);
    double d2 = sqrt((masses_array[1].x - masses_array[0].x) * (masses_array[1].x - masses_array[0].x) + (masses_array[1].y - masses_array[0].y) * (masses_array[1].y - masses_array[0].y));
    if (fabs(d1 - l) > 1e-12 || fabs(d2 - l) > 1e-12) {
        printf(stderr, "Constraint not satisfied at init!\n");
    }
    
    
    
    int steps = 1000;
    double* results = malloc(steps * nom * 2 * sizeof(double));;
    if (!results) {
        fprintf(stderr, "Error: malloc failed in results array\n");
        return 1;
    }
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
                fprintf(stderr, "NaN or Inf detected in q_ddot[%d]!\n", k);
                return 1;
            
            }
        }
        for (int j = 0;j < nom;j++) {
            results[step * nom * 2 + 2*j] = masses_array[j].x;
            results[step * nom * 2 + 2*j + 1] = masses_array[j].y;
            masses_array[j].ddx = gsl_vector_get(ddr, 2*j);
            masses_array[j].ddy = gsl_vector_get(ddr, 2*j + 1);
            masses_array[j].dx = masses_array[j].dx + masses_array[j].ddx * 0.01;
            masses_array[j].dy = masses_array[j].dy + masses_array[j].ddy * 0.01;
            masses_array[j].x = masses_array[j].x + masses_array[j].dx * 0.01;
            masses_array[j].y = masses_array[j].y + masses_array[j].dy * 0.01;

        }
        free(J);
        free(M);
        free(A);
        free(B);
        gsl_vector_free(ddr);
        gsl_permutation_free(p);
        
    }
    
    
    
    visualize_double_pendulum(results, steps, 16);
    return 0;
}
