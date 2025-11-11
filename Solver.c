#include "global_parameters.h"
double timer = 0;
double timer2 = 0;
#define timing(a) \
    do { \
        clock_t start = clock(); \
        a; \
        clock_t stop = clock(); \
        timer += (double)(stop - start) / CLOCKS_PER_SEC; \
    } while (0)
#define timing2(a) \
    do { \
        clock_t start = clock(); \
        a; \
        clock_t stop = clock(); \
        timer2 += (double)(stop - start) / CLOCKS_PER_SEC; \
    } while (0)

void buildM(double* arr) {

    int j = -1;
    for (int i = 0;i < n * nom;i++) {
        int index = i * n * nom;
        if (i % n == 0) {
            j++;
        }
        arr[index + i] = masses_array[j].m;
    }
}

void buildJ(double* arr) {
    int index = 0;
    for (int i = 0; i < n_phi;i++) {
        if (i == 0) {
            arr[index] = 2 * masses_array[i].x;
            arr[index + 1] = 2 * masses_array[i].y;
            arr[index + 2] = 2 * masses_array[i].z;
        }
        else {
            arr[i * nom * n + index] = -2 * (masses_array[i].x - masses_array[i - 1].x);
            arr[i * nom * n + index + 1] = -2 * (masses_array[i].y - masses_array[i - 1].y);
            arr[i * nom * n + index + 2] = -2 * (masses_array[i].z - masses_array[i - 1].z);
            arr[i * nom * n + index + 3] = 2 * (masses_array[i].x - masses_array[i - 1].x);
            arr[i * nom * n + index + 4] = 2 * (masses_array[i].y - masses_array[i - 1].y);
            arr[i * nom * n + index + 5] = 2 * (masses_array[i].z - masses_array[i - 1].z);
            index += n;
        }

    }
}

void combine_matrixA(double* A, double* M, double* J, double *JT) {
    int n_dof = nom * n;
    int total = n_dof + n_phi;
    int i, j;
    row_to_col_major(J, JT, n_phi, n_dof);
    omp_set_num_threads(6);
#pragma omp parallel for collapse(2) private(i, j)
    for (i = 0; i < total; ++i) {
        for (j = 0; j < total; ++j) {
            if (i < n_dof && j < n_dof) {
                // Top-left block: M
                A[i * total + j] = M[i * n_dof + j];
            }
            else if (i < n_dof && j >= n_dof) {
                // Top-right block: JT (already transposed)
                int col = j - n_dof;
                A[i * total + j] = JT[i * n_phi + col];
            }
            else if (i >= n_dof && j < n_dof) {
                // Bottom-left block: J
                int row = i - n_dof;
                A[i * total + j] = J[row * n_dof + j];
            }
            else {
                // Bottom-right block: zeros
                A[i * total + j] = 0.0;
            }
        }
    }
}

void buildB(double* ret) {

    for (int i = 1;i < nom * n;i += n) {
        ret[i] = -m * g;
    }

    for (int i = nom * n, j = 0;i < nom * n + n_phi;i++, j++) {
        if (i == nom * n) {
            ret[i] = -2 * (pow(masses_array[0].dx, 2) + pow(masses_array[0].dy, 2) + pow(masses_array[0].dz, 2));
        }
        else {
            ret[i] = -2 * (pow(masses_array[j].dx - masses_array[j - 1].dx, 2) + 
                           pow(masses_array[j].dy - masses_array[j - 1].dy, 2) + 
                           pow(masses_array[j].dz - masses_array[j - 1].dz, 2));
        }
    }
}


void solve(double* x, double* M, double* J, double* A, double* B, double *JT) {
     
    buildM(M);
    buildJ(J);                                // Build Matrices 
    timing(combine_matrixA(A, M, J,JT));
    buildB(B);

    MKL_INT N = n * nom + n_phi;
    cblas_dcopy(N, B, 1, x, 1);

    MKL_INT* ipiv = (MKL_INT*)malloc(N * sizeof(MKL_INT));
    if (ipiv == NULL) {
        fprintf(stderr, "Error: Failed to allocate pivot array.\n");
        return;
    }
    mkl_set_num_threads(12);
    //MKL_INT info = 
    timing2(LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, 1, A, N, ipiv, x, 1));

   /* if (info > 0) {
        fprintf(stderr, "Matrix A is singular to working precision.\n");
    }
    else if (info < 0) {
        fprintf(stderr, "Argument %d had an illegal value.\n", -info);
    }*/

    free(ipiv);
}


double* rk4(int steps, double dt, double* ddr_temp1, double* ddr_temp2, double* ddr_temp3, double* ddr_temp4, double* temp, double* results, double* M, double* J, double* A, double* B,double* JT, double* energy) {
    double* timers = malloc(2 * sizeof(double));
    int energy_index = 0;
    static double last_printed = -0.1;
    for (int step = 0; step < steps; ++step) {
        double kinetic = 0;
        double potential = 0;
        double percent = (double)step * 100.0 / steps;
        if (percent - last_printed >= 0.1) {
            printf("%.1f%%\n", percent);
            last_printed = percent;
        }

        solve(ddr_temp1, M, J, A, B, JT);                                                             //solve
        for (int j = 0;j < nom;j++) {
            int index = j * 2*n;
            kinetic += 0.5 * (pow(masses_array[j].dx, 2) + pow(masses_array[j].dy, 2));
            potential += g * masses_array[j].y;
            temp[index    ] = masses_array[j].dx;
            temp[index + 1] = masses_array[j].dy;
            temp[index + 2] = masses_array[j].dz;
            temp[index + 3] = masses_array[j].x;
            temp[index + 4] = masses_array[j].y;
            temp[index + 5] = masses_array[j].z;
            results[step * nom * n + n * j] = masses_array[j].x;                                      // Integrate...
            results[step * nom * n + n * j + 1] = masses_array[j].y;
            results[step * nom * n + n * j + 2] = masses_array[j].z;
            masses_array[j].ddx = ddr_temp1[n * j];
            masses_array[j].ddy = ddr_temp1[n * j + 1];
            masses_array[j].ddz = ddr_temp1[n * j + 2];
            masses_array[j].dx = masses_array[j].dx + masses_array[j].ddx * (dt / 2);
            masses_array[j].dy = masses_array[j].dy + masses_array[j].ddy * (dt / 2);
            masses_array[j].dz = masses_array[j].dz + masses_array[j].ddz * (dt / 2);
            masses_array[j].x = masses_array[j].x + masses_array[j].dx * (dt / 2);
            masses_array[j].y = masses_array[j].y + masses_array[j].dy * (dt / 2);
            masses_array[j].z = masses_array[j].z + masses_array[j].dz * (dt / 2);
        }

        solve(ddr_temp2, M, J, A, B, JT);
        for (int j = 0;j < nom;j++) {
            int index = j *2*n;
            masses_array[j].ddx = ddr_temp2[n * j];
            masses_array[j].ddy = ddr_temp2[n * j + 1];
            masses_array[j].ddz = ddr_temp2[n * j + 2];
            masses_array[j].dx = temp[index] + masses_array[j].ddx * (dt / 2);
            masses_array[j].dy = temp[index + 1] + masses_array[j].ddy * (dt / 2);
            masses_array[j].dz = temp[index + 2] + masses_array[j].ddz * (dt / 2);
            masses_array[j].x = temp[index + 3] + masses_array[j].dx * (dt / 2);
            masses_array[j].y = temp[index + 4] + masses_array[j].dy * (dt / 2);
            masses_array[j].z = temp[index + 5] + masses_array[j].dz * (dt / 2);
        }

        solve(ddr_temp3, M, J, A, B, JT);
        for (int j = 0;j < nom;j++) {
            int index = j * 2*n;
            masses_array[j].ddx = ddr_temp3[n * j];;
            masses_array[j].ddy = ddr_temp3[n * j + 1];
            masses_array[j].ddz = ddr_temp3[n * j + 2];
            masses_array[j].dx  = temp[index] + masses_array[j].ddx * dt;
            masses_array[j].dy  = temp[index + 1] + masses_array[j].ddy * dt;
            masses_array[j].dz  = temp[index + 2] + masses_array[j].ddz * dt;
            masses_array[j].x   = temp[index + 3] + masses_array[j].dx * dt;
            masses_array[j].y   = temp[index + 4] + masses_array[j].dy * dt;
            masses_array[j].z   = temp[index + 5] + masses_array[j].dz * dt;
        }


        solve(ddr_temp4, M, J, A, B, JT);
        for (int j = 0;j < nom;j++) {
            int index = j * 2*n;
            masses_array[j].ddx = (ddr_temp1[n * j] + 2 * ddr_temp2[n * j] + 2 * masses_array[j].ddx + ddr_temp4[n * j]) / 6;
            masses_array[j].ddy = (ddr_temp1[n * j + 1] + 2 * ddr_temp2[n * j + 1] + 2 * masses_array[j].ddy + ddr_temp4[n * j + 1]) / 6;
            masses_array[j].ddz = (ddr_temp1[n * j + 2] + 2 * ddr_temp2[n * j + 2] + 2 * masses_array[j].ddz + ddr_temp4[n * j + 2]) / 6;
            masses_array[j].dx  = temp[index    ] + masses_array[j].ddx * dt;
            masses_array[j].dy  = temp[index + 1] + masses_array[j].ddy * dt;
            masses_array[j].dz  = temp[index + 2] + masses_array[j].ddz * dt;
            masses_array[j].x   = temp[index + 3] + masses_array[j].dx * dt;
            masses_array[j].y   = temp[index + 4] + masses_array[j].dy * dt;
            masses_array[j].z   = temp[index + 5] + masses_array[j].dz * dt;
        }

        /*if (step % (steps/250) == 0) {
            energy[energy_index] = kinetic + potential;
            energy_index++;
        }*/
    }
    timers[0] = timer;
    timers[1] = timer2;
    return timers;
}