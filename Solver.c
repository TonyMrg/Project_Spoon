#include "global_parameters.h"

void buildM(double* arr) {

    int j = -1;
    for (int i = 0;i < n * nom;i++) {
        int index = i * n * nom;
        if (i % 2 == 0) {
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
        }
        else {
            arr[i * nom * n + index] = -2 * (masses_array[i].x - masses_array[i - 1].x);
            arr[i * nom * n + index + 1] = -2 * (masses_array[i].y - masses_array[i - 1].y);
            arr[i * nom * n + index + 2] = 2 * (masses_array[i].x - masses_array[i - 1].x);
            arr[i * nom * n + index + 3] = 2 * (masses_array[i].y - masses_array[i - 1].y);
            index += 2;
        }

    }
}

void combine_matrixA(double* A, double* M, double* J) {
    int n_dof = nom * n;
    int total = n_dof + n_phi;
    for (int i = 0; i < total * total; i++) A[i] = 0.0;    //zero A

    for (int i = 0; i < n_dof; ++i)
        for (int j = 0; j < n_dof; ++j)
            A[i * total + j] = M[i * n_dof + j];

    double* JT = malloc(n_dof * n_phi * sizeof(double)); if (!JT) { return;}
    
    row_to_col_major(J, JT, n_phi, n_dof);
    for (int i = 0; i < n_dof; ++i)
        for (int j = 0; j < n_phi; ++j)
            A[i * total + (n_dof + j)] = JT[i * n_phi + j];
    free(JT);

    for (int i = 0; i < n_phi; ++i)
        for (int j = 0; j < n_dof; ++j)
            A[(n_dof + i) * total + j] = J[i * n_dof + j];

}

void buildB(double* ret) {

    for (int i = 1;i < nom * n;i += 2) {
        ret[i] = -m * g;
    }

    for (int i = nom * n, j = 0;i < nom * n + n_phi;i++, j++) {
        if (i == nom * n) {
            ret[i] = -2 * (pow(masses_array[0].dx, 2) + pow(masses_array[0].dy, 2));
        }
        else {
            ret[i] = -2 * (pow(masses_array[j].dx - masses_array[j - 1].dx, 2) + pow(masses_array[j].dy - masses_array[j - 1].dy, 2));
        }
    }
}

void solve(double* x, double* M, double* J, double* A, double* B, gsl_permutation* p) {
    buildM(M);
    buildJ(J);
    combine_matrixA(A, M, J);
    buildB(B);
    gsl_matrix_view A_view = gsl_matrix_view_array(A, n * nom + n_phi, n * nom + n_phi); 
    gsl_vector_view B_view = gsl_vector_view_array(B, n * nom + n_phi); 
    int signum; 
    gsl_linalg_LU_decomp(&A_view.matrix, p, &signum); 
    gsl_linalg_LU_solve(&A_view.matrix, p, &B_view.vector, x);
}



void rk4(int steps,double dt, gsl_vector* ddr_temp1, gsl_vector* ddr_temp2, gsl_vector* ddr_temp3, gsl_vector* ddr_temp4, double* temp, double* results, double* M, double* J, double* A, double* B, gsl_permutation* p,double* energy) {
    static double last_printed = -0.1;
    int energy_index = 0;
    for (int step = 0; step < steps; ++step) {
        double kinetic = 0;
        double potential = 0;
        double percent = (double)step * 100.0 / steps;
        if (percent - last_printed >= 0.1) {
            printf("%.1f%%\n", percent);
            last_printed = percent;
        }
        solve(ddr_temp1, M, J, A, B, p);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            kinetic += 0.5 * (pow(masses_array[j].dx, 2) + pow(masses_array[j].dy, 2));
            potential += g * masses_array[j].y;
            temp[index] = masses_array[j].dx;
            temp[index + 1] = masses_array[j].dy;
            temp[index + 2] = masses_array[j].x;
            temp[index + 3] = masses_array[j].y;
            results[step * nom * 2 + 2 * j] = masses_array[j].x;
            results[step * nom * 2 + 2 * j + 1] = masses_array[j].y;
            masses_array[j].ddx = gsl_vector_get(ddr_temp1, 2 * j);
            masses_array[j].ddy = gsl_vector_get(ddr_temp1, 2 * j + 1);
            masses_array[j].dx = masses_array[j].dx + masses_array[j].ddx * (dt / 2);
            masses_array[j].dy = masses_array[j].dy + masses_array[j].ddy * (dt / 2);
            masses_array[j].x = masses_array[j].x + masses_array[j].dx * (dt / 2);
            masses_array[j].y = masses_array[j].y + masses_array[j].dy * (dt / 2);
        }


        solve(ddr_temp2, M, J, A, B, p);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            masses_array[j].ddx = (gsl_vector_get(ddr_temp2, 2 * j));
            masses_array[j].ddy = (gsl_vector_get(ddr_temp2, 2 * j + 1));
            masses_array[j].dx = temp[index] + masses_array[j].ddx * (dt / 2);
            masses_array[j].dy = temp[index + 1] + masses_array[j].ddy * (dt / 2);
            masses_array[j].x = temp[index + 2] + masses_array[j].dx * (dt / 2);
            masses_array[j].y = temp[index + 3] + masses_array[j].dy * (dt / 2);
        }


        solve(ddr_temp3, M, J, A, B, p);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            masses_array[j].ddx = (gsl_vector_get(ddr_temp3, 2 * j));
            masses_array[j].ddy = (gsl_vector_get(ddr_temp3, 2 * j + 1));
            masses_array[j].dx = temp[index] + masses_array[j].ddx * dt;
            masses_array[j].dy = temp[index + 1] + masses_array[j].ddy * dt;
            masses_array[j].x = temp[index + 2] + masses_array[j].dx * dt;
            masses_array[j].y = temp[index + 3] + masses_array[j].dy * dt;
        }


        solve(ddr_temp4, M, J, A, B, p);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            masses_array[j].ddx = (gsl_vector_get(ddr_temp1, 2 * j) + 2 * gsl_vector_get(ddr_temp2, 2 * j) + 2 * masses_array[j].ddx + gsl_vector_get(ddr_temp4, 2 * j)) / 6;
            masses_array[j].ddy = (gsl_vector_get(ddr_temp1, 2 * j + 1) + 2 * gsl_vector_get(ddr_temp2, 2 * j + 1) + 2 * masses_array[j].ddy + gsl_vector_get(ddr_temp3, 2 * j + 1)) / 6;
            masses_array[j].dx = temp[index] + masses_array[j].ddx * dt;
            masses_array[j].dy = temp[index + 1] + masses_array[j].ddy * dt;
            masses_array[j].x = temp[index + 2] + masses_array[j].dx * dt;
            masses_array[j].y = temp[index + 3] + masses_array[j].dy * dt;
        }

        if (step % (steps/250) == 0) {
            energy[energy_index] = kinetic + potential;
            energy_index++;
        }
    }
}