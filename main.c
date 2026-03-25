#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL.h>
#include <stdarg.h>

//// Simulation parameters
//const double g = 9.81;
//const int nom = 20;
//const double m = 1.0;
//const int n_phi = 20;  //number of constraints
//const int n = 2;
//const double l = 1.0;
//
//double* L;
//
//struct mass* masses_array;
//
//double* buildB();
//double* buildJ();
//double* combine_matrixA(double* M, double* J);
//double* buildM();
//void row_to_col_major(double* src, double* dst, int rows, int cols);

//struct mass {
//    double x;
//    double y;
//    double dx;
//    double dy;
//    double ddx;
//    double ddy;
//    double m;
//};


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
        double kinetic = 0;
        double potential = 0;
        double* temp = malloc(4 * nom * sizeof(double));
        if (!temp) {
            fprintf(stderr, "Error: malloc failed for temp array\n");
            return 1;
        }
        
          
        gsl_vector* ddr_temp1 = gsl_vector_alloc(n * nom + n_phi);
        solve(ddr_temp1);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            kinetic += 0.5 * (pow(masses_array[j].dx, 2)+pow(masses_array[j].dy,2));
            potential += g * masses_array[j].y;
            temp[index] = masses_array[j].dx;
            temp[index + 1] = masses_array[j].dy;
            temp[index + 2] = masses_array[j].x;
            temp[index + 3] = masses_array[j].y;
            results[step * nom * 2 + 2*j] = masses_array[j].x;
            results[step * nom * 2 + 2*j + 1] = masses_array[j].y;
            masses_array[j].ddx = gsl_vector_get(ddr_temp1, 2*j);
            masses_array[j].ddy = gsl_vector_get(ddr_temp1, 2*j + 1);
            masses_array[j].dx = masses_array[j].dx + masses_array[j].ddx * 0.0005;
            masses_array[j].dy = masses_array[j].dy + masses_array[j].ddy * 0.0005;
            masses_array[j].x = masses_array[j].x + masses_array[j].dx * 0.0005;
            masses_array[j].y = masses_array[j].y + masses_array[j].dy * 0.0005;
        }

        gsl_vector* ddr_temp2 = gsl_vector_alloc(n * nom + n_phi);
        solve(ddr_temp2);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            masses_array[j].ddx = (gsl_vector_get(ddr_temp2, 2 * j));
            masses_array[j].ddy = (gsl_vector_get(ddr_temp2, 2 * j + 1));
            masses_array[j].dx = temp[index] + masses_array[j].ddx * 0.0005;
            masses_array[j].dy = temp[index+1] + masses_array[j].ddy * 0.0005;
            masses_array[j].x = temp[index+2] + masses_array[j].dx * 0.0005;
            masses_array[j].y = temp[index+3] + masses_array[j].dy * 0.0005;
        }

        gsl_vector* ddr_temp3 = gsl_vector_alloc(n * nom + n_phi);
        solve(ddr_temp3);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            masses_array[j].ddx = (gsl_vector_get(ddr_temp3, 2 * j));
            masses_array[j].ddy = (gsl_vector_get(ddr_temp3, 2 * j + 1));
            masses_array[j].dx = temp[index] + masses_array[j].ddx * 0.001;
            masses_array[j].dy = temp[index + 1] + masses_array[j].ddy * 0.001;
            masses_array[j].x = temp[index + 2] + masses_array[j].dx * 0.001;
            masses_array[j].y = temp[index + 3] + masses_array[j].dy * 0.001;
        }

        gsl_vector* ddr_temp4 = gsl_vector_alloc(n * nom + n_phi);
        solve(ddr_temp4);
        for (int j = 0;j < nom;j++) {
            int index = j * 4;
            masses_array[j].ddx = (gsl_vector_get(ddr_temp1, 2 * j) + 2* gsl_vector_get(ddr_temp2, 2 * j)+ 2 * masses_array[j].ddx + gsl_vector_get(ddr_temp4, 2 * j) ) / 6;
            masses_array[j].ddy = (gsl_vector_get(ddr_temp1, 2 * j + 1) + 2 * gsl_vector_get(ddr_temp2, 2 * j + 1) + 2 * masses_array[j].ddy + gsl_vector_get(ddr_temp3, 2 * j + 1)) / 6;
            masses_array[j].dx = temp[index] + masses_array[j].ddx * 0.001;
            masses_array[j].dy = temp[index + 1] + masses_array[j].ddy * 0.001;
            masses_array[j].x = temp[index + 2] + masses_array[j].dx * 0.001;
            masses_array[j].y = temp[index + 3] + masses_array[j].dy * 0.001;
        }
        
        
        
        free(temp);
        gsl_vector_free(ddr_temp1);
        gsl_vector_free(ddr_temp2);
        gsl_vector_free(ddr_temp3);
        gsl_vector_free(ddr_temp4);


        if (step % 40 == 0) {
            energy[energy_index] = kinetic + potential;
            energy_index++;
        }
    }

    save_energy(energy, 250);
    visualize_double_pendulum(results, steps, 1,nom);
    return 0;
}
