#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL.h>
#include <stdarg.h>
#include <string.h>
#include <mkl.h>
#include <time.h>
#include <omp.h>


#ifndef GLOBALS_H
#define GLOBALS_H

// Simulation parameters
extern double g;
extern int nom;
extern double m;
extern int n_phi;
extern int n;
extern double l;
extern struct mass* masses_array;

struct mass {
    double x;
    double y;
    double dx;
    double dy;
    double ddx;
    double ddy;
    double m;
};

#endif