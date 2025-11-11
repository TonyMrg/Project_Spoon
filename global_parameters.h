#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL.h>
#include <stdarg.h>
#include <string.h>
#include <mkl.h>
#include <time.h>
#include <omp.h>
#include <GL/gl.h>
#include <GL/glu.h>



#ifndef GLOBALS_H
#define GLOBALS_H

extern double g;
extern int nom;
extern double m;                            //Parameters
extern int n_phi;
extern int n;
extern double l;
extern struct mass* masses_array;

struct mass {
    double x;
    double y;
    double z;
    double dx;
    double dy;
    double dz;
    double ddx;
    double ddy;
    double ddz;
    double m;
};

#endif