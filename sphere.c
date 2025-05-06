#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define mp 1000
#define mr 2048
// number of particles
#define np 200
#define R 7.0

double boxl, rc;
double del;
double beta, alpha;
int np1;
long iseed;

float ranf(long *idum);
double iniconf(double xc[], double yc[],double zc[], double d, double rc);

void config(double xs[], double ys[],double zs[], double xi, double yi,double zi);

int main() {
double xc[np], yc[np], zc[np];
double xp[np], yp[np], zp[np];
np1=100;
double xs[np1], ys[np1], zs[np1];
double xi, yi, zi;
int nn;
     FILE *initconf;    
    initconf = fopen("init_config.dat", "r"); // 'r' means read
    if (initconf == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }
    for (int i = 0; i < np; i++) {
        fscanf(initconf, "%lf %lf %lf", &xp[i], &yp[i], &zp[i]);
    }

FILE *init;
    init = fopen("init_config", "w"); // Open file for writing ("w" mode)

    if (init == NULL) {
        perror("Error opening file");
        return -1;
    }
    
    for (int i=0; i<np; i++) {
        fprintf(init,"%.8f %.8f %.8f\n", xp[i], yp[i], zp[i]); // Write x(i), y(i), z(i) to file
   }
    
    for (int j=0; j<np; j++) {
    xi=xp[j];
     yi=yp[j];
      zi=zp[j];
config(xs, ys, zs, xi, yi, zi);
for (int i=0; i<np1; i++) {
        fprintf(init,"%.8f %.8f %.8f\n", xs[i], ys[i], zs[i]); // Write x(i), y(i), z(i) to file
   }
     }
   fclose(init);
return 0;
}

void config(double xs[], double ys[],double zs[], double xi, double yi,double zi) {
double phi[np1], theta[np1];
double pi, gold_r;
//double xi, yi, zi;
       pi=4.0*atan(1.0);
      gold_r=(1.0+sqrt(5.0))/2.0;
      for (int i=0; i<np1; i++) {
        phi[i]=2.0*pi*i/gold_r;
        theta[i]=acos(1.0-2.0*i/(np1-1.0));
        xs[i]=R*sin(theta[i])*cos(phi[i])+xi;
        ys[i]=R*sin(theta[i])*sin(phi[i])+yi;
        zs[i]=R*cos(theta[i])+zi;
      }
}
