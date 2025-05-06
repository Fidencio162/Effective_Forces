#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define mp 1000
#define mr 2048
// number of particles
#define np 200
#define R 6.5

double boxl, rc;
double del;
double beta, alpha;
int np1;
long iseed;

float ranf(long *idum);
void iniconf(double xc[], double yc[],double zc[], double d, double rc);

void config(double xs[], double ys[],double zs[], double xi, double yi,double zi);

int main() {
double xc[np], yc[np], zc[np];
    double diam1, diam2, xi1, xi2, rho1, rho2, dv1, dv2, dv12;
    //double aux1, aux2, aux3;
    int nz1, nz2, np2;
    iseed = 123456789;
    double rho = 0.000025; // rediced density
    diam1 = 10.0; // diameter of the cation
    diam2 = 1.0;  // diameter of the anion
    xi1=1.0/101.0;  // cation ionic concentration
    xi2=1.0-xi1;  // anion ionic concentration
    boxl = pow(np / rho, 1.0 / 3.0); // box length
    rc = boxl / 2.0; // cut-off radius
    //beta = 1.0/dT;
    nz1 = 100; // cation charge
    nz2 = -1;  // anion charge
    np1=-nz2*np/(nz1-nz2); // positive particles
    np2=np-np1;  // negative particles
    rho1 = rho*xi1; // reduced cation density
    rho2 = rho*xi2;  // reduced anion density
    double phi1=M_PI*rho*pow(diam1,3)/6.0;
    double phi2=M_PI*rho2*pow(diam2,3)/6.0;
    double phi = phi1+phi2; // packing fraction
    alpha = 4.0/boxl;  // Wolf splitting parameter
    double dr = rc / mr;
    printf("rho_M = %.8f\n", rho1);
    printf("phi_M = %.8f\n", phi1);
    printf("phi = %f\n", phi);
    printf("dr = %f\n", dr);
    double d=pow(rho, -1.0 / 3.0);
    printf("Mean interparticle distance: %f\n", d);
    printf("n1 = %d\n", np1);
    printf("n2 = %d\n", np2);
    printf("boxl = %f\n", boxl);

FILE *init;
    init = fopen("init_config.dat", "w"); // Open file for writing ("w" mode)

    if (init == NULL) {
        perror("Error opening file");
        return -1;
    }
    iniconf( xc, yc, zc, d, rc);
    for (int i = 0; i < np; ++i) {
        fprintf(init, "%.8f %.8f %.8f\n", xc[i], yc[i], zc[i]); // Write x(i), y(i), z(i) to file
    }
   fclose(init);
return 0;
}

void iniconf(double xc[], double yc[],double zc[], double d, double rc)
{
  xc[0]=d/2.0-rc;
  yc[0]=d/2.0-rc;
  zc[0]=d/2.0-rc;
  for (int i=1; i<np; i++) {
    xc[i]=xc[i-1]+d;
    yc[i]=yc[i-1];
    zc[i]=zc[i-1];
    if (xc[i] > rc)
    {
      xc[i]=xc[0];
      yc[i]=yc[i-1]+d;
      zc[i]=zc[i-1];
      if (yc[i] > rc)
      {
        xc[i]=xc[0];
        yc[i]=yc[0];
        zc[i]=zc[i-1]+d;
      }
    }
  }
}
