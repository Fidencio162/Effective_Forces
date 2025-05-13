#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define mp 1000
#define mr 2048
// number of particles
#define np 8200

//Random number generator parameters
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
// parameters of the potential (WCA: lambda, T*)
#define dl  50.0
#define dT  0.55954475
//////////////////////////////////////

double boxl, rc;
double del;
double beta, alpha;
int np1;
long iseed;

float ranf(long *idum);
double iniconf(double xc[], double yc[],double zc[], double d, double rc);
void mcmove(double x[], double y[], double z[], double sigma[], double zz[], double *ener, int *nattemp, int *nacc, double del, double boxl);
void adjust(int nattemp, int nacc, double *del);
void average(double x[], double y[], double z[], double g11[], double g12[], double g22[], int *nattemp, int *nacc, double del, double dr);

int main() {
    double xc[np], yc[np], zc[np];
    double diam[np], zz[np];
    double g11[mr], g12[mr], g22[mr];
    double diam1, diam2, xi1, xi2, rho1, rho2, dv1, dv2, dv12;
    double aux1, aux2, aux3;
    int nz1, nz2, np2, nn;
    iseed = 123456789;
    double rho = 0.00078304232; // rediced density
    diam1 = 10.0; // diameter of the cation
    diam2 = 1.0;  // diameter of the anion
    xi1=1.0/41.0;  // cation ionic concentration
    xi2=1.0-xi1;  // anion ionic concentration
    boxl = pow(np / rho, 1.0 / 3.0); // box length
    rc = boxl / 2.0; // cut-off radius
    beta = 1.0/dT;
    nz1 = 40; // cation charge
    nz2 = -1;  // anion charge
    np1=-nz2*np/(nz1-nz2); // positive particles
    np2=np-np1;  // negative particles
    rho1 = rho*xi1; // reduced cation density
    rho2 = rho*xi2;  // reduced anion density
    double phi1=M_PI*rho1*pow(diam1,3)/6.0;
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
    //printf("seed: %f\n", ranf(&iseed));    
   //abort();
        ////////////////////////////////////////////////
    // Allocated of the charge and size
    for (int l = 0; l < np1; l++) {
       zz[l] = (double)nz1;
       diam[l] = diam1;
    }
    for (int l = np1; l < np; l++) {
       zz[l] = (double)nz2;
       diam[l] = diam2;
    }
    //printf("zn1 = %f\n", zz[np1]);
    //iniconf( xc, yc, zc, d, rc);
        FILE *initconf;    
    initconf = fopen("init_config_no_overlap", "r"); // 'r' means read
    if (initconf == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }
    for (int i = 0; i < np; i++) {
        fscanf(initconf, "%d %lf %lf %lf %lf %lf %lf", &nn, &xc[i], &yc[i], &zc[i], &aux1, &aux2, &aux3);
    }
    
    FILE *fp;
    fp = fopen("iniconf.dat", "w"); // Open file for writing ("w" mode)

    if (fp == NULL) {
        perror("Error opening file");
        return -1;
    }

    for (int i = 0; i < np; ++i) {
        fprintf(fp, "%.8f %.8f %.8f\n", xc[i], yc[i], zc[i]); // Write x(i), y(i), z(i) to file
    }

    fclose(fp); // Close the file
    
    FILE *file;    
    file = fopen("energy-Mie.dat", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }
    
    FILE *fc;
    fc = fopen("confequi.dat", "w"); // Open file for writing ("w" mode)

    if (fc == NULL) {
        perror("Error opening file");
        return -1;
    }
    
    del = 0.1;
int nattemp = 0;
int nacc = 1;
double ener = 0.0;

///// Thermalization process////////////////////////////////////////
   // Main loop
    for (int i = 1; i <= 30000000; i++) {
        mcmove(xc, yc, zc, diam, zz, &ener, &nattemp, &nacc, del, boxl);
        adjust(nattemp, nacc, &del);

        if (i % 100 == 0) {
            fprintf(file, "%.1f %.8f\n", (double)i, ener / np);
        }
//        if (i % 1000 == 0) {
//            printf("%d %.16f %.16f\n", i, del, ener / np);
//        }
    }
    
    for (int i = 0; i < np; ++i) {
        fprintf(fc, "%.8f %.8f %.8f\n", xc[i], yc[i], zc[i]); // Configuration in equilibruim
    }
    fclose(fc);
    fclose(file);
    printf("The system has thermalized\n");
    abort();
////////////////////////////////////////////////////////////////////    
    // Calculation of thermodynamics properties
    
    FILE *f_gr;    
    f_gr = fopen("gr1.dat", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return -2;
    }
    
    double r[mr];
    // Loop to initialize elements of g(r) to 0.0
    for (int i = 0; i <= mr; i++) {
        g11[i] = 0.0;
        g12[i] = 0.0;
        g22[i] = 0.0;
    }
    int nav=0; //initialize number to average
        for (int i = 1; i <= 30000000; i++) {
        mcmove(xc, yc, zc, diam, zz, &ener, &nattemp, &nacc, del, boxl);
        adjust(nattemp, nacc, &del);

        if (i % np == 0) {
            //Calculating g(r)
            average(xc, yc, zc, g11, g12, g22, &nattemp, &nacc, del, dr);
            nav++;
        }
    }
    
    printf("Average number = %i\n", nav);
    for (int i = 2; i <= mr; i++) {
        r[i] = (i-1) * dr;
        dv1 = (4.0 * M_PI * pow(r[i], 2) * dr) * rho1;
        dv12 = (4.0 * M_PI * pow(r[i], 2) * dr) * rho2;
        dv2 = (4.0 * M_PI * pow(r[i], 2) * dr) * rho2;
        g11[i] = g11[i] / (np1 * nav * dv1);
        g12[i] = g12[i] / (np1 * nav * dv12);
        g22[i] = g22[i] / (np2 * nav * dv2);
        fprintf(f_gr,"%.8f %.8f %.8f %.8f\n", r[i], g11[i], g12[i], g22[i]);
    }
    fclose(f_gr);

    return 0;
}

double iniconf(double xc[], double yc[],double zc[], double d, double rc)
{
  xc[0]=d/2.0-rc;
  yc[0]=d/2.0-rc;
  zc[0]=d/2.0-rc;
  for (int i=1; i<=50; i++) {
    xc[i]=xc[i-1]+12;
    yc[i]=yc[i-1];
    zc[i]=zc[i-1];
    if (xc[i] > rc)
    {
      xc[i]=xc[0];
      yc[i]=yc[i-1]+12;
      zc[i]=zc[i-1];
      if (yc[i] > rc)
      {
        xc[i]=xc[0];
        yc[i]=yc[0];
        zc[i]=zc[i-1]+12;
      }
    }
  }
  
  for (int i=51; i<np; i++) {
    xc[i]=xc[i-1]+0.913*d;
    yc[i]=yc[i-1];
    zc[i]=zc[i-1];
    if (xc[i] > rc)
    {
      xc[i]=xc[0];
      yc[i]=yc[i-1]+0.913*d;
      zc[i]=zc[i-1];
      if (yc[i] > rc)
      {
        xc[i]=xc[0];
        yc[i]=yc[0];
        zc[i]=zc[i-1]+0.913*d;
      }
    }
  }
}
///////////////
//This configuration calculates the pair potential between particles i & j
void potential(double rij, double diam[], double zz[], double *uij, int i, int j) {
    double erfcc, erfcij, c2api, Bcut, Aij, sigmaij, Bij;

    erfcij = erfc(alpha*rij);
    erfcc = erfc(alpha*rc)/rc;
    c2api = 2.0*alpha/sqrt(M_PI);
    Bcut = erfcc/rc + c2api*exp(-pow(alpha*rc,2))/rc;
    Aij = zz[i]*zz[j];
    sigmaij = 0.5*(diam[i]+diam[j]);
    /////////////////////////////////////////////////////////////
    Bij = fabs(Aij)*(pow(sigmaij,8)*erfc(alpha*sigmaij)+c2api*pow(sigmaij,9)*exp(-pow(alpha*sigmaij,2))-pow(sigmaij,10)*Bcut)/9.0;
    
    *uij = Bij/pow(rij,9)+Aij*(erfcij/rij-erfcc);
}

void energy(double x[],double y[],double z[], double xj, double yj, double zj, double diam[], double zz[], double *ener, int j) {
    double ener_local = 0.0;
    double xij, yij, zij, rij, rij2, uij;
	
    for (int i = 0; i < j; ++i) {
        xij = x[i] - xj;
        yij = y[i] - yj;
        zij = z[i] - zj;

        // Apply periodic boundary conditions
        xij -= boxl * round(xij / boxl);
        yij -= boxl * round(yij / boxl);
        zij -= boxl * round(zij / boxl);

        rij2 = xij * xij + yij * yij + zij * zij;
        rij = sqrt(rij2);

        if (rij < rc) {
                potential(rij, diam, zz, &uij, i, j); // Call potential subroutine
            ener_local += uij;
        }
    }

    for (int i = j+1; i < np; ++i) {
        xij = x[i] - xj;
        yij = y[i] - yj;
        zij = z[i] - zj;
        
        // Apply periodic boundary conditions
        xij -= boxl * round(xij / boxl);
        yij -= boxl * round(yij / boxl);
        zij -= boxl * round(zij / boxl);

        rij2 = xij * xij + yij * yij + zij * zij;
        rij = sqrt(rij2);

        if (rij < rc) {
                potential(rij, diam, zz, &uij, i, j);
            ener_local += uij;
        }
    }

    *ener = ener_local;
    //printf("ener = %Lf\n", *ener);
}
//This subroutine displace the system to a new configuration
void mcmove(double x[], double y[], double z[], double diam[], double zz[], double *ener, int *nattemp, int *nacc, double del, double boxl) {
   //double xx[np], yy[np], zz[np];
    *nattemp = *nattemp + 1;
    int no = (int)(ranf(&iseed) * np); // rand() returns a pseudo-random integer between 0 and RAND_MAX
    double xnij, ynij, znij, rnij2, distance, xmij, ymij, zmij, rmij2, r12_old, r12_new;
    double enero_no, enero_nj,enern_no, enern_nj, enero12, enern12, enero, enern, x12, y12, z12;
    double xo = x[no];
    double yo = y[no];
    double zo = z[no];
    
    ///// Chantal algorithm
    int neighbors[np];
    int neighbor_count = 0;
    for (int j = 0; j < np; j++) {
        if (j == no) continue;
        xnij = x[j] - xo;
        ynij = y[j] - yo;
        znij = z[j] - zo;

        rnij2 = xnij * xnij + ynij * ynij + znij * znij;
        distance = sqrt(rnij2);
        if (distance <= 0.08*boxl) {
            neighbors[neighbor_count++] = j;
        }
    }

    if (neighbor_count == 0) return;
    //printf("%d \n", neighbor_count); //

    // Pick random neighbor
    int nj = neighbors[rand() % neighbor_count];
    //printf("%d %d \n", no, nj);
        xmij = x[nj] - xo;
        ymij = y[nj] - yo;
        zmij = z[nj] - zo;
        rmij2 = xmij * xmij + ymij * ymij + zmij * zmij;
        r12_old = sqrt(rmij2);
/////

    energy(x,y,z,xo,yo,zo, diam, zz, &enero, no);
    enero_no=enero;
    energy(x,y,z,x[nj],y[nj],z[nj], diam, zz, &enero, nj);  // Arreglar cálculo de las energías enero1+enero2
    enero_nj=enero;
    enero12=enero_no+enero_nj;

    x[no] = x[no] + (ranf(&iseed) - 0.5) * del;
    y[no] = y[no] + (ranf(&iseed) - 0.5) * del;
    z[no] = z[no] + (ranf(&iseed) - 0.5) * del;
    ///////////////////////////////////////////
    x[nj] = x[nj] + (ranf(&iseed) - 0.5) * del;
    y[nj] = y[nj] + (ranf(&iseed) - 0.5) * del;
    z[nj] = z[nj] + (ranf(&iseed) - 0.5) * del;

    // Apply periodic boundary conditions
    x[no] -= boxl * round(x[no] / boxl);
    y[no] -= boxl * round(y[no] / boxl);
    z[no] -= boxl * round(z[no] / boxl);
    ////////////////////////////////////
    x[nj] -= boxl * round(x[nj] / boxl);
    y[nj] -= boxl * round(y[nj] / boxl);
    z[nj] -= boxl * round(z[nj] / boxl);
    
    // Aquí
    // Calcular r12_new
    x12 = x[no]-x[nj];
    y12 = y[no]-y[nj];
    z12 = z[no]-z[nj];
    
    r12_new = x12*x12 + y12*y12 + z12*z12;

    energy(x,y,z,x[no],y[no],z[no],diam,zz, &enern, no);
      enern_no=enern;
    energy(x,y,z,x[nj],y[nj],z[nj],diam,zz, &enern, nj);
      enern_nj=enern;
      enern12=enern_no+enern_nj;
    if (ranf(&iseed) < pow(r12_new / r12_old, 2)*exp(-beta*(enern12 - enero12))) {
        *ener = *ener + 0.5*(enern12 - enero12);
        *nacc = *nacc + 1;
    } else {
        x[no] = xo;
        y[no] = yo;
        z[no] = zo;
        
        x[nj] = xmij;
        y[nj] = ymij;
        z[nj] = zmij;
    }
}
// This subroutine calculates the g(r)
void average(double x[], double y[], double z[], double g11[], double g12[], double g22[], int *nattemp, int *nacc, double del, double dr) {
    (*nattemp)++;

    // Variables for loop indices and distance calculations
    int i, j, nbin;
    double xij, yij, zij, rij2, rij;

    // Calculating the g(r)
    //g11
    for (i = 0; i < np1 - 1; i++) {
        for (j = i + 1; j < np1; j++) {
            xij = x[j] - x[i];
            yij = y[j] - y[i];
            zij = z[j] - z[i];

            xij -= boxl * nearbyint(xij / boxl);
            yij -= boxl * nearbyint(yij / boxl);
            zij -= boxl * nearbyint(zij / boxl);

            rij2 = xij * xij + yij * yij + zij * zij;
            rij = sqrt(rij2);
            if (rij < rc){
              nbin = (int)(rij / dr);
              if (nbin < mr) {
                g11[nbin] += 2.0;
              }
            }
        }
    }
    
    //g12
    for (i = 0; i < np1 ; i++) {
        for (j = np1; j < np; j++) {
            xij = x[j] - x[i];
            yij = y[j] - y[i];
            zij = z[j] - z[i];

            xij -= boxl * nearbyint(xij / boxl);
            yij -= boxl * nearbyint(yij / boxl);
            zij -= boxl * nearbyint(zij / boxl);

            rij2 = xij * xij + yij * yij + zij * zij;
            rij = sqrt(rij2);
            if (rij < rc){
              nbin = (int)(rij / dr);
              if (nbin < mr) {
                g12[nbin] += 1.0;
              }
            }
        }
    }
    
        //g22
    for (i = np1; i < np-1 ; i++) {
        for (j = i+1; j < np; j++) {
            xij = x[j] - x[i];
            yij = y[j] - y[i];
            zij = z[j] - z[i];

            xij -= boxl * nearbyint(xij / boxl);
            yij -= boxl * nearbyint(yij / boxl);
            zij -= boxl * nearbyint(zij / boxl);

            rij2 = xij * xij + yij * yij + zij * zij;
            rij = sqrt(rij2);
            if (rij < rc){
              nbin = (int)(rij / dr);
              if (nbin < mr) {
                g22[nbin] += 2.0;
              }
            }
        }
    }
    (*nacc)++;
}
//Random generator algorithm
float ranf(long *idum) {
    static long iv[NTAB];
    static long iy = 0;
    static long idum2 = 123456789;
    float temp;
    int j;
    long k;

    if (*idum <= 0) {
       if (-(*idum) < 1) *idum=1;
       else *idum = -(*idum);
       idum2=(*idum);
       for (j=NTAB+7;j>=0;j--) {
           k=(*idum)/IQ1;
           *idum=IA1*(*idum-k*IQ1)-k*IR1;
           if (*idum < 0) *idum += IM1;
           if (j < NTAB) iv[j] = *idum;
       }
       iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

//This subroutine adjusts the displacement of particles
void adjust(int nattemp, int nacc, double *dr) {
    if (nattemp % nacc == 0) {
        double ratio = (double)nacc / (double)nattemp;
        if (ratio > 0.5) {
            *dr *= 1.05;
        } else {
            *dr *= 0.95;
        }
    }
}
