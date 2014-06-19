#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fftw3.h"
#define nmedia 3
int ix,iy,iz;
double dx,dy,dz;
void read_shotrec();

int main(){
    int i,j,k,ijk;
    int iseabed, NN;
    double type[100], *sig;
    double sigmax,sigmin;
    double center;
    int xstart, xend;
    int ystart, yend;
    int zstart, zend;
    FILE *ofm1,*ofm2,*ofm3,*ofm4;

    read_shotrec();

    NN=ix*iy*iz;
    sig = (double *)malloc(sizeof(double) * NN);
    ofm1 = fopen("./data/conductivity.dat","w");
    ofm2 = fopen("./data/minmax.dat","w");
    ofm3 = fopen("./data/waku.dat","w");
    ofm4 = fopen("./data/iseabed.dat","w");

// id = 0 test
    type[0] = 0.45f;
// id = 1 seawater
    type[1] = 3.3f;
// id = 2 basalt
    type[2] = 5.0f;
// id = 3 submarine massive sulphide
    type[3] = 10.0f;

// seawater and basalt's boundary
    iseabed = 15;
    sigmax = type[0];
    sigmin = type[0];
    sigmax = 7.f;
    sigmin = 0.4f;
    for(i=0;i<nmedia;i++){
      if(type[i] > sigmax)  sigmax = type[i];
      if(type[i] < sigmin)  sigmin = type[i];
    }

    center =16;
    xstart = center -2;
    xend   = center +2;
    ystart = center -2;
    yend   = center +2;
    zstart = iseabed-3;
    zend   = iseabed;
    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk = k*iy*ix + j*ix + i;
          sig[ijk]=type[0];
          if(k>iseabed) sig[ijk] = type[1];
          if(i >= xstart-1 && i < xend \
          && j >= ystart-1 && j < yend \
          && k >= zstart && k <= zend) sig[ijk] = type[2];
          fprintf(ofm1,"%12lf   %12lf   %12lf   %15.5e\n", i*dx, j*dy,k*dz, sig[ijk]);
        }
      }
    }
    fprintf(ofm2,"%lf\n",sigmin);
    fprintf(ofm2,"%lf\n",sigmax);

    fprintf(ofm3,"%12lf    %12lf\n", xstart*dx, zstart*dz);
    fprintf(ofm3,"%12lf    %12lf\n", xend  *dx, zstart*dz);
    fprintf(ofm3,"%12lf    %12lf\n", xend  *dx, zend  *dz);
    fprintf(ofm3,"%12lf    %12lf\n", xstart*dx, zend  *dz);
    fprintf(ofm3,"%12lf    %12lf\n", xstart*dx, zstart*dz);
    fprintf(ofm4,"%12d\n",iseabed);
    fclose(ofm1);
    fclose(ofm2);
    fclose(ofm3);
    fclose(ofm4);
}
///////////////////////////////////////////////
void read_shotrec()
{
  int i,itmp;
  double dtmp;
  FILE *ifp1;
  ifp1 = fopen("./data/model_env.dat","r");

  fscanf(ifp1,"%d",&ix);
  fscanf(ifp1,"%d",&iy);
  fscanf(ifp1,"%d",&iz);
  fscanf(ifp1,"%lf",&dx);
  fscanf(ifp1,"%lf",&dy);
  fscanf(ifp1,"%lf",&dz);
  fscanf(ifp1,"%lf",&dtmp);
  fscanf(ifp1,"%d",&itmp);
  fscanf(ifp1,"%d",&itmp);
  fscanf(ifp1,"%d",&itmp);
  fscanf(ifp1,"%lf",&dtmp);
  fscanf(ifp1,"%lf",&dtmp);
  fscanf(ifp1,"%lf",&dtmp);
  fscanf(ifp1,"%d",&itmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%d",&itmp);

}
///////////////////////////////////////////////
