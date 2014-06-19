#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fwi3d_cpu_alloex1_v1.4.h"

void lattice_time_3d()
{
    int i,j,k,ijk;
    double kk,tmax,t0;
    double sigmax,sigmin;
// total number of lattice
    nx = ix - 1;
    ny = iy - 1;
    nz = iz - 1;

// time step and courant condition
    //courant = 1.0/CC/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    //dt = courant*dt_ratio;
    ////courant_ori = 1.0/CC/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    ////vc = sqrt( 2.f*omega_0 / MU0 /0.7f);
    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk = k*ix*iy + j*ix + i;
          vel[ijk] = sqrt(2.f*omega_0 /MU0 /sig[ijk]);
        }
      }
    }
    cmax = vel[0];
    cmin = vel[0];
    for(i=0;i<ix*iy*iz;i++){
      if(vel[i] > cmax) cmax = vel[i];
      if(vel[i] < cmin) cmin = vel[i];
    }
    cmax = sqrt(2.f*omega_0 /MU0 /minval);
    cmin = sqrt(2.f*omega_0 /MU0 /maxval);

    courant = 1.0f/cmax/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    dt = courant*6.f/7.f*0.999f;

    printf("inputted CELL length: dx: %lf, dy: %lf, dz: %lf \n",dx,dy,dz);
    printf("delta t : %e\n", dt);
    printf("courant : %e\n", courant);
    printf("check dt: %e < %e \n",dt,courant_ori);
    if(dt>courant){
        printf(" **** violated courant stability condition ***\n");
    }
    printf("\n");
    printf("total time : %e\n",it*dt);
    printf("total step : %d\n",it);
    printf("lattice size ix: %d,  iy: %d,  iz: %d\n",ix,iy,iz);
    printf("cmax       : %e\n",cmax);
    printf("cmin       : %e\n",cmin);
}

