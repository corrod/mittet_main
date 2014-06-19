#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fwi3d_cpu_alloex1_v1.1.h"

void lattice_time_3d()
{
    double kk,tmax,t0;
// total number of lattice
    nx = ix - 1;
    ny = iy - 1;
    nz = iz - 1;

// time step and courant condition
    //courant = 1.0/CC/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    //dt = courant*dt_ratio;
    ////courant_ori = 1.0/CC/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    ////vc = sqrt( 2.f*omega_0 / MU0 /0.7f);
    cmin = sqrt(2.f*omega_0/MU0/1.f);
    cmax = sqrt(2.f*omega_0/MU0/1.f);
    //kk = 1.f/sqrt(3.f);
    //courant = kk*dx/cmax;
    courant = 1.0f/cmax/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    dt = courant*0.8f;
    //tmax = 1.0*(double)ix*dx/cmin;
    //t0   = M_PI/fmax_w;
    //it   = (int)((tmax+2.f*t0)/dt);

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
    printf("cmax       ; %f\n",cmax);
    printf("cmin       ; %f\n",cmin);
}

