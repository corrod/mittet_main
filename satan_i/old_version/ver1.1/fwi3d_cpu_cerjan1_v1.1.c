#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fwi3d_cpu_alloex1_v1.1.h"
#include "fwi3d_cpu_cerjan1_v1.1.h"
//////////////////////////////////////////
void cerjan_e(double *EX,double *EY,double *EZ)
{
    int i,j,k,ijk;
    int beta = 20;
    double alpha = 5.5e-2;

#ifdef _OPENMP
//#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
//    #pragma omp for
    #endif
    for(i=0;i<beta;i++){
        for(j=0;j<iy;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                EX[ijk] *= exp(-alpha*pow(beta-i,2.f));
                EY[ijk] *= exp(-alpha*pow(beta-i,2.f));
                EZ[ijk] *= exp(-alpha*pow(beta-i,2.f));
            }
        }
    }
    #ifdef _OPENMP
//    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=0;j<beta;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                EX[ijk] *= exp(-alpha*pow(beta-j,2.f));
                EY[ijk] *= exp(-alpha*pow(beta-j,2.f));
                EZ[ijk] *= exp(-alpha*pow(beta-j,2.f));
            }
        }
    }
    #ifdef _OPENMP
//    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=0;j<iy;j++){
            for(k=0;k<beta;k++){
                ijk = k*ix*iy + j*ix + i;
                EX[ijk] *= exp(-alpha*pow(beta-k,2.f));
                EY[ijk] *= exp(-alpha*pow(beta-k,2.f));
                EZ[ijk] *= exp(-alpha*pow(beta-k,2.f));
            }
        }
    }
    #ifdef _OPENMP
//    #pragma omp for
    #endif
    for(i=ix-beta;i<ix;i++){
        for(j=0;j<iy;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                EX[ijk] *= exp(-alpha*pow(beta-(ix-i)+1,2.f));
                EY[ijk] *= exp(-alpha*pow(beta-(ix-i)+1,2.f));
                EZ[ijk] *= exp(-alpha*pow(beta-(ix-i)+1,2.f));
            }
        }
    }
    #ifdef _OPENMP
//    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=iy-beta;j<iy;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                EX[ijk] *= exp(-alpha*pow(beta-(iy-j)+1,2.f));
                EY[ijk] *= exp(-alpha*pow(beta-(iy-j)+1,2.f));
                EZ[ijk] *= exp(-alpha*pow(beta-(iy-j)+1,2.f));
            }
        }
    }
    #ifdef _OPENMP
//    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=0;j<iy;j++){
            for(k=iz-beta;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                EX[ijk] *= exp(-alpha*pow(beta-(iz-k)+1,2.f));
                EY[ijk] *= exp(-alpha*pow(beta-(iz-k)+1,2.f));
                EZ[ijk] *= exp(-alpha*pow(beta-(iz-k)+1,2.f));
            }
        }
    }
}
}
//////////////////////////////////////////
void cerjan_h(double *HX,double *HY,double *HZ)
{
    int i,j,k,ijk;
    int beta = 20;
    double alpha = 5.5e-2;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(i=0;i<beta;i++){
        for(j=0;j<iy;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                HX[ijk] *= exp(-alpha*pow(beta-i,2.f));
                HY[ijk] *= exp(-alpha*pow(beta-i,2.f));
                HZ[ijk] *= exp(-alpha*pow(beta-i,2.f));
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=0;j<beta;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                HX[ijk] *= exp(-alpha*pow(beta-j,2.f));
                HY[ijk] *= exp(-alpha*pow(beta-j,2.f));
                HZ[ijk] *= exp(-alpha*pow(beta-j,2.f));
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=0;j<iy;j++){
            for(k=0;k<beta;k++){
                ijk = k*ix*iy + j*ix + i;
                HX[ijk] *= exp(-alpha*pow(beta-k,2.f));
                HY[ijk] *= exp(-alpha*pow(beta-k,2.f));
                HZ[ijk] *= exp(-alpha*pow(beta-k,2.f));
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(i=ix-beta;i<ix;i++){
        for(j=0;j<iy;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                HX[ijk] *= exp(-alpha*pow(beta-(ix-i)+1,2.f));
                HY[ijk] *= exp(-alpha*pow(beta-(ix-i)+1,2.f));
                HZ[ijk] *= exp(-alpha*pow(beta-(ix-i)+1,2.f));
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=iy-beta;j<iy;j++){
            for(k=0;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                HX[ijk] *= exp(-alpha*pow(beta-(iy-j)+1,2.f));
                HY[ijk] *= exp(-alpha*pow(beta-(iy-j)+1,2.f));
                HZ[ijk] *= exp(-alpha*pow(beta-(iy-j)+1,2.f));
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(i=0;i<ix;i++){
        for(j=0;j<iy;j++){
            for(k=iz-beta;k<iz;k++){
                ijk = k*ix*iy + j*ix + i;
                HX[ijk] *= exp(-alpha*pow(beta-(iz-k)+1,2.f));
                HY[ijk] *= exp(-alpha*pow(beta-(iz-k)+1,2.f));
                HZ[ijk] *= exp(-alpha*pow(beta-(iz-k)+1,2.f));
            }
        }
    }
}
}
//////////////////////////////////////////
void cerjan2_e(double *EX,double *EY,double *EZ)
{
    int i,j,k,n,ijk,jx,jy,jz;
    int beta = 20;
    double alpha = 1.0e-3;
    FILE *otita,*otita2;
    otita=fopen("./data/otita.dat","w");
    otita2=fopen("./data/otita2.dat","w");

    #ifdef _OPENMP
    #pragma omp parallel for private(jx,jy,jz)
    #endif
    for(i=0;i<ix*iy*iz;i++){
        jz = i / (ix*iy);
        jy = (i-ix*iy*jz)/ix;
        jx = i-ix*iy*jz - jy*ix;

        if(jx>=0 && jx<beta && \
           jy>=0 && jy<iy   && \
           jz>=0 && jz<iz   ){
             EX[i] *= exp(-alpha*pow((double)(beta-jx),2.f));
             EY[i] *= exp(-alpha*pow((double)(beta-jx),2.f));
             EZ[i] *= exp(-alpha*pow((double)(beta-jx),2.f));
//        fprintf(otita,"%03d  %03d  %03d   %14.3e\n", jx,jy,jz,exp(-alpha*pow((double)(beta-jx),2.f)));
        }
        if(jx>=0 && jx<ix   && \
           jy>=0 && jy<beta && \
           jz>=0 && jz<iz   ){
             EX[i] *= exp(-alpha*pow((double)(beta-jy),2.f));
             EY[i] *= exp(-alpha*pow((double)(beta-jy),2.f));
             EZ[i] *= exp(-alpha*pow((double)(beta-jy),2.f));
        }
        if(jx>=0 && jx<ix   && \
           jy>=0 && jy<iy   && \
           jz>=0 && jz<beta ){
             EX[i] *= exp(-alpha*pow((double)(beta-jz),2.f));
             EY[i] *= exp(-alpha*pow((double)(beta-jz),2.f));
             EZ[i] *= exp(-alpha*pow((double)(beta-jz),2.f));
        }
        if(jx>=ix-beta && jx<ix && \
           jy>=0       && jy<iy && \
           jz>=0       && jz<iz ){
             EX[i] *= exp(-alpha*pow((double)(beta-(ix-jx)+1),2.f));
             EY[i] *= exp(-alpha*pow((double)(beta-(ix-jx)+1),2.f));
             EZ[i] *= exp(-alpha*pow((double)(beta-(ix-jx)+1),2.f));
//        fprintf(otita2,"%03d  %03d  %03d   %14.3e\n", jx,jy,jz,exp(-alpha*pow(beta-(ix-jx)+1,2.f)) );
        }
        if(jx>=0       && jx<ix && \
           jy>=iy-beta && jy<iy && \
           jz>=0       && jz<iz ){
             EX[i] *= exp(-alpha*pow((double)(beta-(iy-jy)+1),2.f));
             EY[i] *= exp(-alpha*pow((double)(beta-(iy-jy)+1),2.f));
             EZ[i] *= exp(-alpha*pow((double)(beta-(iy-jy)+1),2.f));
        }
        if(jx>=0       && jx<ix && \
           jy>=0       && jy<iy && \
           jz>=iz-beta && jz<iz ){
             EX[i] *= exp(-alpha*pow((double)(beta-(iz-jz)+1),2.f));
             EY[i] *= exp(-alpha*pow((double)(beta-(iz-jz)+1),2.f));
             EZ[i] *= exp(-alpha*pow((double)(beta-(iz-jz)+1),2.f));
        }
    }
}
//////////////////////////////////////////
void cerjan2_h(double *EX,double *EY,double *EZ)
{
    int i,j,k,n,ijk,jx,jy,jz;
    int beta = 20;
    double alpha = 1.0e-3;

    #ifdef _OPENMP
    #pragma omp parallel for private(jx,jy,jz)
    #endif
    for(i=0;i<ix*iy*iz;i++){
        jz = i / (ix*iy);
        jy = (i-ix*iy*jz)/ix;
        jx = i-ix*iy*jz - jy*ix;

        if(jx>=0 && jx<beta && \
           jy>=0 && jy<iy   && \
           jz>=0 && jz<iz   ){
             HX[i] *= exp(-alpha*pow((double)(beta-jx),2.f));
             HY[i] *= exp(-alpha*pow((double)(beta-jx),2.f));
             HZ[i] *= exp(-alpha*pow((double)(beta-jx),2.f));
        }
        if(jx>=0 && jx<ix   && \
           jy>=0 && jy<beta && \
           jz>=0 && jz<iz   ){
             HX[i] *= exp(-alpha*pow((double)(beta-jy),2.f));
             HY[i] *= exp(-alpha*pow((double)(beta-jy),2.f));
             HZ[i] *= exp(-alpha*pow((double)(beta-jy),2.f));
        }
        if(jx>=0 && jx<ix   && \
           jy>=0 && jy<iy   && \
           jz>=0 && jz<beta ){
             HX[i] *= exp(-alpha*pow((double)(beta-jz),2.f));
             HY[i] *= exp(-alpha*pow((double)(beta-jz),2.f));
             HZ[i] *= exp(-alpha*pow((double)(beta-jz),2.f));
        }
        if(jx>=ix-beta && jx<ix && \
           jy>=0       && jy<iy && \
           jz>=0       && jz<iz ){
             HX[i] *= exp(-alpha*pow((double)(beta-(ix-jx)+1),2.f));
             HY[i] *= exp(-alpha*pow((double)(beta-(ix-jx)+1),2.f));
             HZ[i] *= exp(-alpha*pow((double)(beta-(ix-jx)+1),2.f));
        }
        if(jx>=0       && jx<ix && \
           jy>=iy-beta && jy<iy && \
           jz>=0       && jz<iz ){
             HX[i] *= exp(-alpha*pow((double)(beta-(iy-jy)+1),2.f));
             HY[i] *= exp(-alpha*pow((double)(beta-(iy-jy)+1),2.f));
             HZ[i] *= exp(-alpha*pow((double)(beta-(iy-jy)+1),2.f));
        }
        if(jx>=0       && jx<ix && \
           jy>=0       && jy<iy && \
           jz>=iz-beta && jz<iz ){
             HX[i] *= exp(-alpha*pow((double)(beta-(iz-jz)+1),2.f));
             HY[i] *= exp(-alpha*pow((double)(beta-(iz-jz)+1),2.f));
             HZ[i] *= exp(-alpha*pow((double)(beta-(iz-jz)+1),2.f));
        }
    }
}
