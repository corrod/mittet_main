#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fwi3d_cpu_alloex1_v1.4.h"
#include "fwi3d_cpu_cpmlabc_v1.4.h"
//////////////////////////////////////////
void init_cpml(){
    int i;
    //double eps2 = 1.f/2.f/omega_0;
    // thickness of the PML layer in grid points
    double delta;
    // Refrection coefficient
    double Rcoef=0.01f; // R should be [10^-2, 10^-12]
    // compute esig_max
    // Rickard &  Georgieva(2003) Problem-independent enhancement of PML ABC ....
    // (n=1, beta=2)
    double nn=3.f;  // nn should be [2,6]
    double order = 0.f; // order should be (0,3]
    double ma = 1.f;
    double esig_max,msig_max;
    double kappa_max = 1.f;  //kappa should be [0,10]
    double aex_max   = 0.2f; //if aex_max = 0.f, CPML changes to UPML
    double optToMax  = 10.0;

    FILE *fpmlx1,*fpmly1,*fpmlz1;
    fpmlx1 = fopen("./data/pml_conductivity_x","w");
    fpmly1 = fopen("./data/pml_conductivity_y","w");
    fpmlz1 = fopen("./data/pml_conductivity_z","w");
    printf("ncpml    : %d\n",ncpml);
    printf("order    : %lf\n",order);

    delta = ncpml*dx;
    esig_max = (nn+order+1.f)*cmax*log(1.f/Rcoef) / (2.f*delta) * optToMax;
    //msig_max = esig_max*MU0/eps2;
    printf("esig_max : %e\n",esig_max);
    //printf("msig_max : %e\n",msig_max);
    fprintf(fpmlx1," i     esig            msig            ekap               mkap    \
            aex           amx           be_x           ce_x           kedx\n");
    for(i=0;i<ix;i++){
      if(i<ncpml){
        esigx[i]   = esig_max * pow((double)(ncpml - i      )/(double)(ncpml-1.f), nn+order);
        msigx[i]   = esig_max * pow((double)(ncpml - i -0.5f)/(double)(ncpml-1.f), nn+order);
        ekappax[i] = 1.f + (kappa_max -1.f)*pow((double)(ncpml -i      )/(double)(ncpml-1.f), nn);
        mkappax[i] = 1.f + (kappa_max -1.f)*pow((double)(ncpml -i -0.5f)/(double)(ncpml-1.f), nn);
        aex[i]     = aex_max * pow((double)(i      )/(double)(ncpml-1.f), ma);
        amx[i]     = aex_max * pow((double)(i +0.5f)/(double)(ncpml-1.f), ma);

        be_x[i]    = exp(-(esigx[i] / ekappax[i] + aex[i] ) * dt);
        bh_x[i]    = exp(-(msigx[i] / mkappax[i] + amx[i] ) * dt);
        ce_x[i]    = esigx[i] * (be_x[i] - 1.f) / ( esigx[i] + ekappax[i]*aex[i]) / ekappax[i];
        ch_x[i]    = msigx[i] * (bh_x[i] - 1.f) / ( msigx[i] + mkappax[i]*amx[i]) / mkappax[i];
        kedx[i]    = ekappax[i] * dx;
        khdx[i]    = mkappax[i] * dx;
      }else if(i > (ix-ncpml-1) ){
        esigx[i]   = esig_max * pow((double)(i -ix +1.f  + ncpml)/(double)(ncpml-1.f), nn+order);
        msigx[i]   = esig_max * pow((double)(i -ix +0.5f + ncpml)/(double)(ncpml-1.f), nn+order);
        ekappax[i] = 1.f + (kappa_max -1.f)*pow((double)(i -ix +1.0f +ncpml)/(double)(ncpml-1.f), nn);
        mkappax[i] = 1.f + (kappa_max -1.f)*pow((double)(i -ix +0.5f +ncpml)/(double)(ncpml-1.f), nn);
        aex[i]     = aex_max * pow((double)(-i +ix - 1.0f)/(double)(ncpml-1.f), ma);
        amx[i]     = aex_max * pow((double)(-i +ix - 0.5f)/(double)(ncpml-1.f), ma);

        be_x[i]    = exp(-(esigx[i] / ekappax[i] + aex[i] ) * dt);
        bh_x[i]    = exp(-(msigx[i] / mkappax[i] + amx[i] ) * dt);
        ce_x[i]    = esigx[i] * (be_x[i] - 1.f) / ( esigx[i] + ekappax[i]*aex[i]) / ekappax[i];
        ch_x[i]    = msigx[i] * (bh_x[i] - 1.f) / ( msigx[i] + mkappax[i]*amx[i]) / mkappax[i];
        kedx[i]    = ekappax[i] * dx;
        khdx[i]    = mkappax[i] * dx;
      }else{
        esigx[i]   = msigx[i]    = 0.f;
        ekappax[i] = mkappax[i]  = 1.f;
        aex[i]     = amx[i]      = 0.f;
        be_x[i]    = ce_x[i]     = 0.f;
        bh_x[i]    = ch_x[i]     = 0.f;
        kedx[i]    = ekappax[i] * dx;
        khdx[i]    = mkappax[i] * dx;
      }
      fprintf(fpmlx1,"%03d    %e    %e    %e    %e    %e    %e    %e    %e    %e\n", \
                       i,esigx[i],msigx[i],ekappax[i],mkappax[i],aex[i],amx[i],be_x[i], ce_x[i], kedx[i]);
    }

    fprintf(fpmly1," i     esig            msig            ekap               mkap    \
            aex           amx           be_x           ce_x           kedx\n");
    delta = ncpml*dy;
    esig_max  = (nn+order+1.f)*cmax*log(1.f/Rcoef) / (2.f*delta);
    for(i=0;i<iy;i++){
      if(i<ncpml){
        esigy[i]   = esig_max * pow((double)(ncpml - i      )/(double)(ncpml-1.f), nn+order);
        msigy[i]   = esig_max * pow((double)(ncpml - i -0.5f)/(double)(ncpml-1.f), nn+order);
        ekappay[i] = 1.f + (kappa_max -1.f)*pow((double)(ncpml -i      )/(double)(ncpml-1.f), nn);
        mkappay[i] = 1.f + (kappa_max -1.f)*pow((double)(ncpml -i -0.5f)/(double)(ncpml-1.f), nn);
        aey[i]     = aex_max * pow((double)(i      )/(double)(ncpml-1.f), ma);
        amy[i]     = aex_max * pow((double)(i +0.5f)/(double)(ncpml-1.f), ma);

        be_y[i]    = exp(-(esigy[i] / ekappay[i] + aey[i] ) * dt);
        bh_y[i]    = exp(-(msigy[i] / mkappay[i] + amy[i] ) * dt);
        ce_y[i]    = esigy[i] * (be_y[i] - 1.f) / ( esigy[i] + ekappay[i]*aey[i]) / ekappay[i];
        ch_y[i]    = msigy[i] * (bh_y[i] - 1.f) / ( msigy[i] + mkappay[i]*amy[i]) / mkappay[i];
        kedy[i]    = ekappay[i] * dy;
        khdy[i]    = mkappay[i] * dy;
      }else if(i > (iy-ncpml-1) ){
        esigy[i]   = esig_max * pow((double)(i -iy +1.f  + ncpml)/(double)(ncpml-1.f), nn+order);
        msigy[i]   = esig_max * pow((double)(i -iy +0.5f + ncpml)/(double)(ncpml-1.f), nn+order);
        ekappay[i] = 1.f + (kappa_max -1.f)*pow((double)(i -iy +1.0f +ncpml)/(double)(ncpml-1.f), nn);
        mkappay[i] = 1.f + (kappa_max -1.f)*pow((double)(i -iy +0.5f +ncpml)/(double)(ncpml-1.f), nn);
        aey[i]     = aex_max * pow((double)(-i +iy - 1.0f)/(double)(ncpml-1.f), ma);
        amy[i]     = aex_max * pow((double)(-i +iy - 0.5f)/(double)(ncpml-1.f), ma);

        be_y[i]    = exp(-(esigy[i] / ekappay[i] + aey[i] ) * dt);
        bh_y[i]    = exp(-(msigy[i] / mkappay[i] + amy[i] ) * dt);
        ce_y[i]    = esigy[i] * (be_y[i] - 1.f) / ( esigy[i] + ekappay[i]*aey[i]) / ekappay[i];
        ch_y[i]    = msigy[i] * (bh_y[i] - 1.f) / ( msigy[i] + mkappay[i]*amy[i]) / mkappay[i];
        kedy[i]    = ekappay[i] * dy;
        khdy[i]    = mkappay[i] * dy;
      }else{
        esigy[i]   = msigy[i]    = 0.f;
        ekappay[i] = mkappay[i]  = 1.f;
        aey[i]     = amy[i]      = 0.f;
        be_y[i]    = ce_y[i]     = 0.f;
        bh_y[i]    = ch_y[i]     = 0.f;
        kedy[i]    = ekappay[i] * dy;
        khdy[i]    = mkappay[i] * dy;
      }
      fprintf(fpmly1,"%03d    %e    %e    %e    %e    %e    %e    %e    %e    %e\n", \
                       i,esigy[i],msigy[i],ekappay[i],mkappay[i],aey[i],amy[i],be_y[i], ce_y[i], kedy[i]);
    }

    fprintf(fpmlz1," i     esig            msig            ekap               mkap    \
            aex           amx           be_x           ce_x           kedx\n");
    delta = ncpml*dz;
    esig_max  = (nn+order+1.f)*cmax*log(1.f/Rcoef) / (2.f*delta);
    for(i=0;i<iz;i++){
      if(i<ncpml){
        esigz[i]   = esig_max * pow((double)(ncpml - i      )/(double)(ncpml-1.f), nn+order);
        msigz[i]   = esig_max * pow((double)(ncpml - i -0.5f)/(double)(ncpml-1.f), nn+order);
        ekappaz[i] = 1.f + (kappa_max -1.f)*pow((double)(ncpml -i      )/(double)(ncpml-1.f), nn);
        mkappaz[i] = 1.f + (kappa_max -1.f)*pow((double)(ncpml -i -0.5f)/(double)(ncpml-1.f), nn);
        aez[i]     = aex_max * pow((double)(i      )/(double)(ncpml-1.f), ma);
        amz[i]     = aex_max * pow((double)(i +0.5f)/(double)(ncpml-1.f), ma);

        be_z[i]    = exp(-(esigz[i] / ekappaz[i] + aez[i] ) * dt);
        bh_z[i]    = exp(-(msigz[i] / mkappaz[i] + amz[i] ) * dt);
        ce_z[i]    = esigz[i] * (be_z[i] - 1.f) / ( esigz[i] + ekappaz[i]*aez[i]) / ekappaz[i];
        ch_z[i]    = msigz[i] * (bh_z[i] - 1.f) / ( msigz[i] + mkappaz[i]*amz[i]) / mkappaz[i];
        kedz[i]    = ekappaz[i] * dz;
        khdz[i]    = mkappaz[i] * dz;
      }else if(i > (iz-ncpml-1) ){
        esigz[i]   = esig_max * pow((double)(i -iz +1.f  + ncpml)/(double)(ncpml-1.f), nn+order);
        msigz[i]   = esig_max * pow((double)(i -iz +0.5f + ncpml)/(double)(ncpml-1.f), nn+order);
        ekappaz[i] = 1.f + (kappa_max -1.f)*pow((double)(i -iz +1.0f +ncpml)/(double)(ncpml-1.f), nn);
        mkappaz[i] = 1.f + (kappa_max -1.f)*pow((double)(i -iz +0.5f +ncpml)/(double)(ncpml-1.f), nn);
        aez[i]     = aex_max * pow((double)(-i +iz - 1.0f)/(double)(ncpml-1.f), ma);
        amz[i]     = aex_max * pow((double)(-i +iz - 0.5f)/(double)(ncpml-1.f), ma);

        be_z[i]    = exp(-(esigz[i] / ekappaz[i] + aez[i] ) * dt);
        bh_z[i]    = exp(-(msigz[i] / mkappaz[i] + amz[i] ) * dt);
        ce_z[i]    = esigz[i] * (be_z[i] - 1.f) / ( esigz[i] + ekappaz[i]*aez[i]) / ekappaz[i];
        ch_z[i]    = msigz[i] * (bh_z[i] - 1.f) / ( msigz[i] + mkappaz[i]*amz[i]) / mkappaz[i];
        kedz[i]    = ekappaz[i] * dz;
        khdz[i]    = mkappaz[i] * dz;
      }else{
        esigz[i]   = msigz[i]    = 0.f;
        ekappaz[i] = mkappaz[i]  = 1.f;
        aez[i]     = amz[i]      = 0.f;
        be_z[i]    = ce_z[i]     = 0.f;
        bh_z[i]    = ch_z[i]     = 0.f;
        kedz[i]    = ekappaz[i] * dz;
        khdz[i]    = mkappaz[i] * dz;
      }
      fprintf(fpmlz1,"%03d    %e    %e    %e    %e    %e    %e    %e    %e    %e\n", \
                       i,esigz[i],msigz[i],ekappaz[i],mkappaz[i],aez[i],amz[i],be_z[i], ce_z[i], kedz[i]);
    }
    fclose(fpmlx1);
    fclose(fpmly1);
    fclose(fpmlz1);
}
//////////////////////////////////////////
void epml_abcs(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){

    xe_pml(EX, EY, EZ, HX, HY, HZ);
    ye_pml(EX, EY, EZ, HX, HY, HZ);
    ze_pml(EX, EY, EZ, HX, HY, HZ);
}
//////////////////////////////////////////
void hpml_abcs(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){

    xh_pml(EX, EY, EZ, HX, HY, HZ);
    yh_pml(EX, EY, EZ, HX, HY, HZ);
    zh_pml(EX, EY, EZ, HX, HY, HZ);
}
//////////////////////////////////////////
void epml_abcs4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){

    xe_pml4(EX, EY, EZ, HX, HY, HZ);
    ye_pml4(EX, EY, EZ, HX, HY, HZ);
    ze_pml4(EX, EY, EZ, HX, HY, HZ);
}
//////////////////////////////////////////
void hpml_abcs4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){

    xh_pml4(EX, EY, EZ, HX, HY, HZ);
    yh_pml4(EX, EY, EZ, HX, HY, HZ);
    zh_pml4(EX, EY, EZ, HX, HY, HZ);
}
//////////////////////////////////////////
void xe_pml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;

    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=1;i<ncpml;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Ezx1[ijk] = be_x[i] * psi_Ezx1[ijk] \
                        + ce_x[i] * (HY[ijk] - HY[ijk - 1] ) /dx;
          psi_Eyx1[ijk] = be_x[i] * psi_Eyx1[ijk] \
                        + ce_x[i] * (HZ[ijk] - HZ[ijk - 1] ) /dx;
          EZ[ijk] += cb_z[ijk] * psi_Ezx1[ijk];
          EY[ijk] -= cb_y[ijk] * psi_Eyx1[ijk];
        }
      }
    }

    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=ix-ncpml+1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Ezx1[ijk] = be_x[i] * psi_Ezx1[ijk] \
                        + ce_x[i] * (HY[ijk] - HY[ijk - 1] ) /dx;
          psi_Eyx1[ijk] = be_x[i] * psi_Eyx1[ijk] \
                        + ce_x[i] * (HZ[ijk] - HZ[ijk - 1] ) /dx;
          EZ[ijk] += cb_z[ijk] * psi_Ezx1[ijk];
          EY[ijk] -= cb_y[ijk] * psi_Eyx1[ijk];
        }
      }
    }

}
//////////////////////////////////////////
void ye_pml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;

    for(k=0;k<iz-1;k++){
      for(j=1;j<ncpml;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exy1[ijk] = be_y[j] * psi_Exy1[ijk] \
                        + ce_y[j] * (HZ[ijk] - HZ[ijk - ix] ) /dy;
          psi_Ezy1[ijk] = be_y[j] * psi_Ezy1[ijk] \
                        + ce_y[j] * (HX[ijk] - HX[ijk - ix] ) /dy;
          EX[ijk] += cb_x[ijk] * psi_Exy1[ijk];
          EZ[ijk] -= cb_z[ijk] * psi_Ezy1[ijk];
        }
      }
    }

    for(k=0;k<iz-1;k++){
      for(j=iy-ncpml+1;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exy1[ijk] = be_y[j] * psi_Exy1[ijk] \
                        + ce_y[j] * (HZ[ijk] - HZ[ijk - ix] ) /dy;
          psi_Ezy1[ijk] = be_y[j] * psi_Ezy1[ijk] \
                        + ce_y[j] * (HX[ijk] - HX[ijk - ix] ) /dy;
          EX[ijk] += cb_x[ijk] * psi_Exy1[ijk];
          EZ[ijk] -= cb_z[ijk] * psi_Ezy1[ijk];
        }
      }
    }

}
//////////////////////////////////////////
void ze_pml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;

    for(k=1;k<ncpml;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exz1[ijk] = be_z[k] * psi_Exz1[ijk] \
                        + ce_z[k] * (HY[ijk] - HY[ijk - ix*iy] ) /dz;
          psi_Eyz1[ijk] = be_z[k] * psi_Eyz1[ijk] \
                        + ce_z[k] * (HX[ijk] - HX[ijk - ix*iy] ) /dz;
          EX[ijk] -= cb_x[ijk] * psi_Exz1[ijk];
          EY[ijk] += cb_y[ijk] * psi_Eyz1[ijk];
        }
      }
    }

    for(k=iz-ncpml+1;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exz1[ijk] = be_z[k] * psi_Exz1[ijk] \
                        + ce_z[k] * (HY[ijk] - HY[ijk - ix*iy] ) /dz;
          psi_Eyz1[ijk] = be_z[k] * psi_Eyz1[ijk] \
                        + ce_z[k] * (HX[ijk] - HX[ijk - ix*iy] ) /dz;
          EX[ijk] -= cb_x[ijk] * psi_Exz1[ijk];
          EY[ijk] += cb_y[ijk] * psi_Eyz1[ijk];
        }
      }
    }

}
//////////////////////////////////////////
void xh_pml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;

    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ncpml-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hyx1[ijk] = bh_x[i] * psi_Hyx1[ijk] \
                        + ch_x[i] * (EZ[ijk + 1] - EZ[ijk] ) /dx;
          psi_Hzx1[ijk] = bh_x[i] * psi_Hzx1[ijk] \
                        + ch_x[i] * (EY[ijk + 1] - EY[ijk] ) /dx;
          HY[ijk] += db_y[ijk] * psi_Hyx1[ijk];
          HZ[ijk] -= db_z[ijk] * psi_Hzx1[ijk];
        }
      }
    }

    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=ix-ncpml;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hyx1[ijk] = bh_x[i] * psi_Hyx1[ijk] \
                        + ch_x[i] * (EZ[ijk + 1] - EZ[ijk] ) /dx;
          psi_Hzx1[ijk] = bh_x[i] * psi_Hzx1[ijk] \
                        + ch_x[i] * (EY[ijk + 1] - EY[ijk] ) /dx;
          HY[ijk] += db_y[ijk] * psi_Hyx1[ijk];
          HZ[ijk] -= db_z[ijk] * psi_Hzx1[ijk];
        }
      }
    }
}
//////////////////////////////////////////
void yh_pml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;

    for(k=0;k<iz-1;k++){
      for(j=0;j<ncpml-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxy1[ijk] = bh_y[j] * psi_Hxy1[ijk] \
                        + ch_y[j] * (EZ[ijk + ix] - EZ[ijk] ) /dy;
          psi_Hzy1[ijk] = bh_y[j] * psi_Hzy1[ijk] \
                        + ch_y[j] * (EX[ijk + ix] - EX[ijk] ) /dy;
          HX[ijk] -= db_x[ijk] * psi_Hxy1[ijk];
          HZ[ijk] += db_z[ijk] * psi_Hzy1[ijk];
        }
      }
    }

    for(k=0;k<iz-1;k++){
      for(j=iy-ncpml;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxy1[ijk] = bh_y[j] * psi_Hxy1[ijk] \
                        + ch_y[j] * (EZ[ijk + ix] - EZ[ijk] ) /dy;
          psi_Hzy1[ijk] = bh_y[j] * psi_Hzy1[ijk] \
                        + ch_y[j] * (EX[ijk + ix] - EX[ijk] ) /dy;
          HX[ijk] -= db_x[ijk] * psi_Hxy1[ijk];
          HZ[ijk] += db_z[ijk] * psi_Hzy1[ijk];
        }
      }
    }

}
//////////////////////////////////////////
void zh_pml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;

    for(k=0;k<ncpml-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxz1[ijk] = bh_z[k] * psi_Hxz1[ijk] \
                        + ch_z[k] * (EY[ijk + ix*iy] - EY[ijk] ) /dz;
          psi_Hyz1[ijk] = bh_z[k] * psi_Hyz1[ijk] \
                        + ch_z[k] * (EX[ijk + ix*iy] - EX[ijk] ) /dz;
          HX[ijk] += db_x[ijk] * psi_Hxz1[ijk];
          HY[ijk] -= db_y[ijk] * psi_Hyz1[ijk];
        }
      }
    }

    for(k=iz-ncpml;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxz1[ijk] = bh_z[k] * psi_Hxz1[ijk] \
                        + ch_z[k] * (EY[ijk + ix*iy] - EY[ijk] ) /dz;
          psi_Hyz1[ijk] = bh_z[k] * psi_Hyz1[ijk] \
                        + ch_z[k] * (EX[ijk + ix*iy] - EX[ijk] ) /dz;
          HX[ijk] += db_x[ijk] * psi_Hxz1[ijk];
          HY[ijk] -= db_y[ijk] * psi_Hyz1[ijk];
        }
      }
    }

}
//////////////////////////////////////////
void xe_pml4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=2;i<ncpml-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Ezx1[ijk] = be_x[i] * psi_Ezx1[ijk] \
                        + ce_x[i] * (c1*HY[ijk] - c1*HY[ijk - 1] +c2*HY[ijk +1]-c2*HY[ijk-2]) /dx;
          psi_Eyx1[ijk] = be_x[i] * psi_Eyx1[ijk] \
                        + ce_x[i] * (c1*HZ[ijk] - c1*HZ[ijk - 1] +c2*HZ[ijk+ 1]-c2*HZ[ijk-2]) /dx;
          EZ[ijk] += cb_z[ijk] * psi_Ezx1[ijk];
          EY[ijk] -= cb_y[ijk] * psi_Eyx1[ijk];
        }
      }
    }

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=ix-ncpml+1;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Ezx1[ijk] = be_x[i] * psi_Ezx1[ijk] \
                        + ce_x[i] * (c1*HY[ijk] - c1*HY[ijk - 1] +c2*HY[ijk +1]-c2*HY[ijk-2]) /dx;
          psi_Eyx1[ijk] = be_x[i] * psi_Eyx1[ijk] \
                        + ce_x[i] * (c1*HZ[ijk] - c1*HZ[ijk - 1] +c2*HZ[ijk+ 1]-c2*HZ[ijk-2]) /dx;
          EZ[ijk] += cb_z[ijk] * psi_Ezx1[ijk];
          EY[ijk] -= cb_y[ijk] * psi_Eyx1[ijk];
        }
      }
    }
}

}
//////////////////////////////////////////
void ye_pml4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=2;j<ncpml;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exy1[ijk] = be_y[j] * psi_Exy1[ijk] \
                        + ce_y[j] * (c1*HZ[ijk] - c1*HZ[ijk - ix] +c2*HZ[ijk+ix]-c2*HZ[ijk-2*ix]) /dy;
          psi_Ezy1[ijk] = be_y[j] * psi_Ezy1[ijk] \
                        + ce_y[j] * (c1*HX[ijk] - c1*HX[ijk - ix] +c2*HX[ijk+ix]-c2*HX[ijk-2*ix]) /dy;
          EX[ijk] += cb_x[ijk] * psi_Exy1[ijk];
          EZ[ijk] -= cb_z[ijk] * psi_Ezy1[ijk];
        }
      }
    }

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=iy-ncpml+1;j<iy-2;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exy1[ijk] = be_y[j] * psi_Exy1[ijk] \
                        + ce_y[j] * (c1*HZ[ijk] - c1*HZ[ijk - ix] +c2*HZ[ijk+ix]-c2*HZ[ijk-2*ix]) /dy;
          psi_Ezy1[ijk] = be_y[j] * psi_Ezy1[ijk] \
                        + ce_y[j] * (c1*HX[ijk] - c1*HX[ijk - ix] +c2*HX[ijk+ix]-c2*HX[ijk-2*ix]) /dy;
          EX[ijk] += cb_x[ijk] * psi_Exy1[ijk];
          EZ[ijk] -= cb_z[ijk] * psi_Ezy1[ijk];
        }
      }
    }
}

}
//////////////////////////////////////////
void ze_pml4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<ncpml;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exz1[ijk] = be_z[k] * psi_Exz1[ijk] \
                        + ce_z[k] * (c1*HY[ijk] - c1*HY[ijk - ix*iy] +c2*HY[ijk+ix*iy] - c2*HY[ijk-2*ix*iy]) /dz;
          psi_Eyz1[ijk] = be_z[k] * psi_Eyz1[ijk] \
                        + ce_z[k] * (c1*HX[ijk] - c1*HX[ijk - ix*iy] +c2*HX[ijk+ix*iy] - c2*HX[ijk-2*ix*iy]) /dz;
          EX[ijk] -= cb_x[ijk] * psi_Exz1[ijk];
          EY[ijk] += cb_y[ijk] * psi_Eyz1[ijk];
        }
      }
    }
}

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=iz-ncpml+1;k<iz-2;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Exz1[ijk] = be_z[k] * psi_Exz1[ijk] \
                        + ce_z[k] * (c1*HY[ijk] - c1*HY[ijk - ix*iy] +c2*HY[ijk+ix*iy] - c2*HY[ijk-2*ix*iy]) /dz;
          psi_Eyz1[ijk] = be_z[k] * psi_Eyz1[ijk] \
                        + ce_z[k] * (c1*HX[ijk] - c1*HX[ijk - ix*iy] +c2*HX[ijk+ix*iy] - c2*HX[ijk-2*ix*iy]) /dz;
          EX[ijk] -= cb_x[ijk] * psi_Exz1[ijk];
          EY[ijk] += cb_y[ijk] * psi_Eyz1[ijk];
        }
      }
    }
}

}
//////////////////////////////////////////
void xh_pml4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=1;i<ncpml-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hyx1[ijk] = bh_x[i] * psi_Hyx1[ijk] \
                        + ch_x[i] * (c1*EZ[ijk + 1] - c1*EZ[ijk] +c2*EZ[ijk+2]-c2*EZ[ijk-1]) /dx;
          psi_Hzx1[ijk] = bh_x[i] * psi_Hzx1[ijk] \
                        + ch_x[i] * (c1*EY[ijk + 1] - c1*EY[ijk] +c2*EY[ijk+2]-c2*EY[ijk-1]) /dx;
          HY[ijk] += db_y[ijk] * psi_Hyx1[ijk];
          HZ[ijk] -= db_z[ijk] * psi_Hzx1[ijk];
        }
      }
    }

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=ix-ncpml;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hyx1[ijk] = bh_x[i] * psi_Hyx1[ijk] \
                        + ch_x[i] * (c1*EZ[ijk + 1] - c1*EZ[ijk] +c2*EZ[ijk+2]-c2*EZ[ijk-1]) /dx;
          psi_Hzx1[ijk] = bh_x[i] * psi_Hzx1[ijk] \
                        + ch_x[i] * (c1*EY[ijk + 1] - c1*EY[ijk] +c2*EY[ijk+2]-c2*EY[ijk-1]) /dx;
          HY[ijk] += db_y[ijk] * psi_Hyx1[ijk];
          HZ[ijk] -= db_z[ijk] * psi_Hzx1[ijk];
        }
      }
    }
}
}
//////////////////////////////////////////
void yh_pml4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
//    double c1 = 1.14443f, c2=-0.04886f;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=1;j<ncpml-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxy1[ijk] = bh_y[j] * psi_Hxy1[ijk] \
                        + ch_y[j] * (c1*EZ[ijk + ix] - c1*EZ[ijk] +c2*EZ[ijk+2*ix] -c2*EZ[ijk-ix]) /dy;
          psi_Hzy1[ijk] = bh_y[j] * psi_Hzy1[ijk] \
                        + ch_y[j] * (c1*EX[ijk + ix] - c1*EX[ijk] +c2*EX[ijk+2*ix] -c2*EX[ijk-ix]) /dy;
          HX[ijk] -= db_x[ijk] * psi_Hxy1[ijk];
          HZ[ijk] += db_z[ijk] * psi_Hzy1[ijk];
        }
      }
    }

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz-1;k++){
      for(j=iy-ncpml;j<iy-2;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxy1[ijk] = bh_y[j] * psi_Hxy1[ijk] \
                        + ch_y[j] * (c1*EZ[ijk + ix] - c1*EZ[ijk] +c2*EZ[ijk+2*ix] -c2*EZ[ijk-ix]) /dy;
          psi_Hzy1[ijk] = bh_y[j] * psi_Hzy1[ijk] \
                        + ch_y[j] * (c1*EX[ijk + ix] - c1*EX[ijk] +c2*EX[ijk+2*ix] -c2*EX[ijk-ix]) /dy;
          HX[ijk] -= db_x[ijk] * psi_Hxy1[ijk];
          HZ[ijk] += db_z[ijk] * psi_Hzy1[ijk];
        }
      }
    }
}

}
//////////////////////////////////////////
void zh_pml4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ){
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<ncpml-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxz1[ijk] = bh_z[k] * psi_Hxz1[ijk] \
                        + ch_z[k] * (c1*EY[ijk + ix*iy] - c1*EY[ijk] +c2*EY[ijk+2*ix*iy]-c2*EY[ijk-ix*iy]) /dz;
          psi_Hyz1[ijk] = bh_z[k] * psi_Hyz1[ijk] \
                        + ch_z[k] * (c1*EX[ijk + ix*iy] - c1*EX[ijk] +c2*EX[ijk+2*ix*iy]-c2*EX[ijk-ix*iy]) /dz;
          HX[ijk] += db_x[ijk] * psi_Hxz1[ijk];
          HY[ijk] -= db_y[ijk] * psi_Hyz1[ijk];
        }
      }
    }

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=iz-ncpml;k<iz-2;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          psi_Hxz1[ijk] = bh_z[k] * psi_Hxz1[ijk] \
                        + ch_z[k] * (c1*EY[ijk + ix*iy] - c1*EY[ijk] +c2*EY[ijk+2*ix*iy]-c2*EY[ijk-ix*iy]) /dz;
          psi_Hyz1[ijk] = bh_z[k] * psi_Hyz1[ijk] \
                        + ch_z[k] * (c1*EX[ijk + ix*iy] - c1*EX[ijk] +c2*EX[ijk+2*ix*iy]-c2*EX[ijk-ix*iy]) /dz;
          HX[ijk] += db_x[ijk] * psi_Hxz1[ijk];
          HY[ijk] -= db_y[ijk] * psi_Hyz1[ijk];
        }
      }
    }
}
}
