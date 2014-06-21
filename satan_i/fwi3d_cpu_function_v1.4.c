#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "fftw3.h"
#pragma comment( lib, "fftw3.lib")
#include "fwi3d_cpu_function_v1.4.h"
#include "fwi3d_cpu_alloex1_v1.4.h"
#include "fwi3d_cpu_alloex2_v1.4.h"
#include "fwi3d_cpu_cpmlabc_v1.4.h"
#include "fwi3d_cpu_lattice_v1.4.h"

///////////////////////////////////////////////
void read_shotrec()
{
  int i,itmp;
  double dtmp;
  FILE *ifp1,*ifp2,*ifp4;
  ifp1 = fopen("./data/model_env.dat","r");
  ifp2 = fopen("./data/minmax.dat","r");
  ifp4 = fopen("./data/iseabed.dat","r");

  fscanf(ifp1,"%d",&ix);
  fscanf(ifp1,"%d",&iy);
  fscanf(ifp1,"%d",&iz);
  fscanf(ifp1,"%lf",&dx);
  fscanf(ifp1,"%lf",&dy);
  fscanf(ifp1,"%lf",&dz);
  fscanf(ifp1,"%lf",&dt_ratio);
  fscanf(ifp1,"%d",&it);
  fscanf(ifp1,"%d",&mstep);
  fscanf(ifp1,"%d",&trstep);
  fscanf(ifp1,"%lf",&omega_0);
  fscanf(ifp1,"%lf",&fmax_w);
  fscanf(ifp1,"%lf",&Glim);
  fscanf(ifp1,"%d",&ncpml);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%d",&startF);
  fscanf(ifp1,"%d",&upperF);
  fscanf(ifp1,"%d",&shot_num);

  shot_px=(int *)malloc(sizeof(int) * shot_num);
  shot_py=(int *)malloc(sizeof(int) * shot_num);
  shot_pz=(int *)malloc(sizeof(int) * shot_num);

  for(i=0;i<shot_num;i++){
      fscanf(ifp1, "%d  %d  %d", &shot_px[i], &shot_py[i], &shot_pz[i]);
      printf("Transmitter : %d, ( %d, %d, %d)\n", i,shot_px[i],shot_py[i],shot_pz[i]);
  }

  fscanf(ifp1,"%d",&rec_num);

  rec_px=(int *)malloc(sizeof(int) * rec_num);
  rec_py=(int *)malloc(sizeof(int) * rec_num);
  rec_pz=(int *)malloc(sizeof(int) * rec_num);

  for(i=0;i<rec_num;i++){
      fscanf(ifp1, "%d  %d  %d", &rec_px[i], &rec_py[i], &rec_pz[i]);
      printf("Receiver : %d, ( %d, %d, %d)\n", i,rec_px[i],rec_py[i],rec_pz[i]);
  }

  fscanf(ifp2,"%lf", &minval);
  fscanf(ifp2,"%lf", &maxval);
  fscanf(ifp4,"%d",&iseabed);

  fclose(ifp1);
  fclose(ifp2);
  fclose(ifp4);

  printf("shot number: %d,  rec number: %d\n",shot_num,rec_num);
}

///////////////////////////////////////////////
void read_trawave()
{
    int i;
    double tmpt;
    FILE *ift1,*ift2,*ift3;
    ift1 = fopen("./data/waveformX.dat","r");
    ift2 = fopen("./data/waveformY.dat","r");
    ift3 = fopen("./data/waveformZ.dat","r");

    signalX = (double *)malloc(sizeof(double) * it);
    signalY = (double *)malloc(sizeof(double) * it);
    signalZ = (double *)malloc(sizeof(double) * it);

    for(i=0;i<it;i++){
        fscanf(ift1,"%lf %lf",&tmpt,&signalX[i]);
        fscanf(ift2,"%lf %lf",&tmpt,&signalY[i]);
        fscanf(ift3,"%lf %lf",&tmpt,&signalZ[i]);
        JX_f[i] = signalX[i];
    }
    fclose(ift1);
    fclose(ift2);
    fclose(ift3);

}
///////////////////////////////////////////////
void read_waveform()
{
    int i;
    double tmpt;
    FILE *ift1,*ift2,*ift3;
    ift1 = fopen("./data/true_waveformX.dat","r");

    signalTrue = (double *)malloc(sizeof(double) * it);

    for(i=0;i<it;i++){
        fscanf(ift1,"%lf %lf",&tmpt,&signalTrue[i]);
    }
    fclose(ift1);

}
//////////////////////////////////////////////
void init_eh_field_3d()
{
// number of lattice
    n = ix*iy*iz;
// memory allocation
    // for field
    EX  = (double *)malloc(sizeof(double)*n);
    EY  = (double *)malloc(sizeof(double)*n);
    EZ  = (double *)malloc(sizeof(double)*n);
    HX  = (double *)malloc(sizeof(double)*n);
    HY  = (double *)malloc(sizeof(double)*n);
    HZ  = (double *)malloc(sizeof(double)*n);
    // for media id
    id  = (int *)malloc(sizeof(int)*n);
    sig = (double *)malloc(sizeof(double)*n);
    sig2= (double *)malloc(sizeof(double)*n);
    mu  = (double *)malloc(sizeof(double)*n);
    vel = (double *)malloc(sizeof(double)*n);

    // e_coeff1
    //cex = (double *)malloc(sizeof(double)*n);
    //cey = (double *)malloc(sizeof(double)*n);
    //cez = (double *)malloc(sizeof(double)*n);
    // e_coeff2
    //cexry = (double *)malloc(sizeof(double)*n);
    //cexrz = (double *)malloc(sizeof(double)*n);
    //ceyrx = (double *)malloc(sizeof(double)*n);
    //ceyrz = (double *)malloc(sizeof(double)*n);
    //cezrx = (double *)malloc(sizeof(double)*n);
    //cezry = (double *)malloc(sizeof(double)*n);
    // h_coeff
    //chxry = (double *)malloc(sizeof(double)*n);
    //chxrz = (double *)malloc(sizeof(double)*n);
    //chyrx = (double *)malloc(sizeof(double)*n);
    //chyrz = (double *)malloc(sizeof(double)*n);
    //chzrx = (double *)malloc(sizeof(double)*n);
    //chzry = (double *)malloc(sizeof(double)*n);
    // pml
    esigx = (double *)malloc(sizeof(double)*ix);
    esigy = (double *)malloc(sizeof(double)*iy);
    esigz = (double *)malloc(sizeof(double)*iz);
    msigx = (double *)malloc(sizeof(double)*ix);
    msigy = (double *)malloc(sizeof(double)*iy);
    msigz = (double *)malloc(sizeof(double)*iz);
    //exy = (double *)malloc(sizeof(double)*ix*iy*iz);
    //exz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //eyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //eyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //ezx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //ezy = (double *)malloc(sizeof(double)*ix*iy*iz);
    //hxy = (double *)malloc(sizeof(double)*ix*iy*iz);
    //hxz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //hyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //hyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //hzx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //hzy = (double *)malloc(sizeof(double)*ix*iy*iz);
    //cexy = (double *)malloc(sizeof(double)*ix*iy*iz);
    //cexz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //ceyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //ceyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //cezx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //cezy = (double *)malloc(sizeof(double)*ix*iy*iz);
    //chxy = (double *)malloc(sizeof(double)*ix*iy*iz);
    //chxz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //chyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //chyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    //chzx = (double *)malloc(sizeof(double)*ix*iy*iz);
    //chzy = (double *)malloc(sizeof(double)*ix*iy*iz);
    // cpml
    ca_x = (double *)malloc(sizeof(double)*ix*iy*iz);
    ca_y = (double *)malloc(sizeof(double)*ix*iy*iz);
    ca_z = (double *)malloc(sizeof(double)*ix*iy*iz);
    cb_x = (double *)malloc(sizeof(double)*ix*iy*iz);
    cb_y = (double *)malloc(sizeof(double)*ix*iy*iz);
    cb_z = (double *)malloc(sizeof(double)*ix*iy*iz);
    da_x = (double *)malloc(sizeof(double)*ix*iy*iz);
    da_y = (double *)malloc(sizeof(double)*ix*iy*iz);
    da_z = (double *)malloc(sizeof(double)*ix*iy*iz);
    db_x = (double *)malloc(sizeof(double)*ix*iy*iz);
    db_y = (double *)malloc(sizeof(double)*ix*iy*iz);
    db_z = (double *)malloc(sizeof(double)*ix*iy*iz);
    kedx = (double *)malloc(sizeof(double)*ix);
    kedy = (double *)malloc(sizeof(double)*iy);
    kedz = (double *)malloc(sizeof(double)*iz);
    khdx = (double *)malloc(sizeof(double)*ix);
    khdy = (double *)malloc(sizeof(double)*iy);
    khdz = (double *)malloc(sizeof(double)*iz);
    psi_Eyx1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Ezx1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Exy1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Ezy1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Exz1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Eyz1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Hyx1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Hzx1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Hxy1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Hzy1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Hxz1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    psi_Hyz1 = (double *)malloc(sizeof(double)*ix*iy*iz);
    be_x = (double *)malloc(sizeof(double)*ix);
    ce_x = (double *)malloc(sizeof(double)*ix);
    be_y = (double *)malloc(sizeof(double)*iy);
    ce_y = (double *)malloc(sizeof(double)*iy);
    be_z = (double *)malloc(sizeof(double)*iz);
    ce_z = (double *)malloc(sizeof(double)*iz);
    bh_x = (double *)malloc(sizeof(double)*ix);
    ch_x = (double *)malloc(sizeof(double)*ix);
    bh_y = (double *)malloc(sizeof(double)*iy);
    ch_y = (double *)malloc(sizeof(double)*iy);
    bh_z = (double *)malloc(sizeof(double)*iz);
    ch_z = (double *)malloc(sizeof(double)*iz);
    ekappax = (double *)malloc(sizeof(double)*ix);
    ekappay = (double *)malloc(sizeof(double)*iy);
    ekappaz = (double *)malloc(sizeof(double)*iz);
    mkappax = (double *)malloc(sizeof(double)*ix);
    mkappay = (double *)malloc(sizeof(double)*iy);
    mkappaz = (double *)malloc(sizeof(double)*iz);
    aex = (double *)malloc(sizeof(double)*ix);
    aey = (double *)malloc(sizeof(double)*iy);
    aez = (double *)malloc(sizeof(double)*iz);
    amx = (double *)malloc(sizeof(double)*ix);
    amy = (double *)malloc(sizeof(double)*iy);
    amz = (double *)malloc(sizeof(double)*iz);
    // Laplace transform
    EX_f = (double *)malloc(sizeof(double)*it);
    EX_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    EY_f = (double *)malloc(sizeof(double)*it);
    EY_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    EZ_f = (double *)malloc(sizeof(double)*it);
    EZ_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    JX_f = (double _Complex *)malloc(sizeof(double _Complex)*it);
    JX_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    JY_f = (double _Complex *)malloc(sizeof(double _Complex)*it);
    JY_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    JZ_f = (double _Complex *)malloc(sizeof(double _Complex)*it);
    JZ_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    GX_t = (double _Complex *)malloc(sizeof(double _Complex)*it);
    GX_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    GY_t = (double _Complex *)malloc(sizeof(double _Complex)*it);
    GY_w = (double _Complex *)malloc(sizeof(double _Complex)*it);
    GZ_t = (double _Complex *)malloc(sizeof(double _Complex)*it);
    GZ_w = (double _Complex *)malloc(sizeof(double _Complex)*it);

    all_recEX = (double *)malloc(sizeof(double)*it*rec_num);
    all_recEY = (double *)malloc(sizeof(double)*it*rec_num);
    all_recEZ = (double *)malloc(sizeof(double)*it*rec_num);

}
///////////////////////////////////////////////
void true_model()
{
    int i,j,k,ijk;
    double itmp1,itmp2,itmp3;
    FILE *ifm;
    ifm = fopen("./data/conductivity.dat","r");

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk = k*ix*iy + j*ix + i;
          fscanf(ifm,"%lf   %lf   %lf   %lf\n", &itmp1, &itmp2, &itmp3, &sig[ijk]);
        }
      }
    }

    sigmax = sigmin = sig[0];
    for(i=0;i<n;i++){
      if(sig[i] > sigmax)  sigmax = sig[i];
      if(sig[i] < sigmin)  sigmin = sig[i];
    }
    printf("\n");
    printf("sigmax     = %lf\n",sigmax);
    printf("sigmin     = %lf\n",sigmin);
    fclose(ifm);

}
///////////////////////////////////////////////
void init_model()
{
    int i,j,k,ijk;
    double itmp1,itmp2,itmp3;
    FILE *ifm;
    ifm = fopen("./data/conductivity.dat","r");

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk = k*ix*iy + j*ix + i;
          fscanf(ifm,"%lf   %lf   %lf   %lf\n", &itmp1, &itmp2, &itmp3, &sig[ijk]);
          if(sig[ijk] >=4.0f) sig[ijk] = 0.45f;
        }
      }
    }

    sigmax = sigmin = sig[0];
    for(i=0;i<n;i++){
      if(sig[i] > sigmax)  sigmax = sig[i];
      if(sig[i] < sigmin)  sigmin = sig[i];
    }
    printf("\n");
    printf("sigmax     = %lf\n",sigmax);
    printf("sigmin     = %lf\n",sigmin);
    fclose(ifm);

}
///////////////////////////////////////////////
void media_coeff_3d()
{
    int i,j,k,ijk;
    double eps2;
    FILE *ofm3,*ofm4,*ofm5;

    ofm3 = fopen("./data/coef1.dat","w");
    ofm4 = fopen("./data/coef2.dat","w");
    ofm5 = fopen("./data/coef3.dat","w");
    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk  = k*iy*ix + j*ix + i;
          eps2 = sig[ijk] /2.f /omega_0;
          //cex[ijk]   = eps[id[ijk]]/(eps[id[ijk]]+sig[id[ijk]]*dt);
          //cey[ijk]   = eps[id[ijk]]/(eps[id[ijk]]+sig[id[ijk]]*dt);
          //cez[ijk]   = eps[id[ijk]]/(eps[id[ijk]]+sig[id[ijk]]*dt);
          ////// dt/eps/(1+(sig*dt)/(2*eps)) / dl
          //cexry[ijk] = dt/(eps[id[ijk]]+sig[id[ijk]]*dt)/dy;
          //cexrz[ijk] = dt/(eps[id[ijk]]+sig[id[ijk]]*dt)/dz;
          //ceyrx[ijk] = dt/(eps[id[ijk]]+sig[id[ijk]]*dt)/dx;
          //ceyrz[ijk] = dt/(eps[id[ijk]]+sig[id[ijk]]*dt)/dz;
          //cezrx[ijk] = dt/(eps[id[ijk]]+sig[id[ijk]]*dt)/dx;
          //cezry[ijk] = dt/(eps[id[ijk]]+sig[id[ijk]]*dt)/dy;
          ////// dt/mu /dl
          //chxry[ijk] = (dt/mu[id[ijk]])/dy;
          //chxrz[ijk] = (dt/mu[id[ijk]])/dz;
          //chyrx[ijk] = (dt/mu[id[ijk]])/dx;
          //chyrz[ijk] = (dt/mu[id[ijk]])/dz;
          //chzrx[ijk] = (dt/mu[id[ijk]])/dx;
          //chzry[ijk] = (dt/mu[id[ijk]])/dy;

          // Coefficient refering Mittet 2010
          //cexy[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
          //            / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));
          //cexz[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
          //            / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
          //ceyx[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
          //            / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
          //ceyz[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
          //            / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
          //cezx[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
          //            / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
          //cezy[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
          //            / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));

          //chxy[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*MU0))) \
          //            / (1.0f + ((msigy[j]*dt)/(2.f*MU0)));
          //chxz[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*MU0))) \
          //            / (1.0f + ((msigz[k]*dt)/(2.f*MU0)));
          //chyx[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*MU0))) \
          //            / (1.0f + ((msigx[i]*dt)/(2.f*MU0)));
          //chyz[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*MU0))) \
          //            / (1.0f + ((msigz[k]*dt)/(2.f*MU0)));
          //chzx[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*MU0))) \
          //            / (1.0f + ((msigx[i]*dt)/(2.f*MU0)));
          //chzy[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*MU0))) \
          //            / (1.0f + ((msigy[j]*dt)/(2.f*MU0)));

          //cex[ijk]   = 1.f;
          //cey[ijk]   = 1.f;
          //cez[ijk]   = 1.f;
          // dt/eps/(1+(sig*dt)/(2*eps)) / dl
          //cexry[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2)) /dy;
          //cexrz[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2)) /dz;
          //ceyrx[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2)) /dx;
          //ceyrz[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2)) /dz;
          //cezrx[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2)) /dx;
          //cezry[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2)) /dy;
          // dt/mu /dl
          //chxry[ijk] = dt/MU0 /(1.f+(msigy[j]*dt)/(2.f*MU0)) /dy;
          //chxrz[ijk] = dt/MU0 /(1.f+(msigz[k]*dt)/(2.f*MU0)) /dz;
          //chyrx[ijk] = dt/MU0 /(1.f+(msigx[i]*dt)/(2.f*MU0)) /dx;
          //chyrz[ijk] = dt/MU0 /(1.f+(msigz[k]*dt)/(2.f*MU0)) /dz;
          //chzrx[ijk] = dt/MU0 /(1.f+(msigx[i]*dt)/(2.f*MU0)) /dx;
          //chzry[ijk] = dt/MU0 /(1.f+(msigy[j]*dt)/(2.f*MU0)) /dy;
          // CPML Coefficient
          // coef1
          ca_x[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
          ca_y[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));
          ca_z[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
          // coef2
          da_x[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*eps2))) \
                      / (1.0f + ((msigx[i]*dt)/(2.f*eps2)));
          da_y[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*eps2))) \
                      / (1.0f + ((msigy[j]*dt)/(2.f*eps2)));
          da_z[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*eps2))) \
                      / (1.0f + ((msigz[k]*dt)/(2.f*eps2)));
          // coef3
          cb_x[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2));
          cb_y[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2));
          cb_z[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2));
          // coef4
          db_x[ijk] = dt/MU0 /(1.f+(msigx[i]*dt)/(2.f*eps2));
          db_y[ijk] = dt/MU0 /(1.f+(msigy[j]*dt)/(2.f*eps2));
          db_z[ijk] = dt/MU0 /(1.f+(msigz[k]*dt)/(2.f*eps2));

          fprintf(ofm3,"%+8e  %+8e  %+8e  %+8e  %+8e  %+8e\n",ca_x[ijk],ca_y[ijk],ca_z[ijk],cb_x[ijk],cb_y[ijk],cb_z[ijk]);
          fprintf(ofm4,"%+8e  %+8e  %+8e  %+8e  %+8e  %+8e\n",da_x[ijk],da_y[ijk],da_z[ijk],db_x[ijk],db_y[ijk],db_z[ijk]);
          fprintf(ofm5,"%8e  %8e  %8e  %8e  %8e  %8e\n",esigx[i],esigy[j],esigz[k],msigx[i],msigy[j],msigz[k]);
        }
      }
    }
    fclose(ofm3);
    fclose(ofm4);
    fclose(ofm5);

}
///////////////////////////////////////////////
void init_FILE(int isource)
{
  int i;
  char strx[256],stry[256],strz[256];

  ofe1 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofe2 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofe3 = (FILE **)malloc(rec_num*sizeof(FILE *));

  for(i=0;i<rec_num;i++){
    sprintf(strx,"./data/forward_x_%03d_%03d.dat",isource,i);
    sprintf(stry,"./data/forward_y_%03d_%03d.dat",isource,i);
    sprintf(strz,"./data/forward_z_%03d_%03d.dat",isource,i);
    ofe1[i] = fopen(strx,"w");
    ofe2[i] = fopen(stry,"w");
    ofe3[i] = fopen(strz,"w");
  }
}

///////////////////////////////////////////////
void close_FILE(void)
{
    fclose(*ofe1); fclose(*ofe2); fclose(*ofe3);
}
///////////////////////////////////////////////
void close_FILE2(void)
{
    fclose(*ofe1); fclose(*ofe2); fclose(*ofe3);
}
///////////////////////////////////////////////
void close_FILE3(void)
{
    fclose(*ofe1); fclose(*ofe2); fclose(*ofe3);
}
///////////////////////////////////////////////
void e_field(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
//  ex
    for(k=1;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
        //for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EX[ijk]=cex[ijk]*EX[ijk] \
                + cexry[ijk]*(HZ[ijk] - HZ[ijk-ix]   )\
                - cexrz[ijk]*(HY[ijk] - HY[ijk-ix*iy]);
        }
      }
    }
//  ey
    for(k=1;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
      //for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EY[ijk]=cey[ijk]*EY[ijk] \
                + ceyrz[ijk]*(HX[ijk] - HX[ijk-ix*iy]) \
                - ceyrx[ijk]*(HZ[ijk] - HZ[ijk-1]    );
        }
      }
    }
//  ez
    for(k=0;k<iz-1;k++){
    //for(k=0;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EZ[ijk]=cez[ijk]*EZ[ijk] \
                + cezrx[ijk]*(HY[ijk] - HY[ijk-1] ) \
                - cezry[ijk]*(HX[ijk] - HX[ijk-ix]);
        }
      }
    }
}

///////////////////////////////////////////////
void read_source_3d(double *EX_r,double *signalX_r,int isource,int step)
{
    int i,j,k;
    int jx,jy,jz,ijk,sourcep;
    double r;

    jx=shot_px[isource];
    jy=shot_py[isource];
    jz=shot_pz[isource];

    sourcep=ix*iy*(jz-1) + ix*(jy-1) + jx-1;

    EX_r[sourcep] -= (signalX_r[step] +signalX_r[step+1])/2.f* \
                     dt*2.f*omega_0/sig[sourcep]/dx/dy/dz;

}
///////////////////////////////////////////////
void h_field(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
//  hx
    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
        //for(i=0;i<ix;i++){
          ijk = k*ix*iy + j*ix + i;
          HX[ijk]=HX[ijk] \
                - chxry[ijk]*(EZ[ijk+ix] - EZ[ijk]   ) \
                + chxrz[ijk]*(EY[ijk+ix*iy] - EY[ijk]);
        }
      }
    }
//  hy
    for(k=0;k<iz-1;k++){
     for(j=0;j<iy-1;j++){
     //for(j=0;j<iy;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HY[ijk]=HY[ijk] \
                - chyrz[ijk]*(EX[ijk+ix*iy] - EX[ijk]) \
                + chyrx[ijk]*(EZ[ijk+1] - EZ[ijk] );
        }
      }
    }
//  hz
    for(k=0;k<iz-1;k++){
    //for(k=0;k<iz;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HZ[ijk]=HZ[ijk] \
                - chzrx[ijk]*(EY[ijk+1]  - EY[ijk]) \
                + chzry[ijk]*(EX[ijk+ix] - EX[ijk]);
        }
     }
    }
}
///////////////////////////////////////////////
void e_field4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    //double c1 = 1.14443f, c2=-0.04886f;
    double c1 = 1.125f, c2=-0.04167f;
//  ex
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-2;k++){
      for(j=2;j<iy-2;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EX[ijk]=cex[ijk]*EX[ijk] \
                + cexry[ijk]*(c1*HZ[ijk] - c1*HZ[ijk-ix]    + c2*HZ[ijk+ix]    - c2*HZ[ijk-2*ix]) \
                - cexrz[ijk]*(c1*HY[ijk] - c1*HY[ijk-ix*iy] + c2*HY[ijk+ix*iy] - c2*HY[ijk-2*ix*iy]);
        }
      }
    }
//
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-2;k++){
      for(j=1;j<iy-1;j++){
        for(i=2;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          EY[ijk]=cey[ijk]*EY[ijk] \
                + ceyrz[ijk]*(c1*HX[ijk] - c1*HX[ijk-ix*iy] + c2*HX[ijk+ix*iy] -c2*HX[ijk-2*ix*iy]) \
                - ceyrx[ijk]*(c1*HZ[ijk] - c1*HZ[ijk-1]     + c2*HZ[ijk+1]     -c2*HZ[ijk-2]);
        }
      }
    }
//  ez
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-1;k++){
      for(j=2;j<iy-2;j++){
        for(i=2;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          EZ[ijk]=cez[ijk]*EZ[ijk] \
                + cezrx[ijk]*(c1*HY[ijk] - c1*HY[ijk-1]  + c2*HY[ijk+1]  -c2*HY[ijk-2]) \
                - cezry[ijk]*(c1*HX[ijk] - c1*HX[ijk-ix] + c2*HX[ijk+ix] -c2*HX[ijk-2*ix]);
        }
      }
    }
}
}

///////////////////////////////////////////////
void h_field4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    //double c1 = 1.14443f, c2=-0.04886f;
    double c1 = 1.125f, c2=-0.04167f;
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
//  hx
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-2;k++){
      for(j=2;j<iy-2;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HX[ijk]=HX[ijk] \
                - chxry[ijk]*(c1*EZ[ijk+ix] - c1*EZ[ijk]    + c2*EZ[ijk+2*ix] - c2*EZ[ijk-ix]) \
                + chxrz[ijk]*(c1*EY[ijk+ix*iy] - c1*EY[ijk] + c2*EY[ijk+2*ix*iy] - c2*EY[ijk-ix*iy]);
        }
      }
    }
//  hy
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-2;k++){
      for(j=1;j<iy-1;j++){
        for(i=2;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HY[ijk]=HY[ijk] \
                - chyrz[ijk]*(c1*EX[ijk+ix*iy] - c1*EX[ijk] + c2*EX[ijk+2*ix*iy] - c2*EX[ijk-ix*iy]) \
                + chyrx[ijk]*(c1*EZ[ijk+1] - c1*EZ[ijk] + c2*EZ[ijk+2] - c2*EZ[ijk-1]);
        }
      }
    }
//  hz
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-1;k++){
      for(j=2;j<iy-2;j++){
        for(i=2;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HZ[ijk]=HZ[ijk] \
                - chzrx[ijk]*(c1*EY[ijk+1]  - c1*EY[ijk] + c2*EY[ijk+2] - c2*EY[ijk-1]) \
                + chzry[ijk]*(c1*EX[ijk+ix] - c1*EX[ijk] + c2*EX[ijk+2*ix] - c2*EX[ijk-ix]);
        }
      }
    }
}
}
///////////////////////////////////////////////
void e_field42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk,jx,jy,jz;
    //double c1 = 1.14443f, c2=-0.04886f;
    double c1 = 1.125f, c2=-0.04167f;
//  ex
    #ifdef _OPENMP
    #pragma omp parallel for private(jx,jy,jz)
    #endif
    for(i=0;i<ix*iy*iz;i++){
        jz = (int)(i / (ix*iy));
        jy = (int)((i-ix*iy*jz) / ix);
        jx = (int)(i-ix*iy*jz - jy*ix);
// ---
        if(jz>=2 && jz<iz-2 && \
           jy>=2 && jy<iy-2 && \
           jx>=1 && jx<ix-1 ){
        //if(jz>=2 && jz<iz-1 && \
        //   jy>=2 && jy<iy-1 && \
           jx>=1 && jx<ix-1 ){
             EX[i]=cex[i]*EX[i] \
                   + cexry[i]*(c1*HZ[i] - c1*HZ[i-ix]    + c2*HZ[i+ix]    - c2*HZ[i-2*ix]) \
                   - cexrz[i]*(c1*HY[i] - c1*HY[i-ix*iy] + c2*HY[i+ix*iy] - c2*HY[i-2*ix*iy]);
        }
// ---
        if(jz>=2 && jz<iz-2 && \
           jy>=1 && jy<iy-1 && \
           jx>=2 && jx<ix-2 ){
             EY[i]=cey[i]*EY[i] \
                   + ceyrz[i]*(c1*HX[i] - c1*HX[i-ix*iy] + c2*HX[i+ix*iy] -c2*HX[i-2*ix*iy]) \
                   - ceyrx[i]*(c1*HZ[i] - c1*HZ[i-1]     + c2*HZ[i+1]     -c2*HZ[i-2]);
        }
// ---
        if(jz>=1 && jz<iz-1 && \
           jy>=2 && jy<iy-2 && \
           jx>=2 && jx<ix-2 ){
             EZ[i]=cez[i]*EZ[i] \
                   + cezrx[i]*(c1*HY[i] - c1*HY[i-1]  + c2*HY[i+1]  -c2*HY[i-2]) \
                   - cezry[i]*(c1*HX[i] - c1*HX[i-ix] + c2*HX[i+ix] -c2*HX[i-2*ix]);
        }
    }
}

///////////////////////////////////////////////
void h_field42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk,jx,jy,jz;
    //double c1 = 1.14443f, c2=-0.04886f;
    double c1 = 1.125f, c2=-0.04167f;

    #ifdef _OPENMP
    #pragma omp parallel for private(jx,jy,jz)
    #endif
    for(i=0;i<ix*iy*iz;i++){
        jz = (int)(i / (ix*iy));
        jy = (int)((i-ix*iy*jz) / ix);
        jx = (int)(i-ix*iy*jz - jy*ix);
//  hx
        if(jz>=2 && jz<iz-2 && \
           jy>=2 && jy<iy-2 && \
           jx>=1 && jx<ix-1 ){
             HX[i]=HX[i] \
                   - chxry[i]*(c1*EZ[i+ix] - c1*EZ[i]    + c2*EZ[i+2*ix] - c2*EZ[i-ix]) \
                   + chxrz[i]*(c1*EY[i+ix*iy] - c1*EY[i] + c2*EY[i+2*ix*iy] - c2*EY[i-ix*iy]);
        }
//  hy
        if(jz>=2 && jz<iz-2 && \
           jy>=1 && jy<iy-1 && \
           jx>=2 && jx<ix-2 ){
             HY[i]=HY[i] \
                   - chyrz[i]*(c1*EX[i+ix*iy] - c1*EX[i] + c2*EX[i+2*ix*iy] - c2*EX[i-ix*iy]) \
                   + chyrx[i]*(c1*EZ[i+1] - c1*EZ[i] + c2*EZ[i+2] - c2*EZ[i-1]);
        }
//  hz
        if(jz>=1 && jz<iz-1 && \
           jy>=2 && jy<iy-2 && \
           jx>=2 && jx<ix-2 ){
             HZ[i]=HZ[i] \
                   - chzrx[i]*(c1*EY[i+1]  - c1*EY[i] + c2*EY[i+2] - c2*EY[i-1]) \
                   + chzry[i]*(c1*EX[i+ix] - c1*EX[i] + c2*EX[i+2*ix] - c2*EX[i-ix]);
        }
}
}
///////////////////////////////////////////////
void output(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ, \
        FILE **ofe1,FILE **ofe2,FILE **ofe3,FILE **ofh1,FILE **ofh2,FILE **ofh3, \
        int step)
{
    int i,j,k,jx,jy,jz,ijk;
    char stt1[256],stt2[256],stt3[256];

    for(i=0;i<rec_num;i++){
        jx=rec_px[i];
        jy=rec_py[i];
        jz=rec_pz[i];
        ijk = ix*iy*(jz-1) + ix*(jy-1) + jx-1;
        all_recEX[i*it + step] = EX[ijk];
        all_recEY[i*it + step] = EY[ijk];
        all_recEZ[i*it + step] = EZ[ijk];
        fprintf(ofe1[i],"%lf    %e\n",(step+0.f)*dt,EX[ijk]);
        fprintf(ofe2[i],"%lf    %e\n",(step+0.f)*dt,EY[ijk]);
        fprintf(ofe3[i],"%lf    %e\n",(step+0.f)*dt,EZ[ijk]);
      }
// output in files at all point
//      if((step%10) == 0){
//        FILE *fptime1,*fptime2,*fptime3;
//        sprintf(stt1, "./data2/time_ex_%05d.dat",step);
//        sprintf(stt2, "./data2/time_ey_%05d.dat",step);
//        sprintf(stt3, "./data2/time_ez_%05d.dat",step);
//        fptime1 = fopen(stt1,"w");
//        fptime2 = fopen(stt2,"w");
//        fptime3 = fopen(stt3,"w");

//        for(k=0;k<iz;k++){
//          for(j=0;j<iy;j++){
//            for(i=0;i<ix;i++){
//              ijk=ix*iy*k + j*ix + i;
//              fprintf(fptime1,"%12lf   %12lf   %12lf   %15.5e\n", i*dx, j*dy,k*dz, EX[ijk]);
//              //fprintf(fptime2,"%12lf   %12lf   %12lf   %15.5e\n", i*dx, j*dy,k*dz, EY[ijk]);
//              //fprintf(fptime3,"%12lf   %12lf   %12lf   %15.5e\n", i*dx, j*dy,k*dz, EZ[ijk]);
//            }
//          }
//        }
//        fclose(fptime1);
//        fclose(fptime2);
//        fclose(fptime3);
//      }
}
///////////////////////////////////////////////
void checking_omp(void)
{
    printf("\n");
    printf("#################   Initializing  ######################\n");
    #ifdef _OPENMP
        printf("#########   Using OpenMP  #########\n");
    #else
        printf("#########   NO OpenMP  #########\n");
    #endif
}

///////////////////////////////////////////////
void init_gradient(void)
{
    int i,n2;
    int direction =3;
    n2 = ix*iy*iz;
    //n2 = (ix-2*ncpml)*(iy-2*ncpml)*(iz-2*ncpml);
    delEx_f   = (double _Complex *)malloc(rec_num*upperF*sizeof(double _Complex));
    delEy_f   = (double _Complex *)malloc(rec_num*upperF*sizeof(double _Complex));
    delEz_f   = (double _Complex *)malloc(rec_num*upperF*sizeof(double _Complex));

    EcalX_b    = (double _Complex *)malloc(sizeof(double _Complex )*n2*it);
    EcalY_b    = (double _Complex *)malloc(sizeof(double _Complex )*n2*it);
    EcalZ_b    = (double _Complex *)malloc(sizeof(double _Complex )*n2*it);
    EcalX_back = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalY_back = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalZ_back = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalX_bb   = (double _Complex *)malloc(sizeof(double _Complex )*n2*rec_num*direction*upperF);
    EcalY_bb   = (double _Complex *)malloc(sizeof(double _Complex )*n2*rec_num*direction*upperF);
    EcalZ_bb   = (double _Complex *)malloc(sizeof(double _Complex )*n2*rec_num*direction*upperF);
    EcalX_ff   = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalY_ff   = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalZ_ff   = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalX_ffL  = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF*shot_num);
    EcalY_ffL  = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF*shot_num);
    EcalZ_ffL  = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF*shot_num);
    EcalX_ff2  = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalY_ff2  = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    EcalZ_ff2  = (double _Complex *)malloc(sizeof(double _Complex )*n2*upperF);
    grad      = (double *)malloc(sizeof(double)*n2);
    grad_p    = (double *)malloc(sizeof(double)*n2);
    Exobs     = (double _Complex *)malloc(sizeof(double _Complex )*rec_num*upperF);
    Eyobs     = (double _Complex *)malloc(sizeof(double _Complex )*rec_num*upperF);
    Ezobs     = (double _Complex *)malloc(sizeof(double _Complex )*rec_num*upperF);
    sense     = (double *)malloc(sizeof(double)*n2);

    err_sum = (double *)malloc(sizeof(double)*MAXITR);

    for(i=0;i<n2;i++)  grad[i] = grad_p[i] = 0.f;
    for(i=0;i<MAXITR;i++) err_sum[i] = 0.f;
}
///////////////////////////////////////////////
void init_FILE2(int isource,int iter)
{
  int i;
  char strx[256],stry[256],strz[256];

  ofe1 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofe2 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofe3 = (FILE **)malloc(rec_num*sizeof(FILE *));

  for(i=0;i<rec_num;i++){
    sprintf(strx,"./data/fex_%03d_%03d_%03d.dat",isource,i,iter);
    sprintf(stry,"./data/fey_%03d_%03d_%03d.dat",isource,i,iter);
    sprintf(strz,"./data/fez_%03d_%03d_%03d.dat",isource,i,iter);
    ofe1[i] = fopen(strx,"w");
    ofe2[i] = fopen(stry,"w");
    ofe3[i] = fopen(strz,"w");
  }
}
///////////////////////////////////////////////
void init_FILE3(int iter)
{
  int inum,idir,irec;
  char strx[256],stry[256],strz[256];

  ofe1 = (FILE **)malloc(3*rec_num*sizeof(FILE *));
  ofe2 = (FILE **)malloc(3*rec_num*sizeof(FILE *));
  ofe3 = (FILE **)malloc(3*rec_num*sizeof(FILE *));

  for(idir=0;idir<3;idir++){
    for(irec=0;irec<rec_num;irec++){
      inum = idir*rec_num + irec;
      sprintf(strx,"./data/rex_%03d_%03d_%03d.dat",idir,irec,iter);
      sprintf(stry,"./data/rey_%03d_%03d_%03d.dat",idir,irec,iter);
      sprintf(strz,"./data/rez_%03d_%03d_%03d.dat",idir,irec,iter);
      ofe1[inum] = fopen(strx,"w");
      ofe2[inum] = fopen(stry,"w");
      ofe3[inum] = fopen(strz,"w");
    }
  }
}
///////////////////////////////////////////////
void read_Eobs(int isource)
{
  int i,j,ij;
  double dtmp, dreal,dimag,dtmp2,dtmp3;
  FILE *ifp_Ex,*ifp_Ey,*ifp_Ez;
  char ch1[256],ch2[256],ch3[256];

  for(i=0; i<rec_num; i++){
      /* storage of observed waveform */
      sprintf(ch1, "./data/spectrum_ex_%03d_%03d.dat",isource,i);
      sprintf(ch2, "./data/spectrum_ey_%03d_%03d.dat",isource,i);
      sprintf(ch3, "./data/spectrum_ez_%03d_%03d.dat",isource,i);
      ifp_Ex = fopen(ch1,"r");
      ifp_Ey = fopen(ch2,"r");
      ifp_Ez = fopen(ch3,"r");
      for(j=0;j<upperF;j++){
        ij = i*upperF + j;
        fscanf(ifp_Ex,"%lf  %lf  %lf  %lf  %lf\n",&dtmp, &dreal, &dimag, &dtmp2,&dtmp3);
        Exobs[ij] = dreal + I*dimag;
        fscanf(ifp_Ey,"%lf  %lf  %lf  %lf  %lf\n",&dtmp, &dreal, &dimag, &dtmp2,&dtmp3);
        Eyobs[ij] = dreal + I*dimag;
        fscanf(ifp_Ez,"%lf  %lf  %lf  %lf  %lf\n",&dtmp, &dreal, &dimag, &dtmp2,&dtmp3);
        Ezobs[ij] = dreal + I*dimag;
      }
      fclose(ifp_Ex);
      fclose(ifp_Ey);
      fclose(ifp_Ez);
  }
}
///////////////////////////////////////////////
void read_3C_Eobs(int isource,int idir)
{
  int i,j,ij;
  double dtmp, dreal,dimag,dtmp2,dtmp3;
  FILE *ifp_Ex,*ifp_Ey,*ifp_Ez;
  char ch1[256],ch2[256],ch3[256];

  for(i=0; i<rec_num; i++){
      /* storage of observed waveform */
      sprintf(ch1, "./data/spectrum_ex_%03d_%03d.dat",isource,i);
      sprintf(ch2, "./data/spectrum_ey_%03d_%03d.dat",isource,i);
      sprintf(ch3, "./data/spectrum_ez_%03d_%03d.dat",isource,i);
      ifp_Ex = fopen(ch1,"r");
      ifp_Ey = fopen(ch2,"r");
      ifp_Ez = fopen(ch3,"r");
      for(j=0;j<upperF;j++){
        ij = i*upperF + j;
        fscanf(ifp_Ex,"%lf  %lf  %lf  %lf  %lf\n",&dtmp, &dreal, &dimag, &dtmp2,&dtmp3);
        Exobs[ij] = dreal + I*dimag;
        fscanf(ifp_Ey,"%lf  %lf  %lf  %lf  %lf\n",&dtmp, &dreal, &dimag, &dtmp2,&dtmp3);
        Eyobs[ij] = dreal + I*dimag;
        fscanf(ifp_Ez,"%lf  %lf  %lf  %lf  %lf\n",&dtmp, &dreal, &dimag, &dtmp2,&dtmp3);
        Ezobs[ij] = dreal + I*dimag;
      }
      fclose(ifp_Ex);
      fclose(ifp_Ey);
      fclose(ifp_Ez);
  }
}
///////////////////////////////////////////////
void set_zero_eh()
{
    int i;
    n=ix*iy*iz;

    for(i=0;i<n;i++){
        EX[i]=EY[i]=EZ[i]=0.f;
        HX[i]=HY[i]=HZ[i]=0.f;
        psi_Eyx1[i] = 0.f;
        psi_Ezx1[i] = 0.f;
        psi_Exy1[i] = 0.f;
        psi_Ezy1[i] = 0.f;
        psi_Exz1[i] = 0.f;
        psi_Eyz1[i] = 0.f;
        psi_Hyx1[i] = 0.f;
        psi_Hzx1[i] = 0.f;
        psi_Hxy1[i] = 0.f;
        psi_Hzy1[i] = 0.f;
        psi_Hxz1[i] = 0.f;
        psi_Hyz1[i] = 0.f;

    }
}
///////////////////////////////////////////////
void e_field4_bp(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{   // ±が反転している
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
//  ex
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-1;k++){
      for(j=2;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EX[ijk]=cex[ijk]*EX[ijk] \
                - cexry[ijk]*(c1*HZ[ijk] - c1*HZ[ijk-ix]    + c2*HZ[ijk+ix]    - c2*HZ[ijk-2*ix]) \
                + cexrz[ijk]*(c1*HY[ijk] - c1*HY[ijk-ix*iy] + c2*HY[ijk+ix*iy] - c2*HY[ijk-2*ix*iy]);
        }
      }
    }
//
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=2;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EY[ijk]=cey[ijk]*EY[ijk] \
                - ceyrz[ijk]*(c1*HX[ijk] - c1*HX[ijk-ix*iy] + c2*HX[ijk+ix*iy] -c2*HX[ijk-2*ix*iy]) \
                + ceyrx[ijk]*(c1*HZ[ijk] - c1*HZ[ijk-1]     + c2*HZ[ijk+1]     -c2*HZ[ijk-2]);
        }
      }
    }
//  ez
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-1;k++){
      for(j=2;j<iy-1;j++){
        for(i=2;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EZ[ijk]=cez[ijk]*EZ[ijk] \
                - cezrx[ijk]*(c1*HY[ijk] - c1*HY[ijk-1]  + c2*HY[ijk+1]  -c2*HY[ijk-2]) \
                + cezry[ijk]*(c1*HX[ijk] - c1*HX[ijk-ix] + c2*HX[ijk+ix] -c2*HX[ijk-2*ix]);
        }
      }
    }
}
}

///////////////////////////////////////////////
void h_field4_bp(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{ // ±が反転している
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
//  hx
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-2;k++){
      for(j=1;j<iy-2;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HX[ijk]=HX[ijk] \
                + chxry[ijk]*(c1*EZ[ijk+ix] - c1*EZ[ijk]    + c2*EZ[ijk+2*ix] - c2*EZ[ijk-ix]) \
                - chxrz[ijk]*(c1*EY[ijk+ix*iy] - c1*EY[ijk] + c2*EY[ijk+2*ix*iy] - c2*EY[ijk-ix*iy]);
        }
      }
    }
//  hy
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-2;k++){
      for(j=1;j<iy-1;j++){
        for(i=1;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HY[ijk]=HY[ijk] \
                + chyrz[ijk]*(c1*EX[ijk+ix*iy] - c1*EX[ijk] + c2*EX[ijk+2*ix*iy] - c2*EX[ijk-ix*iy]) \
                - chyrx[ijk]*(c1*EZ[ijk+1] - c1*EZ[ijk] + c2*EZ[ijk+2] - c2*EZ[ijk-1]);
        }
      }
    }
//  hz
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-1;k++){
      for(j=1;j<iy-2;j++){
        for(i=1;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HZ[ijk]=HZ[ijk] \
                + chzrx[ijk]*(c1*EY[ijk+1]  - c1*EY[ijk] + c2*EY[ijk+2] - c2*EY[ijk-1]) \
                - chzry[ijk]*(c1*EX[ijk+ix] - c1*EX[ijk] + c2*EX[ijk+2*ix] - c2*EX[ijk-ix]);
        }
      }
    }
}
}
///////////////////////////////////////////////
void read_backwave_3d(double *EX_r,double *delEx_r,int step,int irec){
    int i,jx,jy,jz,ijk;

    jx      = rec_px[irec];
    jy      = rec_py[irec];
    jz      = rec_pz[irec];

    ijk     = ix*iy*(jz-1) + ix*(jy-1) + (jx-1);
    EX_r[ijk] -= (delEx_r[step] +delEx_r[step+1])/2.f* \
                     dt*2.f*omega_0/sig[ijk]/dx/dy/dz;

}
///////////////////////////////////////////////
void e_field_cpml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
// ex
    for(k=1;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EX[ijk] = ca_x[ijk] *   EX[ijk] \
                  + cb_x[ijk] * ((HZ[ijk] - HZ[ijk - ix   ]) / kedy[j] \
                               - (HY[ijk] - HY[ijk - ix*iy]) / kedz[k]);
        }
      }
    }
// ey
    for(k=1;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EY[ijk] = ca_y[ijk] *   EY[ijk] \
                  + cb_y[ijk] * ((HX[ijk] - HX[ijk - ix*iy]) / kedz[k] \
                               - (HZ[ijk] - HZ[ijk - 1    ]) / kedx[i]);
        }
      }
    }
// ez
    for(k=0;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EZ[ijk] = ca_z[ijk] *   EZ[ijk] \
                  + cb_z[ijk] * ((HY[ijk] - HY[ijk - 1 ]) / kedx[i] \
                               - (HX[ijk] - HX[ijk - ix]) / kedy[j]);
        }
      }
    }
}
///////////////////////////////////////////////
void h_field_cpml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
// hx
    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HX[ijk] = da_x[ijk] *   HX[ijk] \
                  - db_x[ijk] * ((EZ[ijk + ix   ] - EZ[ijk]) / khdy[j] \
                               - (EY[ijk + ix*iy] - EY[ijk]) / khdz[k]);
        }
      }
    }
// hy
    for(k=0;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HY[ijk] = da_y[ijk] *   HY[ijk] \
                  - db_y[ijk] * ((EX[ijk + ix*iy] - EX[ijk]) / khdz[k] \
                               - (EZ[ijk + 1    ] - EZ[ijk]) / khdx[i]);
        }
      }
    }
// hz
    for(k=1;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HZ[ijk] = da_z[ijk] *   HZ[ijk] \
                  - db_z[ijk] * ((EY[ijk + 1 ] - EY[ijk]) / khdx[i] \
                               - (EX[ijk + ix] - EX[ijk]) / khdy[j]);
        }
      }
    }
}
///////////////////////////////////////////////
void e_field_cpml42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
// ex
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-1;k++){
      for(j=2;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EX[ijk] = ca_x[ijk] * EX[ijk] \
                  + cb_x[ijk] * ((c1*HZ[ijk] - c1*HZ[ijk-ix   ] + c2*HZ[ijk+ix   ] -c2*HZ[ijk-2*ix   ]) / kedy[j] \
                              -  (c1*HY[ijk] - c1*HY[ijk-ix*iy] + c2*HY[ijk+ix*iy] -c2*HY[ijk-2*ix*iy]) / kedz[k]);
        }
      }
    }
// ey
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=2;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=2;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EY[ijk] = ca_y[ijk] *  EY[ijk] \
                  + cb_y[ijk] * ((c1*HX[ijk] - c1*HX[ijk-ix*iy] + c2*HX[ijk+ix*iy] - c2*HX[ijk-2*ix*iy]) / kedz[k] \
                              -  (c1*HZ[ijk] - c1*HZ[ijk-1    ] + c2*HZ[ijk+1    ] - c2*HZ[ijk-2      ]) / kedx[i]);
        }
      }
    }
// ez
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-1;k++){
      for(j=2;j<iy-1;j++){
        for(i=2;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EZ[ijk] = ca_z[ijk] *  EZ[ijk] \
                  + cb_z[ijk] * ((c1*HY[ijk] - c1*HY[ijk-1 ] + c2*HY[ijk+1 ] - c2*HY[ijk-2   ]) / kedx[i] \
                              -  (c1*HX[ijk] - c1*HX[ijk-ix] + c2*HX[ijk+ix] - c2*HX[ijk-2*ix]) / kedy[j]);
        }
      }
    }
}
}
///////////////////////////////////////////////
void h_field_cpml42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk)
#endif
{
// hx
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-2;k++){
      for(j=1;j<iy-2;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HX[ijk] = da_x[ijk] *   HX[ijk] \
                  - db_x[ijk] * ((c1*EZ[ijk+ix   ] - c1*EZ[ijk] + c2*EZ[ijk+2*ix   ] - c2*EZ[ijk-ix   ]) / khdy[j] \
                              -  (c1*EY[ijk+ix*iy] - c1*EY[ijk] + c2*EY[ijk+2*ix*iy] - c2*EY[ijk-ix*iy]) / khdz[k]);
        }
      }
    }
// hy
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-2;k++){
      for(j=1;j<iy-1;j++){
        for(i=1;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HY[ijk] = da_y[ijk] *   HY[ijk] \
                  - db_y[ijk] * ((c1*EX[ijk+ix*iy] - c1*EX[ijk] + c2*EX[ijk+2*ix*iy] - c2*EX[ijk-ix*iy]) / khdz[k] \
                              -  (c1*EZ[ijk+1    ] - c1*EZ[ijk] + c2*EZ[ijk+2      ] - c2*EZ[ijk-1    ]) / khdx[i]);
        }
      }
    }
// hz
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=1;k<iz-1;k++){
      for(j=1;j<iy-2;j++){
        for(i=1;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HZ[ijk] = da_z[ijk] *   HZ[ijk] \
                  - db_z[ijk] * ((c1*EY[ijk+1 ] - c1*EY[ijk] + c2*EY[ijk+2   ] - c2*EY[ijk-1 ]) / khdx[i] \
                              -  (c1*EX[ijk+ix] - c1*EX[ijk] + c2*EX[ijk+2*ix] - c2*EX[ijk-ix]) / khdy[j]);
        }
      }
    }
}
}
///////////////////////////////////////////////
void e_field_cpml42bp(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
 // ±が反転している
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;
// ex
    for(k=2;k<iz-1;k++){
      for(j=2;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EX[ijk] = ca_x[ijk] * EX[ijk] \
                  - cb_x[ijk] * ((c1*HZ[ijk] - c1*HZ[ijk-ix   ] + c2*HZ[ijk+ix   ] -c2*HZ[ijk-2*ix   ]) / kedy[j] \
                              -  (c1*HY[ijk] - c1*HY[ijk-ix*iy] + c2*HY[ijk+ix*iy] -c2*HY[ijk-2*ix*iy]) / kedz[k]);
        }
      }
    }
// ey
    for(k=2;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=2;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EY[ijk] = ca_y[ijk] *  EY[ijk] \
                  - cb_y[ijk] * ((c1*HX[ijk] - c1*HX[ijk-ix*iy] + c2*HX[ijk+ix*iy] - c2*HX[ijk-2*ix*iy]) / kedz[k] \
                              -  (c1*HZ[ijk] - c1*HZ[ijk-1    ] + c2*HZ[ijk+1    ] - c2*HZ[ijk-2      ]) / kedx[i]);
        }
      }
    }
// ez
    for(k=1;k<iz-1;k++){
      for(j=2;j<iy-1;j++){
        for(i=2;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          EZ[ijk] = ca_z[ijk] *  EZ[ijk] \
                  - cb_z[ijk] * ((c1*HY[ijk] - c1*HY[ijk-1 ] + c2*HY[ijk+1 ] - c2*HY[ijk-2   ]) / kedx[i] \
                              -  (c1*HX[ijk] - c1*HX[ijk-ix] + c2*HX[ijk+ix] - c2*HX[ijk-2*ix]) / kedy[j]);
        }
      }
    }
}
///////////////////////////////////////////////
void h_field_cpml42bp(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
 // ±が反転している
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;
// hx
    for(k=1;k<iz-2;k++){
      for(j=1;j<iy-2;j++){
        for(i=1;i<ix-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HX[ijk] = da_x[ijk] *   HX[ijk] \
                  + db_x[ijk] * ((c1*EZ[ijk+ix   ] - c1*EZ[ijk] + c2*EZ[ijk+2*ix   ] - c2*EZ[ijk-ix   ]) / khdy[j] \
                              -  (c1*EY[ijk+ix*iy] - c1*EY[ijk] + c2*EY[ijk+2*ix*iy] - c2*EY[ijk-ix*iy]) / khdz[k]);
        }
      }
    }
// hy
    for(k=1;k<iz-2;k++){
      for(j=1;j<iy-1;j++){
        for(i=1;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HY[ijk] = da_y[ijk] *   HY[ijk] \
                  + db_y[ijk] * ((c1*EX[ijk+ix*iy] - c1*EX[ijk] + c2*EX[ijk+2*ix*iy] - c2*EX[ijk-ix*iy]) / khdz[k] \
                              -  (c1*EZ[ijk+1    ] - c1*EZ[ijk] + c2*EZ[ijk+2      ] - c2*EZ[ijk-1    ]) / khdx[i]);
        }
      }
    }
// hz
    for(k=1;k<iz-1;k++){
      for(j=1;j<iy-2;j++){
        for(i=1;i<ix-2;i++){
          ijk = k*ix*iy + j*ix + i;
          HZ[ijk] = da_z[ijk] *   HZ[ijk] \
                  + db_z[ijk] * ((c1*EY[ijk+1 ] - c1*EY[ijk] + c2*EY[ijk+2   ] - c2*EY[ijk-1 ]) / khdx[i] \
                              -  (c1*EX[ijk+ix] - c1*EX[ijk] + c2*EX[ijk+2*ix] - c2*EX[ijk-ix]) / khdy[j]);
        }
      }
    }
}
///////////////////////////////////////////////
int time_check()
{
    double S=0.6f;
    int Nt;

    Nt=(int)(S*nx/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2))) \
       * sqrt(maxval/minval);

    printf("Nt     :  %d [step]\n",Nt);
    if(Nt >it){
        printf("******* time step tt is violated *******\n");
        exit(1);
    }
    printf("it     :  %d [step] ( <= Your time step)\n", it);

    return 0;
}
///////////////////////////////////////////////
int fmax_check()
{
    double fmax;

    fmax = cmin /Glim /maxin3(dx,dy,dz);

    printf("fmax   :  %lf \n",fmax);
    if(fmax <fmax_w){
        printf("******* fmax is violated *******\n");
        exit(1);
    }
    printf("fmax_w :  %lf ( <= Your fmax_w)\n",fmax_w);
    return 0;
}
///////////////////////////////////////////////
int maxin3(double ddx,double ddy,double ddz)
{
    double tmp;
    if(ddx > ddy){
        tmp = ddx;
    }else{
        tmp = ddy;
    }

    if(tmp > ddz){
        tmp = tmp;
    }else{
        tmp = ddz;
    }
    return tmp;
}
///////////////////////////////////////////////
void start_fdtd()
{
    printf("\n");
    printf("#################   Start FDTD  ######################\n");
}
///////////////////////////////////////////////
void laplace_fft(int isource)
{
  int i,j,k,n,ijk,jx,jy,jz;
  double om,t0,beta;
  size_t mem_size = sizeof(fftw_complex) * it;
  FILE *invGJ, *invGE, *invGG, *invSE, *invSG, *invSJ;
  char str1[256],str2[256],str3[256];
  char sts1[256],sts2[256],sts3[256];

  printf("\n");
  printf("#################   Laplace transform  ######################\n");

  om   = 2.f*M_PI/it /dt;
  t0   = M_PI/fmax_w;
  beta = M_PI*pow(fmax_w,2.f);

  for(i=0;i<rec_num;i++){
    sprintf(sts1,"./data/spectrum_j_%03d_%03d.dat",isource,i);
    sprintf(sts2,"./data/spectrum_e_%03d_%03d.dat",isource,i);
    sprintf(sts3,"./data/spectrum_g_%03d_%03d.dat",isource,i);
    sprintf(str1,"./data/inverse_j_%03d_%03d.dat",isource,i);
    sprintf(str2,"./data/inverse_e_%03d_%03d.dat",isource,i);
    sprintf(str3,"./data/inverse_g_%03d_%03d.dat",isource,i);
    invSJ = fopen(sts1,"w");
    invSE = fopen(sts2,"w");
    invSG = fopen(sts3,"w");
    invGJ = fopen(str1,"w");
    invGE = fopen(str2,"w");
    invGG = fopen(str3,"w");

    fftw_complex *in1  = NULL;
    fftw_complex *in2  = NULL;
    fftw_complex *in3  = NULL;
    fftw_complex *out1 = NULL;
    fftw_complex *out2 = NULL;
    fftw_complex *out3 = NULL;
    fftw_plan p1       = NULL;
    fftw_plan p2       = NULL;
    fftw_plan p3       = NULL;
    for(j=0;j<it;j++){
      EX_f[j] = all_recEX[j + i*it];
    }

    for(n=0;n<it;n++){
      EX_w[n] = JX_w[n] = 0.f;
      for(k=0;k<it;k++){
      //GX_w[n] += EX_f[k] / JX_f[k] *dt;
      //EX_w[n] += EX_f[k] *dt\
      //           *cexp(I*sqrt(omega_0*om*n)*k*dt);
      EX_w[n] += EX_f[k] *dt*\
                 exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
      //JX_w[n] += JX_f[k] * dt *\
                 exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
      JX_w[n] += csqrt(-2.f*omega_0/I/om/n)*JX_f[k] * dt *\
                 exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
      //if(n==1000) printf("%d    %e    %e\n",k,creal(JX_f[k]*exp(-sqrt(omega_0*om*n)*k*dt)),exp(-sqrt(omega_0*om*n)*k*dt));
      }
      JX_w[0] = 2.f*omega_0;
      //JX_w[n] = csqrt(-I*om*n/2.f/omega_0) /JX_w[n];
      //JX_w[n] = 2.f*omega_0 * exp(-sqrt(omega_0*om*n)*t0) * cexp(I*sqrt(omega_0*om*n)*t0)* \
                cexp(-I*omega_0*om*n/2.f/beta);     // J(x,w)
      //GX_w[n] = -I*(I+1.f)*sqrt(om*n*omega_0)*exp(-sqrt(om*n*omega_0)*t0) * \
                cexp(I*sqrt(om*n*omega_0)*t0) * cexp(-I*om*n*omega_0/2.f/beta); // J'(x,w')
      //GX_w[n] = EX_w[n] / JX_w[n];
      GX_w[n] = EX_w[n] / JX_w[n];
    }

    for(n=0;n<it;n++){
      fprintf(invSE,"%10e %10e  %10e  %10e   %10e\n", \
              n*om,creal(EX_w[n]),creal(I*EX_w[n]), creal(EX_w[n])*creal(EX_w[n])+creal(I*EX_w[n])*creal(I*EX_w[n]), atan(creal(EX_w[n]*I)/creal(EX_w[n])));
      fprintf(invSJ,"%10e %10e  %10e  %10e  %10e\n", \
              n*om,creal(JX_w[n]),creal(I*JX_w[n]), creal(JX_w[n])*creal(JX_w[n])+creal(I*JX_w[n])*creal(I*JX_w[n]), atan(creal(JX_w[n]*I)/creal(JX_w[n])));
      fprintf(invSG,"%10e %10e  %10e  %10e  %10e\n", \
            n*om,creal(GX_w[n]),creal(I*GX_w[n]), creal(GX_w[n])*creal(GX_w[n])+creal(I*GX_w[n])*creal(I*GX_w[n]), atan(creal(GX_w[n]*I)/creal(GX_w[n])));
    }

// initialize <fftw3>
    in1  = (fftw_complex*)fftw_malloc( mem_size );
    in2  = (fftw_complex*)fftw_malloc( mem_size );
    in3  = (fftw_complex*)fftw_malloc( mem_size );
    out1 = (fftw_complex*)fftw_malloc( mem_size );
    out2 = (fftw_complex*)fftw_malloc( mem_size );
    out3 = (fftw_complex*)fftw_malloc( mem_size );

    if( !in1 || !out1 ){
      fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
      exit(1);
    }

    p1 = fftw_plan_dft_1d( it, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE );
    p2 = fftw_plan_dft_1d( it, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE );
    p3 = fftw_plan_dft_1d( it, in3, out3, FFTW_FORWARD, FFTW_ESTIMATE );

    for(j=0;j<it;j++){
      in1[j] = JX_w[j];
      in2[j] = EX_w[j];
      in3[j] = GX_w[j];
    }

    fftw_execute(p1);
    fftw_execute(p2);
    fftw_execute(p3);

    for(k=0;k<it;k++){
      out1[k] = out1[k]/it/dt*2.f;
      out2[k] = out2[k]/it/dt*2.f;
      out3[k] = out3[k]/it/dt*2.f;
      fprintf(invGJ,"%e    %e   %e\n",k*dt, creal(out1[k]),  creal(I*out1[k]));
      fprintf(invGE,"%e    %e   %e\n",k*dt, creal(out2[k]),  creal(I*out2[k]));
      fprintf(invGG,"%e    %e   %e\n",k*dt, creal(out3[k]),  creal(I*out3[k]));
      GX_t[k] = out3[k];
    }

    for(j=0;j<it;j++){
      all_recGT[j + i*it] = GX_t[j];
    }

    if(p1  ) fftw_destroy_plan(p1);
    if(p2  ) fftw_destroy_plan(p2);
    if(p3  ) fftw_destroy_plan(p3);
    if(in1 ) fftw_free(in1);
    if(in2 ) fftw_free(in2);
    if(in3 ) fftw_free(in3);
    if(out1) fftw_free(out1);
    if(out2) fftw_free(out2);
    if(out3) fftw_free(out3);

    fclose(invGE);fclose(invGG);fclose(invGJ);
    fclose(invSE);fclose(invSG);fclose(invSJ);
  }

  printf("\n");
  printf("#################   Convolution  ######################\n");
  printf("\n");
}
///////////////////////////////////////////////
void convolution_GJ_to_E(double _Complex *GT,double *JT,double _Complex *ET,int itotal,int isource,int irec)
{
  int i,j;
  double scale = 1.f/itotal;
  double om;
  size_t mem_size = sizeof(fftw_complex) * itotal;
  FILE *spec_G, *spec_J, *spec_E;
  char sts1[256],sts2[256],sts3[256];
  sprintf(sts1,"./data/spectrum_confirm_j_%03d_%03d.dat",isource,irec);
  sprintf(sts2,"./data/spectrum_confirm_e_%03d_%03d.dat",isource,irec);
  sprintf(sts3,"./data/spectrum_confirm_g_%03d_%03d.dat",isource,irec);
  spec_J = fopen(sts1,"w");
  spec_E = fopen(sts2,"w");
  spec_G = fopen(sts3,"w");

  fftw_complex *in_G   = NULL;
  fftw_complex *in_J   = NULL;
  fftw_complex *in_EF  = NULL;
  fftw_complex *out_G  = NULL;
  fftw_complex *out_J  = NULL;
  fftw_complex *out_ET = NULL;
  fftw_plan p_G        = NULL;
  fftw_plan p_J        = NULL;
  fftw_plan p_E        = NULL;

  in_G   = (fftw_complex*)fftw_malloc( mem_size );
  in_J   = (fftw_complex*)fftw_malloc( mem_size );
  in_EF  = (fftw_complex*)fftw_malloc( mem_size );
  out_G  = (fftw_complex*)fftw_malloc( mem_size );
  out_J  = (fftw_complex*)fftw_malloc( mem_size );
  out_ET = (fftw_complex*)fftw_malloc( mem_size );

  if( !in_G || !out_G ){
    fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
    exit(1);
  }

  for(i=0;i<itotal;i++){
    in_G[i] = GT[i];
    in_J[i] = JT[i];
  }

  p_G = fftw_plan_dft_1d( itotal, in_G,  out_G,  FFTW_BACKWARD,  FFTW_ESTIMATE );
  p_J = fftw_plan_dft_1d( itotal, in_J,  out_J,  FFTW_BACKWARD,  FFTW_ESTIMATE );

  fftw_execute(p_G);
  fftw_execute(p_J);

  for(i=0;i<itotal;i++){
    om = 2.f*M_PI/itotal /dt;
    in_EF[i] = out_G[i]*out_J[i];
    fprintf(spec_J,"%10e %10e  %10e  %10e   %10e\n", \
              i*om/2.f/M_PI,creal(out_J[i]),creal(I*out_J[i]), creal(out_J[i])*creal(out_J[i])+creal(I*out_J[i])*creal(I*out_J[i]), atan(creal(out_J[i]*I)/creal(out_J[i])));
    fprintf(spec_G,"%10e %10e  %10e  %10e   %10e\n", \
              i*om/2.f/M_PI,creal(out_G[i]),creal(I*out_G[i]), creal(out_G[i])*creal(out_G[i])+creal(I*out_G[i])*creal(I*out_G[i]), atan(creal(out_G[i]*I)/creal(out_G[i])));
    fprintf(spec_E,"%10e %10e  %10e  %10e   %10e\n", \
              i*om/2.f/M_PI,creal(in_EF[i]),creal(I*in_EF[i]), creal(in_EF[i])*creal(in_EF[i])+creal(I*in_EF[i])*creal(I*in_EF[i]), atan(creal(in_EF[i]*I)/creal(in_EF[i])));
  }

  if(p_G ) fftw_destroy_plan(p_G);
  if(p_J ) fftw_destroy_plan(p_J);
  if(p_E ) fftw_destroy_plan(p_E);
  if(in_G ) fftw_free(in_G);
  if(in_J ) fftw_free(in_J);
  if(in_EF) fftw_free(in_EF);
  if(out_G ) fftw_free(out_G);
  if(out_J ) fftw_free(out_J);
  if(out_ET) fftw_free(out_ET);
  fclose(spec_G);fclose(spec_J);fclose(spec_E);
}
///////////////////////////////////////////////
void output_to_rec(double _Complex *ET,int itotal,int irec,int isource)
{
  int i,j;
  char str_rec[256];
  FILE *rec_file;

  sprintf(str_rec,"./data/rex_%03d_%03d.dat",isource,irec);
  rec_file = fopen(str_rec,"w");
  for(j=0;j<itotal;j++){
    fprintf(rec_file,"%e    %e   %e\n",j*dt, creal(ET[j]),  creal(I*ET[j]));
  }
  fclose(rec_file);

}
///////////////////////////////////////////////
void copytoEcal(double *EX_r,double *EY_r,double *EZ_r,int step)
{
  int i,j,k,ijk,ijks;
  int jx,jy,jz,ijkr,irec;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk,ijks)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
  for(k=0;k<iz;k++){
    for(j=0;j<iy;j++){
      for(i=0;i<ix;i++){
        ijk  = k*ix*iy + j*ix +i;
        ijks = step*ix*iy*iz + k*ix*iy + j*ix +i;

        EcalX_b[ijks] = EX_r[ijk];
        EcalY_b[ijks] = EY_r[ijk];
        EcalZ_b[ijks] = EZ_r[ijk];
      }
    }
  }
  //for(irec=0;irec<rec_num;irec++){
  //  jx=rec_px[irec];
  //  jy=rec_py[irec];
  //  jz=rec_pz[irec];
  //  ijkr  = (jz-1)*ix*iy + (jy-1)*ix +jx-1;
  //  fprintf(ofe1[irec],"%lf    %e\n",(step+0.f)*dt,EX_r[ijkr]);
  //  fprintf(ofe2[irec],"%lf    %e\n",(step+0.f)*dt,EY_r[ijkr]);
  //  fprintf(ofe3[irec],"%lf    %e\n",(step+0.f)*dt,EZ_r[ijkr]);
  //}
}
}
///////////////////////////////////////////////
void copytoEcal_b(double *EX_r,double *EY_r,double *EZ_r,int step,int idir, int irec,FILE **ofe1,FILE **ofe2,FILE **ofe3)
{
  int i,j,k,ijk,ijks,inum;
  int jx,jy,jz,iposi;

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk,ijks)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
  for(k=0;k<iz;k++){
    for(j=0;j<iy;j++){
      for(i=0;i<ix;i++){
        ijk  = k*ix*iy + j*ix +i;
        ijks = step*ix*iy*iz + k*ix*iy + j*ix +i;

        EcalX_b[ijks] = EX_r[ijk];
        EcalY_b[ijks] = EY_r[ijk];
        EcalZ_b[ijks] = EZ_r[ijk];
      }
    }
  }
}

  //jx = 13;
  //jy = 11;
  //jz = 11;
  //iposi = (jz-1)*ix*iy + (jy-1)*ix + (jx-1);
  //inum = idir*rec_num + irec;

  //fprintf(ofe1[inum],"%lf    %e\n", step*dt, EX_r[iposi]);
  //fprintf(ofe2[inum],"%lf    %e\n", step*dt, EY_r[iposi]);
  //fprintf(ofe3[inum],"%lf    %e\n", step*dt, EZ_r[iposi]);

}
///////////////////////////////////////////////
void efld_dip(int isource, int irec, int N)
{
    int  i,t;
    double moment = 1.0;
    double sig_f = 1.0f;
    double coef1, coef2, theta,time,r,theta_r, erfc, erf;
    double EXanal,EYanal,EZanal,EX_p, d_EX;
    double sx = shot_px[isource]*dx;
    double sy = shot_py[isource]*dy;
    double sz = shot_pz[isource]*dz;
    double rx = rec_px[irec]*dx;
    double ry = rec_py[irec]*dy;
    double rz = rec_pz[irec]*dz;
    double rel_rx = rx - sx;
    double rel_ry = ry - sy;
    double rel_rz = rz - sz;
    char str[256];
    FILE *of1;

    sprintf(str,"./data/anal1_%03d_%03d.dat",isource, irec);
    of1 = fopen(str,"wt");

    r = sqrt(pow((sx-rx),2)+pow((sy-ry),2)+pow((sz-rz),2));
    printf("Target offset   =  %f  [m]\n",r);

    EXanal=EX_p=0.f;
    fprintf(of1,"%e   %e   %e   %e   %e\n", 0.f, 0.f, 0.f,0.f,0.f);
    for(t=1;t<N;t++){
      time = t*dt;
      theta = sqrt(MU0*sig_f/4.f/time);
      theta_r = theta * r;
      erfc = 1.f-erff(theta_r);
      coef1 = moment/(4.f *M_PI *sig_f *pow(r,3)) * \
              ((4.f/sqrt(M_PI)*pow(theta_r,3) + 6.f/sqrt(M_PI)*theta_r)*exp(-pow(theta_r,2)) +3.f*erfc);
      coef2 = moment/(4.f *M_PI *sig_f *pow(r,3)) * \
              ((4.f/sqrt(M_PI)*pow(theta_r,3) + 2.f/sqrt(M_PI)*theta_r)*exp(-pow(theta_r,2)) +erfc);

      EXanal = coef1 * pow(rel_rx/r,2) - coef2*1.f;
      EYanal = coef1 * rel_rx*rel_ry/r/r;
      EZanal = coef1 * rel_rx*rel_rz/r/r;

      d_EX = (EXanal- EX_p)/dt;
      EX_p = EXanal;
      fprintf(of1,"%e   %e   %e   %e\n", time, d_EX, EXanal, erfc);

    }
    printf("\n");
    fclose(of1);
}
//////////////////////////////////////////
void efld_fic(int isource, int irec, int N)
{
    int t;
    double moment = 1.0;
    double sig_f = 1.0f;
    double coef1, coef2, time,r,t1, t0,beta;
    double EXanal,gamma1, d_gamma, d2_gamma;
    FILE *of2;
    double sx = shot_px[isource]*dx;
    double sy = shot_py[isource]*dy;
    double sz = shot_pz[isource]*dz;
    double rx = rec_px[irec]*dx;
    double ry = rec_py[irec]*dy;
    double rz = rec_pz[irec]*dz;
    double rel_rx = rx - sx;
    double rel_ry = ry - sy;
    double rel_rz = rz - sz;
    char str[256];

    sprintf(str,"./data/anal2_%03d_%03d.dat",isource, irec);
    of2 = fopen(str,"wt");

    beta=M_PI*pow(fmax_w,2);
    t0  =M_PI/fmax_w;

    r = sqrt(pow((sx-rx),2)+pow((sy-ry),2)+pow((sz-rz),2));

    EXanal=0.f;
    for(t=0;t<N;t++){
        time= t*dt;
        t1 = time-r/cmax;
        coef1 = (MU0/4.f/M_PI/r) * (pow(rel_rx/r,2) -1.f);
        coef2 = (MU0/4.f/M_PI/r) * (3.f*pow(rel_rx/r,2)-1.f);

        gamma1    = sqrt(beta/M_PI)*exp(-beta*pow(t1-t0,2));
        d_gamma  = -2.f*beta*(t1-t0)*sqrt(beta/M_PI)*exp(-beta*pow(t1-t0,2));
        d2_gamma = sqrt(beta/M_PI)*exp(-beta*pow(t1-t0,2))*(4.f*pow(beta,2)*pow(t1-t0,2)-2.f*beta);

        EXanal = coef1*d2_gamma + coef2*(cmax/r*d_gamma + pow(cmax/r,2)*gamma1);
        fprintf(of2,"%e    %e\n",time,EXanal);
    }
    fclose(of2);
}

///////////////////////////////////////////////
void laplaceToFreq_back(int irec,int idir)
{
  int i,j,k,n,ijk,ijks,ijkr,jx,jy,jz,w,step;
  double om,t0,beta;
  double _Complex *ctmp_j;
  double _Complex *ctmp_ex, *ctmp_ey, *ctmp_ez;
  double _Complex *ctmp_gx, *ctmp_gy, *ctmp_gz;
  double _Complex *ctmp_f1, *ctmp_f2, *ctmp_f3;

  printf("\n");
  printf("#################   Laplace transform for receiver ######################\n");

  om   = 2.f*M_PI/it /dt;
  t0   = M_PI/fmax_w;
  beta = M_PI*pow(fmax_w,2.f);
  printf("r_check1\n");

  ctmp_f1 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_f2 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_f3 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_ex = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_ey = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_ez = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_j  = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gx = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gy = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gz = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  printf("r_check2\n");

  for(jz=0;jz<iz;jz++){
    for(jy=0;jy<iy;jy++){
      for(jx=0;jx<ix;jx++){
        ijk  = jz*ix*iy + jy*ix + jx;

        for(step=0;step<it;step++){
          ijks = step*ix*iy*iz + jz*ix*iy + jy*ix + jx;
          ctmp_f1[step] = EcalX_b[ijks];
          ctmp_f2[step] = EcalY_b[ijks];
          ctmp_f3[step] = EcalZ_b[ijks];
        }


#ifdef _OPENMP
#pragma omp parallel private(k,w,ijkr)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
        // Laplace transform
        for(w=startF;w<upperF;w++){
          //w=7;
          ctmp_ex[w] = ctmp_ey[w] = ctmp_ez[w] = 0.f;
          ctmp_j[w] = 0.f;
          for(k=0;k<it;k++){
            ctmp_ex[w] += ctmp_f1[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_ey[w] += ctmp_f2[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_ez[w] += ctmp_f3[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_j[w]  += csqrt(-2.f*omega_0/I/om/w) * JX_f[k] * dt *\
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
          }
          ctmp_j[0]    = 2.f*omega_0;
          ctmp_gx[w] = ctmp_ex[w] / ctmp_j[w];
          ctmp_gy[w] = ctmp_ey[w] / ctmp_j[w];
          ctmp_gz[w] = ctmp_ez[w] / ctmp_j[w];
          ijkr = irec*3*upperF*ix*iy*iz + idir*upperF*ix*iy*iz + w*ix*iy*iz + jz*ix*iy + jy*ix + jx;

          EcalX_bb[ijkr] = ctmp_gx[w];
          EcalY_bb[ijkr] = ctmp_gy[w];
          EcalZ_bb[ijkr] = ctmp_gz[w];

        } // End of Laplace transform

}
      }
    }
  }
  printf("r_check3\n");

  free(ctmp_j);
  free(ctmp_f1);free(ctmp_ex);free(ctmp_gx);
  free(ctmp_f2);free(ctmp_ey);free(ctmp_gy);
  free(ctmp_f3);free(ctmp_ez);free(ctmp_gz);

  printf("r_check4\n");
}
///////////////////////////////////////////////
void output_backwave(int irec,int idir)
{
  int i,jx,jy,jz,w,ijk,ijkr;
  FILE *ofb1,*ofb2,*ofb3;
  char stb1[256],stb2[256],stb3[256];
  double om = 2.f*M_PI/it /dt;

  sprintf(stb1, "./data/backwave_check_x_%03d_%03d.dat", irec,idir);
  sprintf(stb2, "./data/backwave_check_y_%03d_%03d.dat", irec,idir);
  sprintf(stb3, "./data/backwave_check_z_%03d_%03d.dat", irec,idir);

  ofb1 = fopen(stb1, "w");
  ofb2 = fopen(stb2, "w");
  ofb3 = fopen(stb3, "w");

  jx= 13;
  jy= 11;
  jz= 11;

  for(w=0;w<upperF;w++){
    ijkr = irec*3*upperF*ix*iy*iz + idir*upperF*ix*iy*iz + w*ix*iy*iz + (jz-1)*ix*iy + (jy-1)*ix + jx-1;
    fprintf(ofb1, "%lf    %e    %e    %e\n", w*om/2.f/M_PI, creal(EcalX_bb[ijkr]),cimag(EcalX_bb[ijkr]), cabs(EcalX_bb[ijkr]) );
    fprintf(ofb2, "%lf    %e    %e    %e\n", w*om/2.f/M_PI, creal(EcalY_bb[ijkr]),cimag(EcalY_bb[ijkr]), cabs(EcalY_bb[ijkr]) );
    fprintf(ofb3, "%lf    %e    %e    %e\n", w*om/2.f/M_PI, creal(EcalZ_bb[ijkr]),cimag(EcalZ_bb[ijkr]), cabs(EcalZ_bb[ijkr]) );
  }

  fclose(ofb1);
  fclose(ofb2);
  fclose(ofb3);

}
///////////////////////////////////////////////
void output_fwdwave(int isource)
{
  int i,jx,jy,jz,w,ijk,ijkw,irec;
  FILE *off1,*off2,*off3;
  char stf1[256],stf2[256],stf3[256];
  double om = 2.f*M_PI/it /dt;

  printf("#################   output FWD  ######################\n");

  for(irec=0;irec<rec_num;irec++){
    sprintf(stf1, "./data/fwd_check_x_%03d_%03d.dat", irec,isource);
    sprintf(stf2, "./data/fwd_check_y_%03d_%03d.dat", irec,isource);
    sprintf(stf3, "./data/fwd_check_z_%03d_%03d.dat", irec,isource);
    off1 = fopen(stf1, "w");
    off2 = fopen(stf2, "w");
    off3 = fopen(stf3, "w");

    jx = rec_px[irec];
    jy = rec_py[irec];
    jz = rec_pz[irec];
    for(w=0;w<upperF;w++){
      ijkw = w*ix*iy*iz + (jz-1)*ix*iy + (jy-1)*ix + jx-1;
      fprintf(off1, "%lf   %e   %e   %e\n", w*om/2.f/M_PI, creal(EcalX_ff[ijkw]),cimag(EcalX_ff[ijkw]), cabs(EcalX_ff[ijkw]) );
      fprintf(off2, "%lf   %e   %e   %e\n", w*om/2.f/M_PI, creal(EcalY_ff[ijkw]),cimag(EcalY_ff[ijkw]), cabs(EcalY_ff[ijkw]) );
      fprintf(off3, "%lf   %e   %e   %e\n", w*om/2.f/M_PI, creal(EcalZ_ff[ijkw]),cimag(EcalZ_ff[ijkw]), cabs(EcalZ_ff[ijkw]) );
    }
    fclose(off1);
    fclose(off2);
    fclose(off3);
  }

}
///////////////////////////////////////////////
void laplaceToFreq(int isource)
{
  int i,j,k,n,ijk,jx,jy,jz;
  double om,t0,beta;
  double _Complex *OBSX_f,*OBSY_f,*OBSZ_f;
  size_t mem_size = sizeof(fftw_complex) * it;
  FILE *invSE, *invSG, *invSJ, *invSC;
  FILE *invSX, *invSY, *invSZ;
  char sts1[256],sts2[256],sts3[256],sts4[256],sts5[256];
  char stx[256], sty[256], stz[256];

  printf("\n");
  printf("#################   Laplace transform  ######################\n");

  om   = 2.f*M_PI/it /dt;
  t0   = M_PI/fmax_w;
  beta = M_PI*pow(fmax_w,2.f);
  OBSX_f = (double _Complex *)malloc(sizeof(double _Complex)*it);
  OBSY_f = (double _Complex *)malloc(sizeof(double _Complex)*it);
  OBSZ_f = (double _Complex *)malloc(sizeof(double _Complex)*it);

  fftw_complex *in_True  = NULL;
  fftw_complex *out_True = NULL;
  fftw_plan     p_True   = NULL;

  in_True  = (fftw_complex*)fftw_malloc( mem_size );
  out_True = (fftw_complex*)fftw_malloc( mem_size );

  for(i=0;i<it;i++) in_True[i] = signalTrue[i];

  p_True   = fftw_plan_dft_1d(it, in_True, out_True, FFTW_BACKWARD,  FFTW_ESTIMATE );
  fftw_execute(p_True);

  for(i=0;i<rec_num;i++){
    sprintf(sts1,"./data/spectrum_j_%03d_%03d.dat",isource,i);
    sprintf(sts2,"./data/spectrum_e_%03d_%03d.dat",isource,i);
    sprintf(sts3,"./data/spectrum_g_%03d_%03d.dat",isource,i);
    sprintf(sts5,"./data/spectrum_c_%03d_%03d.dat",isource,i);
    sprintf(stx, "./data/spectrum_ex_%03d_%03d.dat",isource,i);
    sprintf(sty, "./data/spectrum_ey_%03d_%03d.dat",isource,i);
    sprintf(stz, "./data/spectrum_ez_%03d_%03d.dat",isource,i);
    invSJ = fopen(sts1,"w");
    invSE = fopen(sts2,"w");
    invSG = fopen(sts3,"w");
    invSC = fopen(sts5,"w");
    invSX = fopen(stx ,"w");
    invSY = fopen(sty ,"w");
    invSZ = fopen(stz ,"w");

    for(j=0;j<it;j++){
      EX_f[j] = all_recEX[j + i*it];
      EY_f[j] = all_recEY[j + i*it];
      EZ_f[j] = all_recEZ[j + i*it];
    }

#ifdef _OPENMP
#pragma omp parallel private(n,k)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(n=0;n<it;n++){
      EX_w[n] = EY_w[n] = EZ_w[n] = JX_w[n] = 0.f;
      for(k=0;k<it;k++){
        EX_w[n] += EX_f[k] *dt*\
                   exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
        EY_w[n] += EY_f[k] *dt*\
                   exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
        EZ_w[n] += EZ_f[k] *dt*\
                   exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
        JX_w[n] += csqrt(-2.f*omega_0/I/om/n)*JX_f[k] * dt *\
                   exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
      }
      JX_w[0] = 2.f*omega_0;
      GX_w[n]  = EX_w[n] / JX_w[n];
      GY_w[n]  = EY_w[n] / JX_w[n];
      GZ_w[n]  = EZ_w[n] / JX_w[n];
      OBSX_f[n] = GX_w[n] * out_True[n];
      OBSY_f[n] = GY_w[n] * out_True[n];
      OBSZ_f[n] = GZ_w[n] * out_True[n];
    }
}

    for(n=0;n<it;n++){
      fprintf(invSE,"%10e %10e  %10e  %10e   %10e\n", \
              n*om/2.f/M_PI,creal(EX_w[n]),cimag(EX_w[n]), cabs(EX_w[n]), carg(EX_w[n]));
      fprintf(invSJ,"%10e %10e  %10e  %10e  %10e\n", \
              n*om/2.f/M_PI,creal(JX_w[n]),cimag(JX_w[n]), cabs(JX_w[n]), carg(JX_w[n]));
      fprintf(invSG,"%10e %10e  %10e  %10e  %10e\n", \
              n*om/2.f/M_PI,creal(GX_w[n]),cimag(GX_w[n]), cabs(GX_w[n]), carg(GX_w[n]));
        fprintf(invSC,"%10e %10e  %10e  %10e  %10e\n", \
              n*om/2.f/M_PI,creal(OBSX_f[n]),cimag(OBSX_f[n]), cabs(OBSX_f[n]), carg(OBSX_f[n]));
      if(n<upperF){
        fprintf(invSX,"%10e %10e  %10e  %10e  %10e\n", \
              n*om/2.f/M_PI,creal(OBSX_f[n]),cimag(OBSX_f[n]), cabs(OBSX_f[n]), carg(OBSX_f[n]));
        fprintf(invSY,"%10e %10e  %10e  %10e  %10e\n", \
              n*om/2.f/M_PI,creal(OBSY_f[n]),cimag(OBSY_f[n]), cabs(OBSY_f[n]), carg(OBSY_f[n]));
        fprintf(invSZ,"%10e %10e  %10e  %10e  %10e\n", \
              n*om/2.f/M_PI,creal(OBSZ_f[n]),cimag(OBSZ_f[n]), cabs(OBSZ_f[n]), carg(OBSZ_f[n]));
      }
    }

    fclose(invSE);fclose(invSG);fclose(invSJ);fclose(invSC);
    fclose(invSX);fclose(invSY);fclose(invSZ);
  }

  if(p_True   ) fftw_destroy_plan(p_True);
  if(in_True  ) fftw_free(in_True);
  if(out_True ) fftw_free(out_True);
  free(OBSX_f);
  free(OBSY_f);
  free(OBSZ_f);
}
///////////////////////////////////////////////
void laplaceToFreq_2(int isource)
{
  int i,j,k,n,ijk,ijkw,ijks,jx,jy,jz,w,tstep,irec,ijkr;
  double om,t0,beta;
  double _Complex *ctmp_j;
  double _Complex *ctmp_ex, *ctmp_ey, *ctmp_ez;
  double _Complex *ctmp_gx, *ctmp_gy, *ctmp_gz;
  double _Complex *ctmp_f1, *ctmp_f2, *ctmp_f3;
  size_t mem_size = sizeof(fftw_complex) * it;

  printf("\n");
  printf("#################   Laplace transform for transmitter ######################\n");

  om   = 2.f*M_PI/it /dt;
  t0   = M_PI/fmax_w;
  beta = M_PI*pow(fmax_w,2.f);

  fftw_complex *in_True  = NULL;
  fftw_complex *out_True = NULL;
  fftw_plan     p_True   = NULL;

  printf("check 1\n");
  in_True  = (fftw_complex*)fftw_malloc( mem_size );
  out_True = (fftw_complex*)fftw_malloc( mem_size );

  for(i=0;i<it;i++) in_True[i] = signalTrue[i];

  p_True   = fftw_plan_dft_1d(it, in_True, out_True, FFTW_BACKWARD,  FFTW_ESTIMATE );
  fftw_execute(p_True);
  printf("check 2\n");

  ctmp_f1 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_f2 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_f3 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_ex = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_ey = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_ez = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_j  = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gx = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gy = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gz = (double _Complex *)malloc(sizeof(double _Complex)*upperF);

  printf("check 3\n");
  for(jz=0;jz<iz;jz++){
    for(jy=0;jy<iy;jy++){
      for(jx=0;jx<ix;jx++){
        ijk  = jz*ix*iy + jy*ix + jx;

        for(tstep=0;tstep<it;tstep++){
          ijks = tstep*ix*iy*iz + jz*ix*iy + jy*ix + jx;
          ctmp_f1[tstep] = EcalX_b[ijks];
          ctmp_f2[tstep] = EcalY_b[ijks];
          ctmp_f3[tstep] = EcalZ_b[ijks];
        }

#ifdef _OPENMP
#pragma omp parallel private(w,k,ijkw,ijks)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
        // Laplace transform
        for(w=startF;w<upperF;w++){
          //w=7;
          ctmp_ex[w] = ctmp_ey[w] = ctmp_ez[w] = 0.f;
          ctmp_j[w] = 0.f;
          for(k=0;k<it;k++){
            ctmp_ex[w] += ctmp_f1[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_ey[w] += ctmp_f2[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_ez[w] += ctmp_f3[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_j[w] += csqrt(-2.f*omega_0/I/om/w) * JX_f[k] * dt *\
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
          }
          ctmp_j[0]    = 2.f*omega_0;
          ctmp_gx[w] = ctmp_ex[w] / ctmp_j[w];
          ctmp_gy[w] = ctmp_ey[w] / ctmp_j[w];
          ctmp_gz[w] = ctmp_ez[w] / ctmp_j[w];

          ijkw = isource*upperF*ix*iy*iz + w*ix*iy*iz + jz*ix*iy + jy*ix + jx;
          ijks = w*ix*iy*iz + jz*ix*iy + jy*ix + jx;
          // convolution G * J
          EcalX_ff[ijks] = ctmp_gx[w] * out_True[w];
          EcalY_ff[ijks] = ctmp_gy[w] * out_True[w];
          EcalZ_ff[ijks] = ctmp_gz[w] * out_True[w];
          EcalX_ffL[ijkw] = ctmp_gx[w] * out_True[w];
          EcalY_ffL[ijkw] = ctmp_gy[w] * out_True[w];
          EcalZ_ffL[ijkw] = ctmp_gz[w] * out_True[w];
        }
}

      }
    }
  }
  printf("check 4\n");

  free(ctmp_j);
  free(ctmp_f1);free(ctmp_ex);free(ctmp_gx);
  free(ctmp_f2);free(ctmp_ey);free(ctmp_gy);
  free(ctmp_f3);free(ctmp_ez);free(ctmp_gz);
  printf("check 5\n");
  if(p_True   ) fftw_destroy_plan(p_True);
  if(in_True  ) fftw_free(in_True);
  if(out_True ) fftw_free(out_True);
  printf("check 6\n");

}
///////////////////////////////////////////////
void residualE(int iter)
{
    int i,irec,jx,jy,jz,ijk,ijks,step,isr,w;
    FILE *ofr1,*ofr2,*ofr3;
    char res1[256],res2[256],res3[256];

    sprintf(res1,"./data/residualEX_%03d.dat",iter);
    sprintf(res2,"./data/residualEY_%03d.dat",iter);
    sprintf(res3,"./data/residualEZ_%03d.dat",iter);

    ofr1 = fopen(res1,"w");
    ofr2 = fopen(res2,"w");
    ofr3 = fopen(res3,"w");

    printf("\n");
    printf("#################   Calculation residual  ######################\n");

    for(irec=0;irec<rec_num;irec++){
      jx  = rec_px[irec];
      jy  = rec_py[irec];
      jz  = rec_pz[irec];
      for(w=0; w<upperF;w++){
        ijks = w*ix*iy*iz +  ix*iy*(jz-1) + ix*(jy-1) + (jx-1);
        isr  = irec*upperF + w;

        delEx_f[isr] = Exobs[isr] - EcalX_ff[ijks];
        delEy_f[isr] = Eyobs[isr] - EcalY_ff[ijks];
        delEz_f[isr] = Ezobs[isr] - EcalZ_ff[ijks];

        fprintf(ofr1,"%d   %e   %e   %e   %e   %e   %e\n", w, creal(delEx_f[isr]),cimag(delEx_f[isr]),creal(EcalX_ff[ijks]),cimag(EcalX_ff[ijks]),creal(Exobs[isr]),cimag(Exobs[isr]));
        fprintf(ofr2,"%d   %e   %e   %e   %e   %e   %e\n", w, creal(delEy_f[isr]),cimag(delEy_f[isr]),creal(EcalY_ff[ijks]),cimag(EcalY_ff[ijks]),creal(Exobs[isr]),cimag(Exobs[isr]));
        fprintf(ofr3,"%d   %e   %e   %e   %e   %e   %e\n", w, creal(delEz_f[isr]),cimag(delEz_f[isr]),creal(EcalZ_ff[ijks]),cimag(EcalZ_ff[ijks]),creal(Exobs[isr]),cimag(Exobs[isr]));

        err_sum[iter] += cabs(delEx_f[isr]);
        err_sum[iter] += cabs(delEy_f[isr]);
        err_sum[iter] += cabs(delEz_f[isr]);
      }

    }
    fclose(ofr1);
    fclose(ofr2);
    fclose(ofr3);
}
///////////////////////////////////////////////
void sensitivity()
{
    int i,j,k,jx,jy,jz,ijk,ijks,ijkr,irec,idir,w;

    printf("\n");
    printf("#################   Calculation Sensitivity  ######################\n");

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk,irec,idir,w,ijks,ijkr)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk = k*ix*iy + j*ix + i;
          sense[ijk] = 0.f;
          for(irec=0;irec<rec_num;irec++){
            for(idir=0;idir<3;idir++){
              for(w=0;w<upperF;w++){
                //w=7;
                ijks = w*ix*iy*iz + k*ix*iy + j*ix + i;
                ijkr = irec*3*upperF*ix*iy*iz + idir*upperF*ix*iy*iz \
                     + w*ix*iy*iz + k*ix*iy + j*ix + i;
                //sense[ijk] += pow(cabs(EcalX_ff[ijks]),2.f);
                //sense[ijk] += pow(cabs(EcalY_ff[ijks]),2.f);
                //sense[ijk] += pow(cabs(EcalZ_ff[ijks]),2.f);
                //sense[ijk] += cabs(EcalX_ff[ijks]);
                //sense[ijk] += cabs(EcalY_ff[ijks]);
                //sense[ijk] += cabs(EcalZ_ff[ijks]);
                sense[ijk] += cabs(EcalX_ff[ijks] * EcalX_bb[ijkr]);
                sense[ijk] += cabs(EcalY_ff[ijks] * EcalY_bb[ijkr]);
                sense[ijk] += cabs(EcalZ_ff[ijks] * EcalZ_bb[ijkr]);
              }
            }
          }

        }
      }
    }
}

}
///////////////////////////////////////////////
void conv_GdelE()
{
    int i,j,k,jx,jy,jz,ijk,ijkr,irec,idir,isr,ijkw,w;

    printf("\n");
    printf("#################   Calculation Convolution  ######################\n");

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          for(w=0;w<upperF;w++){
            ijkw = w*ix*iy*iz + k*ix*iy + j*ix + i;
            EcalX_back[ijkw] = EcalY_back[ijkw] = EcalZ_back[ijkw] = 0.f;
          }
        }
      }
    }

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,irec,idir,w,ijkw,ijkr,isr)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          for(irec=0;irec<rec_num;irec++){
            for(idir=0;idir<3;idir++){
              for(w=0;w<upperF;w++){
                ijkw = w*ix*iy*iz + k*ix*iy + j*ix + i;
                ijkr = irec*3*upperF*ix*iy*iz + idir*upperF*ix*iy*iz \
                     + w*ix*iy*iz + k*ix*iy + j*ix + i;
                isr  = irec*upperF + w;

                if(idir==0){
                  EcalX_back[ijkw] += EcalX_bb[ijkr] * conj(delEx_f[isr]);
                  EcalY_back[ijkw] += EcalY_bb[ijkr] * conj(delEx_f[isr]);
                  EcalZ_back[ijkw] += EcalZ_bb[ijkr] * conj(delEx_f[isr]);
                }else if(idir==1){
                  EcalX_back[ijkw] += EcalX_bb[ijkr] * conj(delEy_f[isr]);
                  EcalY_back[ijkw] += EcalY_bb[ijkr] * conj(delEy_f[isr]);
                  EcalZ_back[ijkw] += EcalZ_bb[ijkr] * conj(delEy_f[isr]);
                }else if(idir==2){
                  EcalX_back[ijkw] += EcalX_bb[ijkr] * conj(delEz_f[isr]);
                  EcalY_back[ijkw] += EcalY_bb[ijkr] * conj(delEz_f[isr]);
                  EcalZ_back[ijkw] += EcalZ_bb[ijkr] * conj(delEz_f[isr]);
                }
              }
            }
          }
        }
      }
    }
}

}
///////////////////////////////////////////////
void calc_gamma()
{
    int i,j,k,jx,jy,jz,ijk,ijkr,irec,ijkw,w;
    double om, coef;
    FILE *of_grad;
    of_grad = fopen("./data/gradient.dat","w");

    printf("\n");
    printf("#################   Calculation Gamma  ######################\n");

#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk,w,ijkw)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=ncpml;k<iz-ncpml;k++){
      for(j=ncpml;j<iy-ncpml;j++){
        for(i=ncpml;i<ix-ncpml;i++){
          ijk = k*ix*iy + j*ix + i;
          for(w=0;w<upperF;w++){
            ijkw = w*ix*iy*iz + k*ix*iy + j*ix + i;
            grad[ijk] += creal( EcalX_ff[ijkw] * EcalX_back[ijkw]  \
                              + EcalY_ff[ijkw] * EcalY_back[ijkw]  \
                              + EcalZ_ff[ijkw] * EcalZ_back[ijkw] ) /sense[ijk];
            fprintf(of_grad,"%lf   %lf   %lf   %e\n", i*dx, j*dy,k*dz, grad[ijk]);
          }
        }
      }
    }
}
  fclose(of_grad);
}
///////////////////////////////////////////////
void init_iterate()
{
    int i,j,k,jx,jy,jz,ijk,ijkw,irec,step,w;

    k1 = k2 = 0;
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ijk,w,ijkw)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk       = k*ix*iy + j*ix + i;
          grad[ijk] = 0.f;
          for(w=0;w<upperF;w++){
            ijkw = w*ix*iy*iz + k*ix*iy + j*ix + i;
            EcalX_back[ijkw] = 0.f;
            EcalY_back[ijkw] = 0.f;
            EcalZ_back[ijkw] = 0.f;
          }
        }
      }
    }
}
}

///////////////////////////////////////////////
void banner1(int idir,int irec)
{

  if(idir==0)  printf("This is %2d loop for X direction RECEIVER \n",irec);
  if(idir==1)  printf("This is %2d loop for Y direction RECEIVER \n",irec);
  if(idir==2)  printf("This is %2d loop for Z direction RECEIVER \n",irec);
  printf("[PROGRESS  =>]  ");

}

///////////////////////////////////////////////
void update_para(int iter)
{
    int i,j,k,ijk;
    char str_gra[256],str_sig[256];
    double alpha2, beta, maxv, b1, b2;
    double bb[ix*iy*iz];
    FILE *Fsig, *Fgra;

    sprintf(str_sig, "./data/sig_%03d.dat",iter);
    sprintf(str_gra, "./data/gradient_%03d.dat",iter);
    Fsig = fopen(str_sig,"w");
    Fgra = fopen(str_gra, "w");

    alpha2 =  sqrt(err_sum[iter]/err_sum[0]);

    maxv=0.f;
    b1 = b2 = 0.f;
    for(k=ncpml;k<iz-ncpml;k++){
      for(j=ncpml;j<iy-ncpml;j++){
        for(i=ncpml;i<ix-ncpml;i++){
          ijk = k*ix *iy + j*ix +i;
          if(maxv <= fabs(grad[ijk]))  maxv = fabs(grad[ijk]);
        }
      }
    }

    if(iter==0){
      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            if(k<=iseabed) sig[ijk] += alpha2 * grad[ijk];
            //if(k<=iseabed) sig[ijk] += alpha2 * grad[ijk]/ maxv;
            if(sig[ijk]<minval) sig[ijk] = minval;
            if(sig[ijk]>maxval) sig[ijk] = maxval;
          }
        }
      }
    }else{
      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            b1 += grad[ijk]   *  (grad[ijk] - grad_p[ijk]);
            b2 += grad_p[ijk] * grad_p[ijk];
          }
        }
      }

      beta = b1 / b2;
      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            bb[ijk]   = -grad[ijk] + beta * bb[ijk];
            if(maxv <= fabs(bb[ijk]))  maxv = fabs(bb[ijk]);
          }
        }
      }

      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            if(k<=iseabed) sig[ijk] += alpha2 * bb[ijk];
            //if(k<=iseabed) sig[ijk] += alpha2 * bb[ijk] / maxv;
            if(sig[ijk]<minval) sig[ijk] = minval;
            if(sig[ijk]>maxval) sig[ijk] = maxval;
          }
        }
      }

    }

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
            ijk = k*ix *iy + j*ix +i;
            fprintf(Fgra, "%10.2lf    %10.2lf    %10.2lf    %e\n", i*dx, j*dy, k*dz, grad[ijk]);
            fprintf(Fsig, "%10.2lf    %10.2lf    %10.2lf    %e\n", i*dx, j*dy, k*dz, sig[ijk]);
            grad_p[ijk] = grad[ijk];
            grad[ijk] = 0.f;
        }
      }
    }

    fclose(Fgra); fclose(Fsig);
}
//////////////////////////////////////////
void backpropagation(int idir, int irec, int step){

  e_field_cpml42(EX,EY,EZ,HX,HY,HZ);
  epml_abcs4(EX,EY,EZ,HX,HY,HZ);
  if(idir==0) read_backwave_3d(EX, signalX,step,irec);
  if(idir==1) read_backwave_3d(EY, signalX,step,irec);
  if(idir==2) read_backwave_3d(EZ, signalX,step,irec);
  h_field_cpml42(EX,EY,EZ,HX,HY,HZ);
  hpml_abcs4(EX,EY,EZ,HX,HY,HZ);

}
//////////////////////////////////////////
void fwdpropagation(int isource, int step){

  e_field_cpml42(EX,EY,EZ,HX,HY,HZ);
  epml_abcs4(EX,EY,EZ,HX,HY,HZ);
  read_source_3d(EZ,signalX,isource,step);
  h_field_cpml42(EX,EY,EZ,HX,HY,HZ);
  hpml_abcs4(EX,EY,EZ,HX,HY,HZ);

}

//////////////////////////////////////////
void ThreeCfwdpropagation(int isource, int step, int idir){

  e_field_cpml42(EX,EY,EZ,HX,HY,HZ);
  epml_abcs4(EX,EY,EZ,HX,HY,HZ);
  if(idir ==0) read_source_3d(EX,signalX,isource,step);
  if(idir ==1) read_source_3d(EY,signalX,isource,step);
  if(idir ==2) read_source_3d(EZ,signalX,isource,step);
  h_field_cpml42(EX,EY,EZ,HX,HY,HZ);
  hpml_abcs4(EX,EY,EZ,HX,HY,HZ);

}

//////////////////////////////////////////
void init(){

  checking_omp();
  read_shotrec();
  init_eh_field_3d();
  init_model();
  read_trawave();
  read_waveform();
  lattice_time_3d();
  init_cpml();
  init_gradient();
}
//////////////////////////////////////////
void calc_phi(){
  int i,j,k,ijk;
  double sigmax2,gradmax;
  FILE *ophi;
  sigmax2 = sig[0];
  gradmax = grad[0];
  for(k=0;k<iz;k++){
    for(j=0;j<iy;j++){
      for(i=0;i<ix;i++){
        ijk = k*ix*iy + j*ix + i;
        if(k<=iseabed){
          if(sig[ijk]  > sigmax2)  sigmax2 = sig[ijk];
          if(grad[ijk] > gradmax)  gradmax = grad[ijk];
        }
      }
    }
  }

  scaler = 0.01f * sigmax2/gradmax;
  //scaler = 0.01f * sigmax2;
  ophi = fopen("./data/sig_tmp.dat","w");
  printf("\n");
  printf("#################   Calculation Phi  ######################\n");
  printf("sigmax :   %lf\n",sigmax2);
  printf("gradmax:   %e\n",gradmax);
  printf("scaler :   %e\n",scaler);

  for(k=0;k<iz;k++){
    for(j=0;j<iy;j++){
      for(i=0;i<ix;i++){
        ijk = k*ix*iy + j*ix + i;
        sig2[ijk] = sig[ijk] - scaler * grad[ijk];
        fprintf(ophi, "%lf   %lf   %lf   %lf\n", i*dx, j*dy, k*dz, sig2[ijk]);
      }
    }
  }
  fclose(ophi);

}
//////////////////////////////////////////
void media_coeff_sig_tmp()
{
    int i,j,k,ijk;
    double eps2;
    FILE *ofm3,*ofm4,*ofm5;

    ofm3 = fopen("./data/coef1.dat","w");
    ofm4 = fopen("./data/coef2.dat","w");
    ofm5 = fopen("./data/coef3.dat","w");

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk  = k*iy*ix + j*ix + i;
          eps2 = sig2[ijk] /2.f /omega_0;
          // CPML Coefficient
          // coef1
          ca_x[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
          ca_y[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));
          ca_z[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
          // coef2
          da_x[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*eps2))) \
                      / (1.0f + ((msigx[i]*dt)/(2.f*eps2)));
          da_y[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*eps2))) \
                      / (1.0f + ((msigy[j]*dt)/(2.f*eps2)));
          da_z[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*eps2))) \
                      / (1.0f + ((msigz[k]*dt)/(2.f*eps2)));
          // coef3
          cb_x[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2));
          cb_y[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2));
          cb_z[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2));
          // coef4
          db_x[ijk] = dt/MU0 /(1.f+(msigx[i]*dt)/(2.f*eps2));
          db_y[ijk] = dt/MU0 /(1.f+(msigy[j]*dt)/(2.f*eps2));
          db_z[ijk] = dt/MU0 /(1.f+(msigz[k]*dt)/(2.f*eps2));

          fprintf(ofm3,"%+8e  %+8e  %+8e  %+8e  %+8e  %+8e\n",ca_x[ijk],ca_y[ijk],ca_z[ijk],cb_x[ijk],cb_y[ijk],cb_z[ijk]);
          fprintf(ofm4,"%+8e  %+8e  %+8e  %+8e  %+8e  %+8e\n",da_x[ijk],da_y[ijk],da_z[ijk],db_x[ijk],db_y[ijk],db_z[ijk]);
          fprintf(ofm5,"%8e  %8e  %8e  %8e  %8e  %8e\n",esigx[i],esigy[j],esigz[k],msigx[i],msigy[j],msigz[k]);
        }
      }
    }
    printf("eps2  :    %e\n",eps2);
    fclose(ofm3);
    fclose(ofm4);
    fclose(ofm5);

}
///////////////////////////////////////////////
void laplaceToFreq_3()
{
  int i,j,k,n,ijk,ijkw,ijks,jx,jy,jz,w,tstep;
  double om,t0,beta;
  double _Complex *ctmp_j;
  double _Complex *ctmp_ex, *ctmp_ey, *ctmp_ez;
  double _Complex *ctmp_gx, *ctmp_gy, *ctmp_gz;
  double _Complex *ctmp_f1, *ctmp_f2, *ctmp_f3;
  size_t mem_size = sizeof(fftw_complex) * it;

  printf("\n");
  printf("#################   Laplace transform for transmitter ######################\n");

  om   = 2.f*M_PI/it /dt;
  t0   = M_PI/fmax_w;
  beta = M_PI*pow(fmax_w,2.f);

  fftw_complex *in_True  = NULL;
  fftw_complex *out_True = NULL;
  fftw_plan     p_True   = NULL;

  in_True  = (fftw_complex*)fftw_malloc( mem_size );
  out_True = (fftw_complex*)fftw_malloc( mem_size );

  for(i=0;i<it;i++) in_True[i] = signalTrue[i];

  p_True   = fftw_plan_dft_1d(it, in_True, out_True, FFTW_BACKWARD,  FFTW_ESTIMATE );
  fftw_execute(p_True);

  ctmp_f1 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_f2 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_f3 = (double _Complex *)malloc(sizeof(double _Complex)*it);
  ctmp_ex = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_ey = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_ez = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_j  = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gx = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gy = (double _Complex *)malloc(sizeof(double _Complex)*upperF);
  ctmp_gz = (double _Complex *)malloc(sizeof(double _Complex)*upperF);

  for(jz=0;jz<iz;jz++){
    for(jy=0;jy<iy;jy++){
      for(jx=0;jx<ix;jx++){
        ijk  = jz*ix*iy + jy*ix + jx;

        for(tstep=0;tstep<it;tstep++){
          ijks = tstep*ix*iy*iz + jz*ix*iy + jy*ix + jx;
          ctmp_f1[tstep] = EcalX_b[ijks];
          ctmp_f2[tstep] = EcalY_b[ijks];
          ctmp_f3[tstep] = EcalZ_b[ijks];
        }

#ifdef _OPENMP
#pragma omp parallel private(w,k,ijkw)
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
        // Laplace transform
        for(w=startF;w<upperF;w++){
          //w=7;
          ctmp_ex[w] = ctmp_ey[w] = ctmp_ez[w] = 0.f;
          ctmp_j[w] = 0.f;
          for(k=0;k<it;k++){
            ctmp_ex[w] += ctmp_f1[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_ey[w] += ctmp_f2[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_ez[w] += ctmp_f3[k] * dt * \
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
            ctmp_j[w] += csqrt(-2.f*omega_0/I/om/w) * JX_f[k] * dt *\
                   exp(-sqrt(omega_0*om*w)*k*dt)*cexp(I*sqrt(omega_0*om*w)*k*dt);
          }
          ctmp_j[0]    = 2.f*omega_0;
          ctmp_gx[w] = ctmp_ex[w] / ctmp_j[w];
          ctmp_gy[w] = ctmp_ey[w] / ctmp_j[w];
          ctmp_gz[w] = ctmp_ez[w] / ctmp_j[w];
          ijkw = w*ix*iy*iz + jz*ix*iy + jy*ix + jx;
          // convolution G * J
          EcalX_ff2[ijkw] = ctmp_gx[w] * out_True[w];
          EcalY_ff2[ijkw] = ctmp_gy[w] * out_True[w];
          EcalZ_ff2[ijkw] = ctmp_gz[w] * out_True[w];
        }
}

      }
    }
  }

  free(ctmp_j);
  free(ctmp_f1);free(ctmp_ex);free(ctmp_gx);
  free(ctmp_f2);free(ctmp_ey);free(ctmp_gy);
  free(ctmp_f3);free(ctmp_ez);free(ctmp_gz);
  if(p_True   ) fftw_destroy_plan(p_True);
  if(in_True  ) fftw_free(in_True);
  if(out_True ) fftw_free(out_True);

}
///////////////////////////////////////////////
void calc_alpha(int isource){
  int jx,jy,jz,w,ijkw,irec,ijks,isr;
  double _Complex ctmp_x,ctmp_y,ctmp_z;

  printf("\n");
  printf("#################   Calculation alpha  ######################\n");

  for(irec=0;irec<rec_num;irec++){
    jx  = rec_px[irec];
    jy  = rec_py[irec];
    jz  = rec_pz[irec];
    for(w=0; w<upperF;w++){
      ijks = isource*upperF*ix*iy*iz + w*ix*iy*iz + (jz-1)*ix*iy + (jy-1)*ix + jx-1;
      ijkw = w*ix*iy*iz + (jz-1)*ix*iy + (jy-1)*ix + jx-1;
      isr  = irec*upperF + w;

      ctmp_x = EcalX_ff2[ijkw] - EcalX_ffL[ijks];
      ctmp_y = EcalY_ff2[ijkw] - EcalY_ffL[ijks];
      ctmp_z = EcalZ_ff2[ijkw] - EcalZ_ffL[ijks];

      k1 += delEx_f[isr] * ctmp_x;
      k1 += delEy_f[isr] * ctmp_y;
      k1 += delEz_f[isr] * ctmp_z;

      k2 += ctmp_x * ctmp_x;
      k2 += ctmp_y * ctmp_y;
      k2 += ctmp_z * ctmp_z;

    }
  }

}
///////////////////////////////////////////////
void update_para2(int iter)
{
    int i,j,k,ijk;
    char str_gra[256],str_sig[256];
    double alpha2, beta, maxv, b1, b2;
    double bb[ix*iy*iz];
    FILE *Fsig, *Fgra;

    sprintf(str_sig, "./data/sig_%03d.dat",iter);
    sprintf(str_gra, "./data/gradient_%03d.dat",iter);
    Fsig = fopen(str_sig,"w");
    Fgra = fopen(str_gra, "w");

    alpha = 3.f * cabs(k1)/cabs(k2);
    printf("k1   : %e   %e\n",creal(k1),cimag(k1));
    printf("k2   : %e   %e\n",creal(k2),cimag(k2));
    printf("alpha : %e\n",alpha);

    for(k=ncpml;k<iz-ncpml;k++){
      for(j=ncpml;j<iy-ncpml;j++){
        for(i=ncpml;i<ix-ncpml;i++){
          ijk = k*ix *iy + j*ix +i;
          if(k<=iseabed) sig[ijk] += alpha * scaler * grad[ijk];
          if(sig[ijk]<minval) sig[ijk] = minval;
          if(sig[ijk]>maxval) sig[ijk] = maxval;
        }
      }
    }

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
            ijk = k*ix *iy + j*ix +i;
            fprintf(Fgra, "%10.2lf    %10.2lf    %10.2lf    %e\n", i*dx, j*dy, k*dz, alpha*scaler*grad[ijk]);
            fprintf(Fsig, "%10.2lf    %10.2lf    %10.2lf    %e\n", i*dx, j*dy, k*dz, sig[ijk]);
            grad_p[ijk] = grad[ijk];
            grad[ijk] = 0.f;
        }
      }
    }

    fclose(Fgra); fclose(Fsig);
}
//////////////////////////////////////////
void update_para3(int iter)
{
    int i,j,k,ijk;
    char str_gra[256],str_sig[256];
    double alpha2, beta, maxv, b1, b2;
    double bb[ix*iy*iz];
    FILE *Fsig, *Fgra;

    sprintf(str_sig, "./data/sig_%03d.dat",iter);
    sprintf(str_gra, "./data/gradient_%03d.dat",iter);
    Fsig = fopen(str_sig,"w");
    Fgra = fopen(str_gra, "w");

    alpha = 3.f * cabs(k1)/cabs(k2);
    printf("k1   : %e   %e\n",creal(k1),cimag(k1));
    printf("k2   : %e   %e\n",creal(k2),cimag(k2));
    printf("alpha : %e\n",alpha);

    b1 = b2 = 0.f;
    if(iter==0){
      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            if(k<=iseabed) sig[ijk] += alpha * scaler * grad[ijk];
            if(sig[ijk]<minval) sig[ijk] = minval;
            if(sig[ijk]>maxval) sig[ijk] = maxval;
          }
        }
      }
    }else{
      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            b1 += grad[ijk] * (grad[ijk] - grad_p[ijk]);
            b2 += grad_p[ijk] * grad_p[ijk];
          }
        }
      }

      beta = b1 /b2;

      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            bb[ijk] = -grad[ijk] + beta * bb[ijk];
          }
        }
      }
      for(k=ncpml;k<iz-ncpml;k++){
        for(j=ncpml;j<iy-ncpml;j++){
          for(i=ncpml;i<ix-ncpml;i++){
            ijk = k*ix *iy + j*ix +i;
            if(k<=iseabed) sig[ijk] += alpha * scaler * bb[ijk];
            if(sig[ijk]<minval) sig[ijk] = minval;
            if(sig[ijk]>maxval) sig[ijk] = maxval;
          }
        }
      }
    }

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
            ijk = k*ix *iy + j*ix +i;
            fprintf(Fgra, "%10.2lf    %10.2lf    %10.2lf    %e\n", i*dx, j*dy, k*dz, alpha*scaler*grad[ijk]);
            fprintf(Fsig, "%10.2lf    %10.2lf    %10.2lf    %e\n", i*dx, j*dy, k*dz, sig[ijk]);
            grad_p[ijk] = grad[ijk];
            grad[ijk] = 0.f;
        }
      }
    }

    fclose(Fgra); fclose(Fsig);
}
//////////////////////////////////////////
void show_error(int iter,FILE *ofer){

  printf("error sum %03d  is   %e\n",iter, err_sum[iter]/err_sum[0] );
  fprintf(ofer,"%03d  %lf\n",iter, err_sum[iter]/err_sum[0] );

}
