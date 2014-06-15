#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fwi3d_cpu_function_v1.1.h"
#include "fwi3d_cpu_alloex1_v1.1.h"
#include "fwi3d_cpu_alloex2_v1.1.h"

///////////////////////////////////////////////
void read_shotrec()
{
  int i,itmp;
  FILE *ifp1;
  ifp1 = fopen("./data/model_env.dat","r");

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

  fclose(ifp1);

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
    }
    fclose(ift1);
    fclose(ift2);
    fclose(ift3);

}
///////////////////////////////////////////////
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
    eps = (double *)malloc(sizeof(double)*100);
    sig = (double *)malloc(sizeof(double)*100);
    mu  = (double *)malloc(sizeof(double)*100);

    // e_coeff1
    cex = (double *)malloc(sizeof(double)*n);
    cey = (double *)malloc(sizeof(double)*n);
    cez = (double *)malloc(sizeof(double)*n);
    // e_coeff2
    cexry = (double *)malloc(sizeof(double)*n);
    cexrz = (double *)malloc(sizeof(double)*n);
    ceyrx = (double *)malloc(sizeof(double)*n);
    ceyrz = (double *)malloc(sizeof(double)*n);
    cezrx = (double *)malloc(sizeof(double)*n);
    cezry = (double *)malloc(sizeof(double)*n);
    // h_coeff
    chxry = (double *)malloc(sizeof(double)*n);
    chxrz = (double *)malloc(sizeof(double)*n);
    chyrx = (double *)malloc(sizeof(double)*n);
    chyrz = (double *)malloc(sizeof(double)*n);
    chzrx = (double *)malloc(sizeof(double)*n);
    chzry = (double *)malloc(sizeof(double)*n);
    // mur
//    ex1_xy = (double *)malloc(sizeof(double)*4*ix*iy);
//    ey1_xy = (double *)malloc(sizeof(double)*4*ix*iy);
//    ey1_yz = (double *)malloc(sizeof(double)*4*iy*iz);
//    ez1_yz = (double *)malloc(sizeof(double)*4*iy*iz);
//    ez1_zx = (double *)malloc(sizeof(double)*4*ix*iz);
//    ex1_zx = (double *)malloc(sizeof(double)*4*ix*iz);
//    ex2_xy = (double *)malloc(sizeof(double)*4*ix*iy);
//    ey2_xy = (double *)malloc(sizeof(double)*4*ix*iy);
//    ey2_yz = (double *)malloc(sizeof(double)*4*iy*iz);
//    ez2_yz = (double *)malloc(sizeof(double)*4*iy*iz);
//    ez2_zx = (double *)malloc(sizeof(double)*4*ix*iz);
//    ex2_zx = (double *)malloc(sizeof(double)*4*ix*iz);
//
//    exz1 = (double *)malloc(sizeof(double)*4*ix*iy);
//    eyz1 = (double *)malloc(sizeof(double)*4*ix*iy);
//    eyx1 = (double *)malloc(sizeof(double)*4*iy*iz);
//    ezx1 = (double *)malloc(sizeof(double)*4*iy*iz);
//    exy1 = (double *)malloc(sizeof(double)*4*ix*iz);
//    ezy1 = (double *)malloc(sizeof(double)*4*ix*iz);
//    exz2 = (double *)malloc(sizeof(double)*4*ix*iy);
//    eyz2 = (double *)malloc(sizeof(double)*4*ix*iy);
//    eyx2 = (double *)malloc(sizeof(double)*4*iy*iz);
//    ezx2 = (double *)malloc(sizeof(double)*4*iy*iz);
//    exy2 = (double *)malloc(sizeof(double)*4*ix*iz);
//    ezy2 = (double *)malloc(sizeof(double)*4*ix*iz);
    // pml
    esigx = (double *)malloc(sizeof(double)*ix);
    esigy = (double *)malloc(sizeof(double)*iy);
    esigz = (double *)malloc(sizeof(double)*iz);
    msigx = (double *)malloc(sizeof(double)*ix);
    msigy = (double *)malloc(sizeof(double)*iy);
    msigz = (double *)malloc(sizeof(double)*iz);
    exy = (double *)malloc(sizeof(double)*ix*iy*iz);
    exz = (double *)malloc(sizeof(double)*ix*iy*iz);
    eyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    eyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    ezx = (double *)malloc(sizeof(double)*ix*iy*iz);
    ezy = (double *)malloc(sizeof(double)*ix*iy*iz);
    hxy = (double *)malloc(sizeof(double)*ix*iy*iz);
    hxz = (double *)malloc(sizeof(double)*ix*iy*iz);
    hyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    hyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    hzx = (double *)malloc(sizeof(double)*ix*iy*iz);
    hzy = (double *)malloc(sizeof(double)*ix*iy*iz);
    cexy = (double *)malloc(sizeof(double)*ix*iy*iz);
    cexz = (double *)malloc(sizeof(double)*ix*iy*iz);
    ceyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    ceyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    cezx = (double *)malloc(sizeof(double)*ix*iy*iz);
    cezy = (double *)malloc(sizeof(double)*ix*iy*iz);
    chxy = (double *)malloc(sizeof(double)*ix*iy*iz);
    chxz = (double *)malloc(sizeof(double)*ix*iy*iz);
    chyx = (double *)malloc(sizeof(double)*ix*iy*iz);
    chyz = (double *)malloc(sizeof(double)*ix*iy*iz);
    chzx = (double *)malloc(sizeof(double)*ix*iy*iz);
    chzy = (double *)malloc(sizeof(double)*ix*iy*iz);
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
}
///////////////////////////////////////////////
void media_coeff_3d()
{
    int i,j,k,ijk;
    double eps2;
    FILE *ofm1,*ofm2,*ofm3,*ofm4,*ofm5;
    nmedia = 1;

// id = 0 test
    eps[0] = EPS0*1.0f;
    sig[0] = 1.0f;
     mu[0] = MU0;
// id = 1 seawater
    eps[1] = EPS0*1.0f;
    sig[1] = 0.02f;
     mu[1] = MU0;
// id = 2 basalt
    eps[2] = EPS0*1.0f;
    sig[2] = 0.45f;
     mu[2] = MU0;
// id = 3 submarine massive sulphide
    eps[3] = EPS0*100.0f;
    sig[3] = 10.0f;
     mu[3] = MU0;

// seawater and basalt's boundary
    int iseabed = 18;

    ofm1 = fopen("./data/conductivity.dat","w");
    ofm2 = fopen("./data/epsilon.dat","w");

    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk = k*iy*ix + j*ix + i;
          id[ijk]=0;
//          if(k>iseabed) id[ijk] = 0;
          fprintf(ofm1,"%3d  %3d  %3d   %8lf\n",i,j,k,sig[id[ijk]]);
          fprintf(ofm2,"%3d  %3d  %3d   %8e\n",i,j,k,eps[id[ijk]]);
        }
      }
    }
    fclose(ofm1);
    fclose(ofm2);

    ofm3 = fopen("./data/coef1.dat","w");
    ofm4 = fopen("./data/coef2.dat","w");
    ofm5 = fopen("./data/coef3.dat","w");
    for(k=0;k<iz;k++){
      for(j=0;j<iy;j++){
        for(i=0;i<ix;i++){
          ijk  = k*iy*ix + j*ix + i;
          eps2 = sig[id[ijk]] /2.f /omega_0;
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
          cexy[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));
          cexz[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
          ceyx[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
          ceyz[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
          cezx[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
          cezy[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));

          chxy[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigy[j]*dt)/(2.f*mu[id[ijk]])));
          chxz[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigz[k]*dt)/(2.f*mu[id[ijk]])));
          chyx[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigx[i]*dt)/(2.f*mu[id[ijk]])));
          chyz[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigz[k]*dt)/(2.f*mu[id[ijk]])));
          chzx[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigx[i]*dt)/(2.f*mu[id[ijk]])));
          chzy[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigy[j]*dt)/(2.f*mu[id[ijk]])));

          cex[ijk]   = 1.f;
          cey[ijk]   = 1.f;
          cez[ijk]   = 1.f;
          // dt/eps/(1+(sig*dt)/(2*eps)) / dl
          cexry[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2)) /dy;
          cexrz[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2)) /dz;
          ceyrx[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2)) /dx;
          ceyrz[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2)) /dz;
          cezrx[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2)) /dx;
          cezry[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2)) /dy;
          // dt/mu /dl
          chxry[ijk] = dt/mu[id[ijk]] /(1.f+(msigy[j]*dt)/(2.f*mu[id[ijk]])) /dy;
          chxrz[ijk] = dt/mu[id[ijk]] /(1.f+(msigz[k]*dt)/(2.f*mu[id[ijk]])) /dz;
          chyrx[ijk] = dt/mu[id[ijk]] /(1.f+(msigx[i]*dt)/(2.f*mu[id[ijk]])) /dx;
          chyrz[ijk] = dt/mu[id[ijk]] /(1.f+(msigz[k]*dt)/(2.f*mu[id[ijk]])) /dz;
          chzrx[ijk] = dt/mu[id[ijk]] /(1.f+(msigx[i]*dt)/(2.f*mu[id[ijk]])) /dx;
          chzry[ijk] = dt/mu[id[ijk]] /(1.f+(msigy[j]*dt)/(2.f*mu[id[ijk]])) /dy;
          // CPML Coefficient
          // coef1
          ca_x[ijk]   = (1.0f - ((esigx[i]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigx[i]*dt)/(2.f*eps2)));
          ca_y[ijk]   = (1.0f - ((esigy[j]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigy[j]*dt)/(2.f*eps2)));
          ca_z[ijk]   = (1.0f - ((esigz[k]*dt)/(2.f*eps2))) \
                      / (1.0f + ((esigz[k]*dt)/(2.f*eps2)));
          // coef2
          da_x[ijk]   = (1.0f - ((msigx[i]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigx[i]*dt)/(2.f*mu[id[ijk]])));
          da_y[ijk]   = (1.0f - ((msigy[j]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigy[j]*dt)/(2.f*mu[id[ijk]])));
          da_z[ijk]   = (1.0f - ((msigz[k]*dt)/(2.f*mu[id[ijk]]))) \
                      / (1.0f + ((msigz[k]*dt)/(2.f*mu[id[ijk]])));
          // coef3
          cb_x[ijk] = dt/eps2 /(1.f+(esigx[i]*dt)/(2.f*eps2));
          cb_y[ijk] = dt/eps2 /(1.f+(esigy[j]*dt)/(2.f*eps2));
          cb_z[ijk] = dt/eps2 /(1.f+(esigz[k]*dt)/(2.f*eps2));
          // coef4
          db_x[ijk] = dt/mu[id[ijk]] /(1.f+(msigx[i]*dt)/(2.f*mu[id[ijk]]));
          db_y[ijk] = dt/mu[id[ijk]] /(1.f+(msigy[j]*dt)/(2.f*mu[id[ijk]]));
          db_z[ijk] = dt/mu[id[ijk]] /(1.f+(msigz[k]*dt)/(2.f*mu[id[ijk]]));

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
  char sthx[256],sthy[256],sthz[256];

  ofe1 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofe2 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofe3 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofh1 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofh2 = (FILE **)malloc(rec_num*sizeof(FILE *));
  ofh3 = (FILE **)malloc(rec_num*sizeof(FILE *));

  for(i=0;i<rec_num;i++){
    sprintf(strx,"./data/rex_%03d_%03d.dat",isource,i);
    sprintf(stry,"./data/rey_%03d_%03d.dat",isource,i);
    sprintf(strz,"./data/rez_%03d_%03d.dat",isource,i);
    sprintf(sthx,"./data/rhx_%03d_%03d.dat",isource,i);
    sprintf(sthy,"./data/rhy_%03d_%03d.dat",isource,i);
    sprintf(sthz,"./data/rhz_%03d_%03d.dat",isource,i);
    ofe1[i] = fopen(strx,"w");
    ofe2[i] = fopen(stry,"w");
    ofe3[i] = fopen(strz,"w");
    ofh1[i] = fopen(sthx,"w");
    ofh2[i] = fopen(sthy,"w");
    ofh3[i] = fopen(sthz,"w");
  }
}

///////////////////////////////////////////////
void close_FILE(void)
{
    fclose(*ofe1); fclose(*ofe2); fclose(*ofe3);
    fclose(*ofh1); fclose(*ofh2); fclose(*ofh3);
}
///////////////////////////////////////////////
void set_zero()
{
    int i;
    n=ix*iy*iz;

    for(i=0;i<n;i++){
        exy[i]=exz[i]=0.f;
        eyx[i]=eyz[i]=0.f;
        ezx[i]=ezy[i]=0.f;
        hxy[i]=hxz[i]=0.f;
        hyx[i]=hyz[i]=0.f;
        hzx[i]=hzy[i]=0.f;
        EX[i]=EY[i]=EZ[i]=0.f;
        HX[i]=HY[i]=HZ[i]=0.f;
        id[i]=0;
    }
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
                     dt*2.f*omega_0/sig[id[sourcep]]/dx/dy/dz;

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
        fprintf(ofe1[i],"%lf    %e\n",(step+1.f)*dt,EX[ijk]);
        fprintf(ofe2[i],"%lf    %e\n",(step+1.f)*dt,EY[ijk]);
        fprintf(ofe3[i],"%lf    %e\n",(step+1.f)*dt,EZ[ijk]);
        fprintf(ofh1[i],"%lf    %e\n",(step+1.f)*dt,HX[ijk]);
        fprintf(ofh2[i],"%lf    %e\n",(step+1.f)*dt,HY[ijk]);
        fprintf(ofh3[i],"%lf    %e\n",(step+1.f)*dt,HZ[ijk]);
      }
// output in files at all point 
      if((step%10) == 0){
        FILE *fptime1,*fptime2,*fptime3;
        sprintf(stt1, "./data2/time_ex_%05d.dat",step);
        sprintf(stt2, "./data2/time_ey_%05d.dat",step);
        sprintf(stt3, "./data2/time_ez_%05d.dat",step);
        fptime1 = fopen(stt1,"w");
        fptime2 = fopen(stt2,"w");
        fptime3 = fopen(stt3,"w");

        for(k=0;k<iz;k++){
          for(j=0;j<iy;j++){
            for(i=0;i<ix;i++){
              ijk=ix*iy*k + j*ix + i;
              fprintf(fptime1,"%12lf   %12lf   %12lf   %15.5e\n", i*dx, j*dy,k*dz, EX[ijk]);
              //fprintf(fptime2,"%12lf   %12lf   %12lf   %15.5e\n", i*dx, j*dy,k*dz, EY[ijk]);
              //fprintf(fptime3,"%12lf   %12lf   %12lf   %15.5e\n", i*dx, j*dy,k*dz, EZ[ijk]);
            }
          }
        }
        fclose(fptime1);
        fclose(fptime2);
        fclose(fptime3);
      }
}
///////////////////////////////////////////////
void checking_omp(void)
{
    #ifdef _OPENMP
        printf("\n");
        printf("#################   Using OpenMP  ######################\n");
    #else
        printf("\n");
        printf("#################   NO OpenMP  ######################\n");
    #endif
}

///////////////////////////////////////////////
void init_gradient(void)
{
    int i;
    fy2 = (double *)malloc( it * sizeof(double));
    fy3 = (double *)malloc( it * sizeof(double));
    fy4 = (double *)malloc( it * sizeof(double));
    fz2 = (double *)malloc( it * sizeof(double));
    fz3 = (double *)malloc( it * sizeof(double));
    fz4 = (double *)malloc( it * sizeof(double));

    grad   = (double *)malloc(ix*iy*iz*sizeof(double));
    grad_p = (double *)malloc(ix*iy*iz*sizeof(double));
    delEx  = (double *)malloc(rec_num*it*sizeof(double));
    delEy  = (double *)malloc(rec_num*it*sizeof(double));
    delEz  = (double *)malloc(rec_num*it*sizeof(double));

    for(i=0;i<ix*iy*iz;i++){
        grad[i] = grad_p[i] = 0.f;
    }
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
void read_Eobs()
{
  int i,j;
  double dtmp;
  FILE *ifp_Ex,*ifp_Ey,*ifp_Ez;
  char ch1[256],ch2[256],ch3[256];

  Exobs = (double *)malloc(it*rec_num*sizeof(double));
  Eyobs = (double *)malloc(it*rec_num*sizeof(double));
  Ezobs = (double *)malloc(it*rec_num*sizeof(double));

  for(i=0; i<rec_num; i++){
      /* storage of observed waveform */
      sprintf(ch1, "./data/rex_%03d_%03d.dat",isource,i);
      sprintf(ch2, "./data/rey_%03d_%03d.dat",isource,i);
      sprintf(ch3, "./data/rez_%03d_%03d.dat",isource,i);
      ifp_Ex = fopen(ch1,"r");
      ifp_Ey = fopen(ch2,"r");
      ifp_Ez = fopen(ch3,"r");
      for(j=0;j<it;j++){
          fscanf(ifp_Ex,"%lf     %lf\n",&dtmp, &Exobs[i*it+j]);
          fscanf(ifp_Ey,"%lf     %lf\n",&dtmp, &Eyobs[i*it+j]);
          fscanf(ifp_Ez,"%lf     %lf\n",&dtmp, &Ezobs[i*it+j]);
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
    }
}
///////////////////////////////////////////////
void copytoarray(double *EX,double *EY,double *EZ,int step,FILE **ofe1,FILE **ofe2,FILE **ofe3)
{
    int i,jx,jy,jz,inum,ijk;

    for(i=0;i<rec_num;i++){
        jx   = rec_px[i];
        jy   = rec_py[i];
        jz   = rec_pz[i];
        inum = it*i*step;
        ijk  = ix*iy*jz + ix*jy + jx;

        rec_ex[inum] = EX[ijk];
        rec_ey[inum] = EY[ijk];
        rec_ez[inum] = EZ[ijk];

        fprintf(ofe1[i], "%e     %e\n",step*dt, EX[ijk]);
        fprintf(ofe2[i], "%e     %e\n",step*dt, EY[ijk]);
        fprintf(ofe3[i], "%e     %e\n",step*dt, EZ[ijk]);
    }
}
///////////////////////////////////////////////
void calc_delE(int isource,int iter)
{
    int i,j,ij;
    FILE *ofs1,*ofs2,*ofs3;
    char ch4[256],ch5[256],ch6[256];

    for(i=0;i<rec_num;i++){
      sprintf(ch4, "./data/delEx_%03d_%03d_%03d.dat",isource,i,iter);
      sprintf(ch5, "./data/delEy_%03d_%03d_%03d.dat",isource,i,iter);
      sprintf(ch6, "./data/delEz_%03d_%03d_%03d.dat",isource,i,iter);
      ofs1 = fopen(ch4,"w");
      ofs2 = fopen(ch5,"w");
      ofs3 = fopen(ch6,"w");
      for(j=0;j<it;j++){
        ij = i*it + j;

        delEx[ij] = Exobs[ij] - rec_ex[ij];
        delEy[ij] = Eyobs[ij] - rec_ey[ij];
        delEz[ij] = Ezobs[ij] - rec_ez[ij];

        fprintf(ofs1,"%e      %e\n",j*dt, delEx[ij]);
        fprintf(ofs2,"%e      %e\n",j*dt, delEy[ij]);
        fprintf(ofs3,"%e      %e\n",j*dt, delEz[ij]);
      }
      fclose(ofs1);
      fclose(ofs2);
      fclose(ofs3);
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
void read_backwave_3d(double *EX,double *delEx,int step){
    int i,jx,jy,jz,ijk,ijkt;

    for(i=0;i<rec_num;i++){
        jx      = rec_px[i];
        jy      = rec_py[i];
        jz      = rec_pz[i];

        ijk     = ix*iy*jz + ix*jy + jx;
        ijkt    = i*it+(it-1-step);
        EX[ijk] = delEx[ijkt];
    }
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
// ex
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
///////////////////////////////////////////////
void h_field_cpml42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
    //double c1 = 1.14443f, c2=-0.04886f;
// hx
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
///////////////////////////////////////////////
void e_field_cpml422(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
// ex
    for(k=ncpml+1;k<iz-ncpml-0;k++){
      for(j=ncpml+1;j<iy-ncpml-0;j++){
        for(i=ncpml+0;i<ix-ncpml-0;i++){
          ijk = k*ix*iy + j*ix + i;
          EX[ijk] = ca_x[ijk] * EX[ijk] \
                  + cb_x[ijk] * ((c1*HZ[ijk] - c1*HZ[ijk-ix   ] + c2*HZ[ijk+ix   ] -c2*HZ[ijk-2*ix   ]) / kedy[j] \
                              -  (c1*HY[ijk] - c1*HY[ijk-ix*iy] + c2*HY[ijk+ix*iy] -c2*HY[ijk-2*ix*iy]) / kedz[k]);
        }
      }
    }
// ey
    for(k=ncpml+1;k<iz-ncpml-0;k++){
      for(j=ncpml+0;j<iy-ncpml-0;j++){
        for(i=ncpml+1;i<ix-ncpml-0;i++){
          ijk = k*ix*iy + j*ix + i;
          EY[ijk] = ca_y[ijk] *  EY[ijk] \
                  + cb_y[ijk] * ((c1*HX[ijk] - c1*HX[ijk-ix*iy] + c2*HX[ijk+ix*iy] - c2*HX[ijk-2*ix*iy]) / kedz[k] \
                              -  (c1*HZ[ijk] - c1*HZ[ijk-1    ] + c2*HZ[ijk+1    ] - c2*HZ[ijk-2      ]) / kedx[i]);
        }
      }
    }
// ez
    for(k=ncpml+0;k<iz-ncpml-0;k++){
      for(j=ncpml+1;j<iy-ncpml-0;j++){
        for(i=ncpml+1;i<ix-ncpml-0;i++){
          ijk = k*ix*iy + j*ix + i;
          EZ[ijk] = ca_z[ijk] *  EZ[ijk] \
                  + cb_z[ijk] * ((c1*HY[ijk] - c1*HY[ijk-1 ] + c2*HY[ijk+1 ] - c2*HY[ijk-2   ]) / kedx[i] \
                              -  (c1*HX[ijk] - c1*HX[ijk-ix] + c2*HX[ijk+ix] - c2*HX[ijk-2*ix]) / kedy[j]);
        }
      }
    }
// ex2
    for(k=1;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          if( i<ncpml+0 || i>=(ix-ncpml-0) \
           || j<ncpml+1 || j>=(iy-ncpml-0) \
           || k<ncpml+1 || k>=(iz-ncpml-0) ){
             ijk = k*ix*iy + j*ix + i;
             EX[ijk] = ca_x[ijk] *   EX[ijk] \
                     + cb_x[ijk] * ((HZ[ijk] - HZ[ijk - ix   ]) / kedy[j] \
                                  - (HY[ijk] - HY[ijk - ix*iy]) / kedz[k]);
          }
        }
      }
    }
// ey2
    for(k=1;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          if( i<ncpml+1 || i>=(ix-ncpml-0) \
           || j<ncpml+0 || j>=(iy-ncpml-0) \
           || k<ncpml+1 || k>=(iz-ncpml-0) ){
             ijk = k*ix*iy + j*ix + i;
             EY[ijk] = ca_y[ijk] *   EY[ijk] \
                     + cb_y[ijk] * ((HX[ijk] - HX[ijk - ix*iy]) / kedz[k] \
                                  - (HZ[ijk] - HZ[ijk - 1    ]) / kedx[i]);
          }
        }
      }
    }
// ez2
    for(k=0;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          if( i<ncpml+1 || i>=(ix-ncpml-0) \
           || j<ncpml+1 || j>=(iy-ncpml-0) \
           || k<ncpml+0 || k>=(iz-ncpml-0) ){
             ijk = k*ix*iy + j*ix + i;
             EZ[ijk] = ca_z[ijk] *   EZ[ijk] \
                     + cb_z[ijk] * ((HY[ijk] - HY[ijk - 1 ]) / kedx[i] \
                                  - (HX[ijk] - HX[ijk - ix]) / kedy[j]);
          }
        }
      }
    }
}
///////////////////////////////////////////////
void h_field_cpml422(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double c1 = 1.125f, c2=-0.04167f;
// hx
    for(k=ncpml+0;k<iz-ncpml-1;k++){
      for(j=ncpml+0;j<iy-ncpml-1;j++){
        for(i=ncpml+0;i<ix-ncpml-0;i++){
          ijk = k*ix*iy + j*ix + i;
          HX[ijk] = da_x[ijk] *   HX[ijk] \
                  - db_x[ijk] * ((c1*EZ[ijk+ix   ] - c1*EZ[ijk] + c2*EZ[ijk+2*ix   ] - c2*EZ[ijk-ix   ]) / khdy[j] \
                              -  (c1*EY[ijk+ix*iy] - c1*EY[ijk] + c2*EY[ijk+2*ix*iy] - c2*EY[ijk-ix*iy]) / khdz[k]);
        }
      }
    }
// hy
    for(k=ncpml+0;k<iz-ncpml-1;k++){
      for(j=ncpml+0;j<iy-ncpml-0;j++){
        for(i=ncpml+0;i<ix-ncpml-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HY[ijk] = da_y[ijk] *   HY[ijk] \
                  - db_y[ijk] * ((c1*EX[ijk+ix*iy] - c1*EX[ijk] + c2*EX[ijk+2*ix*iy] - c2*EX[ijk-ix*iy]) / khdz[k] \
                              -  (c1*EZ[ijk+1    ] - c1*EZ[ijk] + c2*EZ[ijk+2      ] - c2*EZ[ijk-1    ]) / khdx[i]);
        }
      }
    }
// hz
    for(k=ncpml+0;k<iz-ncpml-0;k++){
      for(j=ncpml+0;j<iy-ncpml-1;j++){
        for(i=ncpml+0;i<ix-ncpml-1;i++){
          ijk = k*ix*iy + j*ix + i;
          HZ[ijk] = da_z[ijk] *   HZ[ijk] \
                  - db_z[ijk] * ((c1*EY[ijk+1 ] - c1*EY[ijk] + c2*EY[ijk+2   ] - c2*EY[ijk-1 ]) / khdx[i] \
                              -  (c1*EX[ijk+ix] - c1*EX[ijk] + c2*EX[ijk+2*ix] - c2*EX[ijk-ix]) / khdy[j]);
        }
      }
    }
// hx2
    for(k=0;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
          if( i<ncpml+0 || i>=(ix-ncpml-0) \
           || j<ncpml+0 || j>=(iy-ncpml-1) \
           || k<ncpml+0 || k>=(iz-ncpml-1) ){
             ijk = k*ix*iy + j*ix + i;
             HX[ijk] = da_x[ijk] *   HX[ijk] \
                     - db_x[ijk] * ((EZ[ijk + ix   ] - EZ[ijk]) / khdy[j] \
                                  - (EY[ijk + ix*iy] - EY[ijk]) / khdz[k]);
          }
        }
      }
    }
// hy2
    for(k=0;k<iz-1;k++){
      for(j=1;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          if( i<ncpml+0  || i>=(ix-ncpml-1) \
           || j<ncpml+0  || j>=(iy-ncpml-0) \
           || k<ncpml+0  || k>=(iz-ncpml-1) ){
             ijk = k*ix*iy + j*ix + i;
             HY[ijk] = da_y[ijk] *   HY[ijk] \
                     - db_y[ijk] * ((EX[ijk + ix*iy] - EX[ijk]) / khdz[k] \
                                  - (EZ[ijk + 1    ] - EZ[ijk]) / khdx[i]);
          }
        }
      }
    }
// hz2
    for(k=1;k<iz-1;k++){
      for(j=0;j<iy-1;j++){
        for(i=0;i<ix-1;i++){
          if( i<ncpml+0 || i>=(ix-ncpml-1) \
           || j<ncpml+0 || j>=(iy-ncpml-1) \
           || k<ncpml+0 || k>=(iz-ncpml-0) ){
             ijk = k*ix*iy + j*ix + i;
             HZ[ijk] = da_z[ijk] *   HZ[ijk] \
                     - db_z[ijk] * ((EY[ijk + 1 ] - EY[ijk]) / khdx[i] \
                                  - (EX[ijk + ix] - EX[ijk]) / khdy[j]);
          }
        }
      }
    }
}
///////////////////////////////////////////////

