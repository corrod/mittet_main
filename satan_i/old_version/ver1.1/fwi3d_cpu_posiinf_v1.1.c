#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define CC 2.99790000e8
#define MU0 12.5663706e-7
#define nx 121
#define ny 121
#define nz 121
int time_check(double, double, double,int);
int fmax_check(double, double, double, double,double,double);
int exp_check(double, double, double, double,double,int,double);
int maxin3(double, double, double);

int main()
{
  int i;
  int rec_num,shot_num;
  int *shotx,*shoty,*shotz;
  int *recx,*recy,*recz;
  int tt,mstep,trstep;  // 何ステップごとにデータを記録するか
  int ncpml;
  double omega_0,f0,fmax_w,Glim;
  double dx,dy,dz;
  double dt_ratio;
  FILE *fp1,*fp2;

  fp1 = fopen("./data/model_env.dat","w");
  fp2 = fopen("./data/posi4gmt.dat","w");

  shot_num = 1;
  rec_num  = 1;

  dx = dy = dz = 5.f;
  tt       = 600;
  mstep    = 4;
  dt_ratio = 100.;
  f0       = 1.e-1;   // 
  omega_0  = 2.f*M_PI*f0;
  fmax_w   = 1.e+1; //送信波の最大周波数
  Glim     = 10.4f;
  ncpml    = 30;
  trstep   = tt/mstep;

  shotx = (int *)malloc(sizeof(int) * shot_num);
  shoty = (int *)malloc(sizeof(int) * shot_num);
  shotz = (int *)malloc(sizeof(int) * shot_num);
  recx  = (int *)malloc(sizeof(int) * rec_num);
  recy  = (int *)malloc(sizeof(int) * rec_num);
  recz  = (int *)malloc(sizeof(int) * rec_num);

  for(i=0;i<shot_num;i++){
    shotx[i] = 51+0*i;
    shoty[i] = 61;
    shotz[i] = 61;
  }
  for(i=0;i<rec_num;i++){
    recx[i]  = 61+0*i;
    recy[i]  = 61;
    recz[i]  = 61;
  }

// Output in model_env.dat 
  fprintf(fp1,"%d\n",nx);
  fprintf(fp1,"%d\n",ny);
  fprintf(fp1,"%d\n",nz);
  fprintf(fp1,"%f\n",dx);
  fprintf(fp1,"%f\n",dy);
  fprintf(fp1,"%f\n",dz);
  fprintf(fp1,"%f\n",dt_ratio);
  fprintf(fp1,"%d\n",tt);
  fprintf(fp1,"%d\n",mstep);
  fprintf(fp1,"%d\n",trstep);
  fprintf(fp1,"%e\n",omega_0);
  fprintf(fp1,"%e\n",fmax_w);
  fprintf(fp1,"%e\n",Glim);
  fprintf(fp1,"%d\n",ncpml);
  fprintf(fp1,"%d\n",shot_num);
  for(i=0;i<shot_num;i++){
    fprintf(fp1,"%4d  %4d  %4d\n", shotx[i],shoty[i],shotz[i]);
    fprintf(fp2,"%4d  %4d  %4d\n", shotx[i],shoty[i],shotz[i]);
  }

  fprintf(fp1,"%d\n",rec_num);
  for(i=0;i<rec_num;i++){
    fprintf(fp1,"%4d  %4d  %4d\n", recx[i],recy[i],recz[i]);
    fprintf(fp2,"%4d  %4d  %4d\n", recx[i],recy[i],recz[i]);
  }

  time_check(dx,dy,dz,tt);
  fmax_check(dx,dy,dz,omega_0,fmax_w,Glim);
//  exp_check(dx,dy,dz,omega_0,fmax_w,tt,Glim);

  fclose(fp1); fclose(fp2);
  free(shotx); free(shoty); free(shotz);
  free(recx);  free(recy);  free(recz);
  return 0;
}

////////////////////////////////////
int time_check(double dx,double dy,double dz, int N)
{
    double S=0.6f;
    int Nt;

    Nt=(int)(S*nx/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2)));

    printf("Nt     =  %d [step]\n",Nt);
    if(Nt >N) printf("******* time step tt is violated *******\n");

    return 0;
}
////////////////////////////////////
int fmax_check(double dx,double dy,double dz,double omega_0,double fmax_w,double Glim)
{
    double fmax, cmin;

    cmin = sqrt(2.f*omega_0/MU0/1.0f);
    fmax = cmin /Glim /maxin3(dx,dy,dz);

    printf("fmax   =  %f \n",fmax);
    if(fmax <fmax_w) printf("******* fmax is violated *******\n");
    return 0;
}
////////////////////////////////////
int exp_check(double dx,double dy,double dz,double omega_0,double fmax_w,int N,double Glim)
{
    int i;
    double f_max,fmax, cmin,om_max,dt,vc,est,om_lim;
    double f0_lim = 0.1f;
    double fmax_w_lim = 1.0f;
    double omega_lim = 2.f*M_PI*f0_lim;
    double vc_p,dt_p,om_max_p,est_p;
    double dx_p=100.f,dy_p=100.f,dz_p=100.f;
    FILE *icheck;

    icheck = fopen("./data/check.dat","wt");

    vc_p = sqrt(2.f * omega_lim / MU0 /1.f);
    dt_p = 1.0f/vc_p/sqrt(1.0/pow(dx_p,2) + 1.0/pow(dy_p,2) + 1.0/pow(dz_p,2));
    om_max_p  = 2.f*M_PI/dt_p;   // omega の最大値

    vc = sqrt(2.f * omega_0 / MU0 /1.f);
    dt = 1.0f/vc/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    om_max  = 2.f*M_PI/dt;   // omega の最大値
    for(i=0;i<N;i++){
      est_p = -sqrt(om_max_p*omega_lim)*i*dt_p;   // expの右上
      est = -sqrt(om_max*omega_0)*i*dt;   // expの右上
      fprintf(icheck,"%e   %e\n",exp(est_p),exp(est));
    }

    printf("basic dt       = %lf \n", dt_p);
    printf("dt             = %lf \n", dt);
    printf("exponent       = %lf \n",est);
    printf("Limit exponent = %lf \n",est_p);
    if(est_p > est)   printf("********  expornt is too small *******\n");

    return 0;
}
////////////////////////////////////
int maxin3(double dx,double dy,double dz)
{
    double tmp;
    if(dx > dy){
        tmp = dx;
    }else{
        tmp = dy;
    }

    if(tmp > dz){
        tmp = tmp;
    }else{
        tmp = dz;
    }
    return tmp;
}
