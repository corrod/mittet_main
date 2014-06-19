#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define CC 2.99790000e8
#define MU0 12.5663706e-7
#define nx 31  // X方向グリッド数
#define ny 31  // Y方向グリッド数
#define nz 31  // Z方向グリッド数
int time_check(double, double, double,int);
int fmax_check(double, double, double, double,double,double);
int maxin3(double, double, double);

int main()
{
  int i;
  int rec_num,shot_num;
  int *shotx,*shoty,*shotz;
  int *recx,*recy,*recz;
  int tt,mstep,trstep;  // 何ステップごとにデータを記録するか
  int ncpml,startF,upperF;
  double omega_0,f0,fmax_w,Glim;
  double dx,dy,dz;
  double dt_ratio;
  double FL,FH,FS,AP,AS;
  FILE *fp1,*fp2;

  fp1 = fopen("./data/model_env.dat","w");
  fp2 = fopen("./data/posi4gmt.dat","w");

  shot_num = 7;
  rec_num  = 7;

  dx = dy = dz = 30.f; // grid幅
  tt       = 3000;     // total のステップ数
  mstep    = 4;        // 不要
  dt_ratio = 100.; // 不要
  f0       = 1.e+0;   // f0が小さいとdtがでかくなる
  omega_0  = 2.f*M_PI*f0;
  fmax_w   = 3.e+0; //送信波の最大周波数
  Glim     = 10.4f;
  startF   = 7;   // inversion用 (ちょっとフォワードでも使う
  upperF   = 8;   // inversion用 (ちょっとフォワードでも使う
  ncpml    = 5;   // Convolutional PMLのgrid数
  FL       = 4.f; // BandPass 用　不要
  FH       = 8.f; // BandPass 用　不要
  FS       = 12.f; // BandPass 用　不要
  AP       = 0.5f; // BandPass 用　不要
  AS       = 5.0f; // BandPass 用　不要

  trstep   = tt/mstep; // iran

  shotx = (int *)malloc(sizeof(int) * shot_num);
  shoty = (int *)malloc(sizeof(int) * shot_num);
  shotz = (int *)malloc(sizeof(int) * shot_num);
  recx  = (int *)malloc(sizeof(int) * rec_num);
  recy  = (int *)malloc(sizeof(int) * rec_num);
  recz  = (int *)malloc(sizeof(int) * rec_num);

  for(i=0;i<shot_num;i++){
    shotx[i] = 10+2*i;
    shoty[i] = 16;
    shotz[i] = 17;
  }
  for(i=0;i<rec_num;i++){
    recx[i]  = 10+2*i;
    recy[i]  = 16;
    recz[i]  = 16;
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
  fprintf(fp1,"%f\n",FL);
  fprintf(fp1,"%f\n",FH);
  fprintf(fp1,"%f\n",FS);
  fprintf(fp1,"%f\n",AP);
  fprintf(fp1,"%f\n",AS);
  fprintf(fp1,"%d\n",startF);
  fprintf(fp1,"%d\n",upperF);
  fprintf(fp1,"%d\n",shot_num);
  for(i=0;i<shot_num;i++){
    fprintf(fp1,"%4d  %4d  %4d\n", shotx[i],shoty[i],shotz[i]);
    fprintf(fp2,"%10.2lf  %10.2lf  %10.2lf\n", shotx[i]*dx,shoty[i]*dy,shotz[i]*dz);
  }

  fprintf(fp1,"%d\n",rec_num);
  for(i=0;i<rec_num;i++){
    fprintf(fp1,"%4d  %4d  %4d\n", recx[i],recy[i],recz[i]);
    fprintf(fp2,"%10.2lf  %10.2lf  %10.2lf\n", recx[i]*dx,recy[i]*dy,recz[i]*dz);
  }

  //time_check(dx,dy,dz,tt);
  //fmax_check(dx,dy,dz,omega_0,fmax_w,Glim);

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
