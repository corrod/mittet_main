#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define LMAX 10000
#define MU0 12.5663706e-7
#define EXP 13
void read_para();
double dt;

int main()
{
  int i,ii;
  double dtmp1, dtmp2,dtmp3, dtmp4, dtmp5;
  double kyori[99];
  FILE *if1;
  char str1[256];

  read_para();
  kyori[ 0] = 100;
  kyori[ 1] = 200;
  kyori[ 2] = 400;
  kyori[ 3] = 800;
  kyori[ 4] = 1000;
  kyori[ 5] = 1200;
  kyori[ 6] = 1400;
  kyori[ 7] = 1600;
  kyori[ 8] = 1800;
  kyori[ 9] = 2000;
  kyori[10] = 3200;
  kyori[11] = 5000;
  kyori[12] = 6400;

  for(i=0;i<13;i++){
    sprintf(str1, "./data_hokan/spectrum_2_%04d", (int)(kyori[i]) );

    if1 = fopen(str1,"r");
    for(ii=0;ii<LMAX;ii++){
      fscanf(if1, "%lf    %le    %lf    %lf    %lf\n",&dtmp1,&dtmp2, &dtmp3, &dtmp4, &dtmp5);
      if(ii==EXP) printf("%04d    %e    %e\n", (int)(kyori[i]), dtmp4, dtmp5);
    }
    fclose(if1);
  }

}
///////////////////////////////////
void read_para(void)
{
  int tmp,i;
  double Glim;
  double dx,dy,dz, dtmp,omega_0,fmax_w;
  double cmax,cmin,courant;
  FILE *ifp1;
  ifp1 = fopen("./data/model_env.dat","rt");

  fscanf(ifp1,"%d", &tmp);
  fscanf(ifp1,"%d", &tmp);
  fscanf(ifp1,"%d", &tmp);
  fscanf(ifp1,"%lf", &dx);
  fscanf(ifp1,"%lf", &dy);
  fscanf(ifp1,"%lf", &dz);
  fscanf(ifp1,"%lf", &dtmp);
  fscanf(ifp1,"%d",  &tmp);
  fscanf(ifp1,"%d",  &tmp);
  fscanf(ifp1,"%d",  &tmp);
  fscanf(ifp1,"%lf", &omega_0);
  fscanf(ifp1,"%lf", &fmax_w);
  fscanf(ifp1,"%lf", &Glim);
  fscanf(ifp1,"%d", &tmp);
  fclose(ifp1);


  //courant = 1.0f/velocity/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  //dt = courant *dt_ratio;
  //courant_ori = 1.0f/velocity/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  cmin = sqrt(2.f*omega_0/MU0/1.f);
  cmax = sqrt(2.f*omega_0/MU0/1.f);
  courant = 1.0f/cmax/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  dt = courant*6.f/7.f*0.999f;
}

