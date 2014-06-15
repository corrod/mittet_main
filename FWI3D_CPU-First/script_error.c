#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define LMAX 40000
#define MU0 12.5663706e-7
void read_para();
double dt;

int main()
{
    int i,ii;
    double dtmp1, dtmp2,dtmp3, dtmp4, dtmp5, err_sum, maxval;
    double ctmp1, ctmp2,ctmp3;
    double dtime[LMAX], real1[LMAX],real2[LMAX], kyori[10];
    FILE *if1,*if2;
    char str1[256],str2[256];

    read_para();
    kyori[0] = 250;
    kyori[1] = 500;
    kyori[2] = 1000;
    kyori[3] = 2000;
    kyori[4] = 3000;
    kyori[5] = 4000;
    kyori[6] = 5000;

    for(i=0;i<5;i++){
      err_sum = 0.f;
      maxval  = 0.f;
      sprintf(str1, "./data_hokan/anal_%04d", (int)(kyori[i]) );
      sprintf(str2, "./data_hokan/inverse_g_%04d", (int)(kyori[i]) );

      if1 = fopen(str1,"r");
      if2 = fopen(str2,"r");
      for(ii=0;ii<LMAX;ii++){
          fscanf(if1, "%lf    %le    %lf    %lf    %lf\n",&dtmp1,&dtmp2, &dtmp3, &dtmp4, &dtmp5);
          fscanf(if2, "%lf    %le    %lf\n",&ctmp1,&ctmp2, &ctmp3);

          if(maxval <= dtmp2) maxval = dtmp2;
          if(ii>2){
            if(ii*dt < 100)  err_sum += fabs(ctmp2 - dtmp2)\
                                     *(log10(ii*dt)-log10((ii-1)*dt));
          }
      }
      printf("err_sum %04d   =   %e\n", (int)(kyori[i]), err_sum /maxval);
      fclose(if1);fclose(if2);
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
  printf(" dt  = %lf\n", dt);
}

