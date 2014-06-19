#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#define MU0 12.5663706e-7
void init_var(void);
void init_emfield(void);
void efld_dip(void);
void efld_fic(void);
int maxin3(int,int,int);

int N,nx,ny,nz;
double dx, dy, dz, dt, omega_0, fmax_w,vc,cmax,cmin;
double FL,FH,FS,AP,AS;
double *EX_f;
double _Complex invJX;
double _Complex *EX_w, *EX_t;
double _Complex *JX_f, *JX_w, *JX_t;
double _Complex *GX_w,*GX_t;

int main()
{
  int ii,i,j,tmp,k,n,nmax;
  double dtmp, courant,beta,t0,om_p,om;
  double kk,tmax;
  FILE *invGJ,*invGE, *invGG, *invSE, *invSG, *invSJ;
  invGJ = fopen("./data/inverse_j.dat","w");
  invGE = fopen("./data/inverse_e.dat","w");
  invGG = fopen("./data/inverse_g.dat","w");
  invSJ = fopen("./data/spectrum_j.dat","w");
  invSE = fopen("./data/spectrum_e.dat","w");
  invSG = fopen("./data/spectrum_g.dat","w");

  init_var();

  t0   = M_PI/fmax_w;
  beta = M_PI*pow(fmax_w,2);
  cmin = sqrt(2.0f * omega_0 / MU0/ 1.0f);
  cmax = sqrt(2.0f * omega_0 / MU0/ 1.0f);
  courant = 1.0f/cmax/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
  dt = courant*6.f/7.f*0.999f;

  efld_fic();
  init_emfield();
  efld_dip();

  printf("dt        = %e\n", dt);
  printf("t0        = %e\n", t0);
  printf("omega_0   = %e\n", omega_0);
  printf("fmax_w    = %e\n", fmax_w);

  for(n=0;n<N;n++){
    om = 2.f*M_PI/N/dt;
    EX_w[n] = JX_w[n] = 0.f;
    for(k=0;k<N;k++){
      EX_w[n] += EX_f[k] * dt*\
                 exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
      //JX_w[n] += JX_f[k] * dt *\
                 exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
//      EX_w[n] += EX_f[k];
//      JX_w[n] += JX_f[k];
      //JX_w[n] += csqrt(-2.f*omega_0/I/om/n)*JX_f[k] * dt *\
      //           exp(-sqrt(omega_0*om*n)*k*dt)*cexp(I*sqrt(omega_0*om*n)*k*dt);
      //if(n==1000) printf("%d    %e    %e\n",k,creal(JX_f[k]*exp(-sqrt(omega_0*om*n)*k*dt)),exp(-sqrt(omega_0*om*n)*k*dt));
    }
    //JX_w[0] = 2.f*omega_0;
    JX_w[n] = 2.f*omega_0 * exp(-sqrt(omega_0*om*n)*t0) * cexp(I*sqrt(omega_0*om*n)*t0)* \
              cexp(-I*omega_0*om*n/2.f/beta);     // J(x,w)
    //GX_w[n] = -I*(I+1.f)*sqrt(om*n*omega_0)*exp(-sqrt(om*n*omega_0)*t0) * \
              cexp(I*sqrt(om*n*omega_0)*t0) * cexp(-I*om*n*omega_0/2.f/beta); // J'(x,w')
    GX_w[n] = EX_w[n] / JX_w[n];
    //GX_w[n] = EX_w[n] / JX_w[n] * csqrt(-I*om*n/2.f/omega_0);

    printf("Cauction!  om is multipuled by 0.17\n");
    fprintf(invSE,"%10e %10e  %10e  %10e   %10e\n", \
            n*om/2.f/M_PI,creal(EX_w[n]),creal(I*EX_w[n]), pow(creal(EX_w[n])*creal(EX_w[n])+creal(I*EX_w[n])*creal(I*EX_w[n]), 0.5f),  atan(creal(EX_w[n]*I)/creal(EX_w[n])));
    fprintf(invSJ,"%10e %10e  %10e  %10e  %10e\n", \
            n*om/2.f/M_PI,creal(JX_w[n]),creal(I*JX_w[n]), pow(creal(JX_w[n])*creal(JX_w[n])+creal(I*JX_w[n])*creal(I*JX_w[n]), 0.5f), atan(creal(JX_w[n]*I)/creal(JX_w[n])));
    fprintf(invSG,"%10e %10e  %10e  %10e  %10e\n", \
            n*om/2.f/M_PI,creal(GX_w[n]),creal(I*GX_w[n]), pow(creal(GX_w[n])*creal(GX_w[n])+creal(I*GX_w[n])*creal(I*GX_w[n]),0.5f), atan(creal(GX_w[n]*I)/creal(GX_w[n])));
  }

  for(k=0;k<N;k++){
    for(n=0;n<N;n++){
      JX_t[k] += JX_w[n] * cexp(-I*2.f*M_PI*k*n/N) /N/dt*2.f;
      EX_t[k] += EX_w[n] * cexp(-I*2.f*M_PI*k*n/N) /N/dt*2.f;
      GX_t[k] += GX_w[n] * cexp(-I*2.f*M_PI*k*n/N) /N/dt*2.f;
    }
    fprintf(invGJ,"%e    %e   %e\n",k*dt, creal(JX_t[k]),  creal(I*JX_t[k]));
    fprintf(invGE,"%e    %e   %e\n",k*dt, creal(EX_t[k]),  creal(I*EX_t[k]));
    fprintf(invGG,"%e    %e   %e\n",k*dt, creal(GX_t[k]),  creal(I*GX_t[k]));
  }
}
//////////////////////////////////////////
void init_var()
{
    int tmp;
    double dtmp;
    FILE *ifm;

    ifm  = fopen("./data/model_env.dat","r");

    fscanf(ifm,"%d", &nx);
    fscanf(ifm,"%d", &ny);
    fscanf(ifm,"%d", &nz);
    fscanf(ifm,"%lf", &dx);
    fscanf(ifm,"%lf", &dy);
    fscanf(ifm,"%lf", &dz);
    fscanf(ifm,"%lf", &dtmp);
    fscanf(ifm,"%d", &N);
    fscanf(ifm,"%d", &tmp);
    fscanf(ifm,"%d", &tmp);
    fscanf(ifm,"%lf", &omega_0);
    fscanf(ifm,"%lf", &fmax_w);
    fscanf(ifm,"%lf", &dtmp);
    fscanf(ifm,"%d", &tmp);
    fscanf(ifm,"%lf\n",&FL);
    fscanf(ifm,"%lf\n",&FH);
    fscanf(ifm,"%lf\n",&FS);
    fscanf(ifm,"%lf\n",&AP);
    fscanf(ifm,"%lf\n",&AS);

    fclose(ifm);
}
//////////////////////////////////////////
void init_emfield()
{
    int i;
    double dtmp, ReF, ImF;
    FILE *ifp1,*ifp2,*ifp3,*ifw1,*ifw2,*ifw3;

    ifw1 = fopen("./data/waveformX.dat","r");
    ifw2 = fopen("./data/waveformY.dat","r");
    ifw3 = fopen("./data/waveformZ.dat","r");
    //ifp1 = fopen("./data/rex_000_000.dat","r");
    ifp1 = fopen("./data/anal2.dat","r");
    //ifp2 = fopen("./data/rey_000_000.dat","r");
    //ifp3 = fopen("./data/rez_000_000.dat","r");

    EX_f = (double *)malloc(sizeof(double)*N);
    EX_w = (double _Complex *)malloc(sizeof(double _Complex)*N);
    EX_t = (double _Complex *)malloc(sizeof(double _Complex)*N);

    JX_f = (double _Complex *)malloc(sizeof(double _Complex)*N);
    JX_w = (double _Complex *)malloc(sizeof(double _Complex)*N);
    JX_t = (double _Complex *)malloc(sizeof(double _Complex)*N);

    GX_t = (double _Complex *)malloc(sizeof(double _Complex)*N);
    GX_w = (double _Complex *)malloc(sizeof(double _Complex)*N);

    ReF=ImF=0.f;
    for(i=0;i<N;i++){
        fscanf(ifp1,"%lf   %lf", &dtmp, &ReF);
        EX_f[i] = ReF;
        fscanf(ifw1,"%lf   %lf", &dtmp, &ReF);
        JX_f[i] = ReF;
    }

    fclose(ifp1); //fclose(ifp2); fclose(ifp3);
    fclose(ifw1); fclose(ifw2); fclose(ifw3);
}
//////////////////////////////////////////
void efld_dip()
{
    int  i,t;
    double sx=0.f, sy=0.f, sz=0.f;
    double rx=2000.f, ry=0.f, rz=0.f;
    double moment = 1.0;
    double sig = 1.f;
    double coef1, coef2, theta,time,r,theta_r, erfc, erf;
    double EX,EY,EZ,EX_p, d_EX;
    FILE *of1;
    setbuf(stdout,NULL);

    of1 = fopen("./data/anal1.dat","wt");

    r = sqrt(pow((sx-rx),2)+pow((sy-ry),2)+pow((sz-rz),2));
    printf("Target offset   =  %f  [m]\n",r);

    printf("\n => ");
    EX=EX_p=0.f;
    fprintf(of1,"%e   %e   %e   %e   %e\n", 0.f, 0.f, 0.f,0.f,0.f);
    for(t=1;t<N;t++){
      time = t*dt;
      theta = sqrt(MU0*sig/4.f/time);
      theta_r = theta*r;
      erfc = 1.f-erff(theta_r);
      coef1 = moment/(4.f*M_PI*sig*pow(r,3)) * \
              ((4.f/sqrt(M_PI)*pow(theta_r,3) + 6.f/sqrt(M_PI)*theta_r)*exp(-pow(theta_r,2)) +3.f*erfc);
      coef2 = moment/(4.f*M_PI*sig*pow(r,3)) * \
              ((4.f/sqrt(M_PI)*pow(theta_r,3) + 2.f/sqrt(M_PI)*theta_r)*exp(-pow(theta_r,2)) +erfc);

      EX = coef1 * pow(rx/r,2) - coef2*1.f;
      EY = coef1 * rx*ry/r/r;
      EZ = coef1 * rx*rz/r/r;

      d_EX = (EX- EX_p)/dt;
      EX_p = EX;
      if(t%(N/20)==N/20-1) printf("#");
      fprintf(of1,"%e   %e   %e   %e   %e\n", time, d_EX, EX,erfc,sqrt(M_PI)*theta_r)*exp(-pow(theta_r,2));
    }
    printf("\n");
    fclose(of1);
}
//////////////////////////////////////////
void efld_fic()
{
    int t;
    double sx=0.f, sy=0.f, sz=0.f;
    double rx=200.f, ry=0.f, rz=0.f;
    double moment = 1.0;
    double sig = 1.f;
    double coef1, coef2, time,r,t1, t0,beta;
    double EX,gamma1, d_gamma, d2_gamma;
    FILE *of2;

    of2=fopen("./data/anal2.dat","wt");

    beta=M_PI*pow(fmax_w,2);
    t0  =M_PI/fmax_w;

    r = sqrt(pow((sx-rx),2)+pow((sy-ry),2)+pow((sz-rz),2));

    EX=0.f;
    for(t=0;t<N;t++){
        time= t*dt;
        t1 = time-r/cmax;
        coef1 = (MU0/4.f/M_PI/r) * (pow(rx/r,2) -1.f);
        coef2 = (MU0/4.f/M_PI/r) * (3.f*pow(rx/r,2)-1.f);

        gamma1    = sqrt(beta/M_PI)*exp(-beta*pow(t1-t0,2));
        d_gamma  = -2.f*beta*(t1-t0)*sqrt(beta/M_PI)*exp(-beta*pow(t1-t0,2));
        d2_gamma = sqrt(beta/M_PI)*exp(-beta*pow(t1-t0,2))*(4.f*pow(beta,2)*pow(t1-t0,2)-2.f*beta);

        EX = coef1*d2_gamma + coef2*(cmax/r*d_gamma + pow(cmax/r,2)*gamma1);
        fprintf(of2,"%e    %e\n",time,EX);
    }
    fclose(of2);
}

////////////////////////////////////
int maxin3(int dx,int dy,int dz)
{
    int tmp;
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
