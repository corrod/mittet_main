#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#define MU0 12.5663706e-7
void init_var();
void init_emfield();
void fourier_trans(FILE *invJT,FILE *invJF,FILE *invJW);
void inv_fourier(FILE *invJT,FILE *invJF,FILE *invJW);

int N,nx,ny,nz;
double dx,dy,dz,dt,omega_0,fmax_w,vc,cmax,cmin,courant;
double _Complex *JX_f,*JX_w,*JX_t;

int main()
{
    int i,j,k;
    double dtmp;
    FILE *invJT,*invJF,*invJW;
    invJT = fopen("./data/transtime_j.dat","w");
    invJF = fopen("./data/transfreq_j.dat","w");
    invJW = fopen("./data/transomeg_j.dat","w");

    init_var();
    init_emfield();

    fourier_trans(invJT,invJF,invJW);
    inv_fourier(invJT,invJF,invJW);
}

//////////////////////////////////////////
void init_var()
{
    int tmp;
    double dtmp;
    double minv,maxv;
    FILE *ifm, *ifm2;

    ifm  = fopen("./data/model_env.dat","r");
    ifm2 = fopen("./data/minmax.dat","r");

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

    fscanf(ifm2, "%lf\n",&minv);
    fscanf(ifm2, "%lf\n",&maxv);

    fclose(ifm);
    fclose(ifm2);

    cmin = sqrt(2.0f * omega_0 / MU0/ maxv);
    cmax = sqrt(2.0f * omega_0 / MU0/ minv);
    courant = 1.0f/cmax/sqrt(1.0/pow(dx,2) + 1.0/pow(dy,2) + 1.0/pow(dz,2));
    dt = courant*6.f/7.f*0.999f;
}
//////////////////////////////////////////
void init_emfield()
{
    int i;
    double dtmp, ReF, ImF;
    FILE *ifw1,*ifw2,*ifw3;

    ifw1 = fopen("./data/waveformX.dat","r");
    ifw2 = fopen("./data/waveformY.dat","r");
    ifw3 = fopen("./data/waveformZ.dat","r");

    JX_t = (double _Complex *)malloc(sizeof(double _Complex)*N); // diffusive-time
    JX_f = (double _Complex *)malloc(sizeof(double _Complex)*N); // diffusive-freq
    JX_w = (double _Complex *)malloc(sizeof(double _Complex)*N); // fictitious-freq

    ReF=ImF=0.f;
    for(i=0;i<N;i++){
        fscanf(ifw1,"%lf   %lf", &dtmp, &ReF);
        JX_t[i] = ReF;
    }

    fclose(ifw1); fclose(ifw2); fclose(ifw3);
}
//////////////////////////////////////////
void fourier_trans(FILE *inv_JT,FILE *inv_JF,FILE *inv_JW)
{
    int n,tt;
    double om;

    for(n=0;n<N;n++){
      om = 2.f*M_PI/N/dt;
      JX_f[n] = 0.f;
      for(tt=0;tt<N;tt++){
        JX_f[n] += JX_t[tt] * cexp(I*2.f*M_PI*tt*n/N);
      }
      JX_w[n] = csqrt(-I*om*n/2.f/omega_0) * JX_f[n];


      fprintf(inv_JF, "%lf     %e      %e     %e    %e\n", \
              om*n, creal(JX_f[n]), creal(I*JX_f[n]), creal(JX_f[n])*creal(JX_f[n])+creal(I*JX_f[n])*creal(I*JX_f[n]),atan(creal(JX_f[n]*I)/creal(JX_f[n])) );
      fprintf(inv_JW, "%lf     %e      %e     %e    %e\n", \
              om*n, creal(JX_w[n]), creal(I*JX_w[n]), creal(JX_w[n])*creal(JX_w[n])+creal(I*JX_w[n])*creal(I*JX_w[n]),atan(creal(JX_w[n]*I)/creal(JX_w[n])) );
    }

}
//////////////////////////////////////////
void inv_fourier(FILE *inv_JT,FILE *inv_JF,FILE *inv_JW)
{
    int n,k;
    double om;

    for(k=0;k<N;k++){
      om = 2.f*M_PI/N/dt;
      for(n=0;n<N;n++){
          JX_t[k] += (I+1)*csqrt(-I/2.f)/4.f/M_PI*JX_f[n] \
                     *cexp(sqrt(om*n*omega_0)*k*dt)*cexp(-I*sqrt(om*n*omega_0)*k*dt);
      }
      fprintf(inv_JT,"%lf     %e    %e\n",dt*k, creal(JX_t[k]),creal(I*JX_t[k]));
    }

}
