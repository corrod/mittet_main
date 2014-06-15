#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define CC 2.99790000e8
#define MU0 12.5663706e-7

int firstderiv_gauss(int, double *, double);
int gausspulse(int,double *,double);
int ricker(int,double *,double);
int kukei(int,double *,double);
int sinw(int,double *,double);
int sinc(int,double *,double);
int zero(int,double *,double);
int transform(void);
int read_para(void);
int output(void);
int maxin3(int,int,int);

int i,tt,tmp,nx,ny,nz,nmax;
double dt,courant,courant_ori;
double kk;
double dx,dy,dz,dt_ratio,vc,cmax,cmin;
double omega_0, fmax_w,tmax,Glim,t0;
double *signalX;
double *signalY;
double *signalZ;

int main()
{
  read_para();
  //t0   = M_PI/fmax_w;

  //courant = 1.0f/CC/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  //dt = courant *dt_ratio;
  //courant_ori = 1.0f/CC/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  cmin = sqrt(2.f*omega_0/MU0/1.f);
  cmax = sqrt(2.f*omega_0/MU0/1.f);
  courant = 1.0f/cmax/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  dt = courant*0.8f;

  //fmax_w = 1.f;
  //fmax_w = cmin/Glim/dx;
//  nmax   = maxin3(nx,ny,nz);
//  tmax   = 1.0*nmax*dx/cmin;
//  tt     = (int)((tmax+2.f*t0) / dt);

// @ input transmitter signal //////
  firstderiv_gauss(tt,signalX, dt);
  //sinw(tt,signalX,dt);
  //gausspulse(tt,signalX,dt);
  zero(tt,signalY,dt);
  zero(tt,signalZ,dt);

// @ transform signal to fictitious domain///////
  //transform();

// @ output transmitter signal ///////
  output();

}

//////////////////////////////////////////////////////
int firstderiv_gauss(int tt,double *signal,double dt)
{
  int i;
  double beta = M_PI*pow(fmax_w,2);
  double t0   = M_PI/fmax_w;
  FILE *of3;

  of3 = fopen("./data/waveform1.dat","wt");
  for(i=0;i<tt;i++){
      signal[i] = -2.0f*beta*(i*dt-t0)*sqrt(beta/M_PI)* \
                  exp(-beta*pow(i*dt - t0,2.f));
      fprintf(of3, "%e    %e\n",i*dt, signal[i]);
  }

  fclose(of3);
  return 0;
}
//////////////////////////////////////////////////////
int gausspulse(int tt,double *signal, double dt)
{
  int i;
  double tau0 = (int)(tt/10);
  double alpha = pow(4.f/tau0,2);

  for(i=0;i<tt;i++){
    if(i<tt/5) signal[i]=exp(-alpha * pow((double)i-tau0,2));
    if(tt/5<i && i<tt) signal[i] = -exp(-alpha*pow((double)(i-(int)(tt/5))-tau0,2));
  }
  return 0;
}
//////////////////////////////////////////////////////
int zero(int tt, double *signal, double dt)
{
  int i;
  for(i=0;i<tt;i++){
    signal[i] = 0.f;
  }
  return 0;
}

//////////////////////////////////////////////////////
int sinw(int tt,double *signal,double dt)
{
    int i;
    double freq=100.f;

    for(i=0;i<tt;i++){
        signal[i] = sin(2.0*M_PI*freq*i*dt);
    }
    return 0;
}

//////////////////////////////////////////////////////
int sinc(int tt,double *signal,double dt)
{
    int i;
    double freq=1000.f;

    for(i=0;i<tt;i++){
        signal[i] = sin(2.0*M_PI*freq*(i-tt/2)*dt)/(2.0*M_PI*freq*i*dt);
    }
    return 0;
}
//////////////////////////////////////////////////////
int kukei(int tt, double *signal,double dt)
{
  int i;
  int fp = 1500;
  int fp_r;

  for(i=0;i<tt;i++){
    fp_r=i%(fp*4);
    if(fp_r < fp) signal[i]=0.0f;
    if(fp*1 <= fp_r && fp_r < fp*2) signal[i]=-1.0f;
    if(fp*2 <= fp_r && fp_r < fp*3) signal[i]= 0.0f;
    if(fp*3 <= fp_r && fp_r < fp*4) signal[i]=+1.0f;
  }

  return 0;
}
/////////////////////////////////////////////////
int ricker(int tt, double *signal,double dt)
{
  int i;
  double alpha;
  double beta;
  double fp = 20.e+0;
  //double fp = 100.e+0;
  double Md = 1.0f;
  double phase = Md/fp;

  for(i=0;i<tt;i++){
    alpha = pow(M_PI*fp,2)*pow((double)(i-200000)*dt - phase,2);
    beta  = 2.0f*alpha - 1.0f;
    signal[i] = beta * exp(-alpha);

//    printf("%e %e %e\n",alpha,beta,signal[i]);
  }

  return 0;
}

/////////////////////////////////////////////////
int read_para()
{
  FILE *ifp1;
  ifp1 = fopen("./data/model_env.dat","rt");

  fscanf(ifp1,"%d", &nx);
  fscanf(ifp1,"%d", &ny);
  fscanf(ifp1,"%d", &nz);
  fscanf(ifp1,"%lf", &dx);
  fscanf(ifp1,"%lf", &dy);
  fscanf(ifp1,"%lf", &dz);
  fscanf(ifp1,"%lf", &dt_ratio);
  fscanf(ifp1,"%d",  &tt);
  fscanf(ifp1,"%d",  &tmp);
  fscanf(ifp1,"%d",  &tmp);
  fscanf(ifp1,"%lf", &omega_0);
  fscanf(ifp1,"%lf", &fmax_w);
  fscanf(ifp1,"%lf", &Glim);

  signalX= (double *)malloc( sizeof(double)*tt );
  signalY= (double *)malloc( sizeof(double)*tt );
  signalZ= (double *)malloc( sizeof(double)*tt );
  for(i=0;i<tt;i++) signalX[i]=signalY[i]=signalZ[i]=0.f;
  fclose(ifp1);
  return 0;
}

/////////////////////////////////////////////////
int output()
{
  FILE *ofp1,*ofp2,*ofp3;
  ofp1 = fopen("./data/waveformX.dat","wt");
  ofp2 = fopen("./data/waveformY.dat","wt");
  ofp3 = fopen("./data/waveformZ.dat","wt");
  for(i=0;i<tt;i++){
    fprintf(ofp1, "%e  %e\n",i*dt, signalX[i]);
    fprintf(ofp2, "%e  %e\n",i*dt, signalY[i]);
    fprintf(ofp3, "%e  %e\n",i*dt, signalZ[i]);
  }

  printf("cmax     = %e\n",  cmax);
  printf("cmin     = %e\n",  cmin);
  printf("tmax     = %lf\n", tmax);
  printf("dx       = %lf\n", dx);
  printf("dy       = %lf\n", dy);
  printf("dz       = %lf\n", dz);
  printf("dt_ratio = %lf\n", dt_ratio);
  printf("tt       = %d\n", tt);
  printf("t0       = %e\n", t0);
  printf("delta t = %lf\n",dt);
  printf("courant_ori = %e\n",courant_ori);
  printf("fmax of wave = %e\n",fmax_w);
  printf("Total calculation time = %e\n", dt*tt);
  printf("freqency band = %4e ~ %4e  \n", 2.f*M_PI/tt/dt, M_PI/dt);
  fclose(ofp1);fclose(ofp2);fclose(ofp3);
  free(signalX);
  free(signalY);
  free(signalZ);
  return 0;
}
/////////////////////////////////////////////////
//int transform()
//{
//    int k,n;
//    double om, *Re,*Im;
//    double _Complex *J_func, *J_freq;
//    FILE *of1,*of2;
//
//    of1 =fopen("./data/waveform2.dat","wt");
//    of2 =fopen("./data/waveform3.dat","wt");
//
//    Re = (double *)malloc(sizeof(double)*tt);
//    Im = (double *)malloc(sizeof(double)*tt);
//    J_freq=(double _Complex *)malloc(sizeof(double _Complex)*tt);
//    J_func=(double _Complex *)malloc(sizeof(double _Complex)*tt);
//
//    for(n=0;n<tt;n++){
//      J_freq[n]=0.f;
//      for(k=0;k<tt;k++){
//        J_freq[n] += signalX[k]*cexp(I*2.f*M_PI*k*n/tt);
//      }
//      fprintf(of2,"%d   %10e   %10e\n",n,creal(J_freq[n]),creal(J_freq[n]*I));
//    }
//
//    for(n=0;n<tt;n++){
//      om = 2.f*M_PI/tt/dt;
//      J_freq[n] = csqrt(-I*om*n/2.f/omega_0)*J_freq[n];
//    }
//
//    for(k=0;k<tt;k++){
//      J_func[k]=0.f;
//      for(n=0;n<tt;n++){
//        J_func[k] += J_freq[n]*cexp(-I*2.f*M_PI*k*n/tt)/tt;
//      }
//      fprintf(of1,"%f   %10e   %10e\n",k*dt,creal(J_func[k]),creal(J_func[k]*I));
//      signalX[k] = J_func[k];
//    }
//
//    fclose(of1);
//    fclose(of2);
//    free(J_freq);  free(J_func);
//  return 0;
//}
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
