#include <stdlib.h>
#include "fftw3.h"
#include <stdio.h>
#include <math.h>
#define velocity 2.99790000e8
#define MU0 12.5663706e-7
#define PI 3.14159265358979
#define HP (PI/2.0)
#define sw 0
void read_para(void);
void init_bandpass(void);
void firstderiv_gauss(int nt,double *signal,double dt);
void zero(int, double *, double);
int sinw(int tt,double *signal,double dt);
void output(double *,double *,double *,double *,int);
int butpas(double *h,int *m,double *gn,int *n,double fl,double fh,double fs,double ap,double as);
int recfil(double *x,double *y,int n,double *h,int nml,double *uv);
int tandem(double *x,double *y,int n,double *h,int m,int nml,double *uv);
int ricker(int tt, double *signal,double dt);
int sin2(int tt, double *signal,double dt);
int gauss_sin(int tt,double *signal,double dt);
int fourier1(int tt, double *signal, double dt);

int nx,ny,nz,nt;
double dt,courant,courant_ori,dx,dy,dz,cmax,cmin,f_max,f_min;
double FL,FH,FS,AP,AS;
double omega_0,fmax_w,tmax, dt_ratio;
double *signalX,*signalY,*signalZ,*signalTrue;
double *after,*past;

int main()
{
  int M,N,i;
  double G0;
  double *keisu,*hkeisu;

  read_para();
  init_bandpass();

  firstderiv_gauss(nt,signalX, dt);
//  gauss_sin(nt,signalX,dt);
  zero(nt,signalY,dt);
  zero(nt,signalZ,dt);

  //sin2(nt, signalTrue, dt);
  ricker(nt, signalTrue, dt);

// BandPass filter  /////////////
//  hkeisu=(double *)calloc(4, sizeof(double));
//  butpas(hkeisu,&M,&G0,&N,FL*dt,FH*dt,FS*dt,AP,AS);
//  keisu=(double *)calloc(M*4, sizeof(double));
//  butpas(keisu,&M,&G0,&N,FL*dt,FH*dt,FS*dt,AP,AS);
//  tandem(signalX,after,nt,keisu,M,1,past);
//  tandem(after,after,nt,keisu,M,-1,past);

//  if(sw == 1){
//  for(i=0;i<nt;i++){
//      signalX[i] = after[i]*G0*G0;
//  }
//  }

  fourier1(nt,signalX,dt);
  output(signalX,signalY,signalZ,signalTrue,nt);
}
//********************************************************
void read_para(void)
{
  int tmp,i;
  double Glim, minv,maxv;
  FILE *ifp1,*ifp2;
  ifp1 = fopen("./data/model_env.dat","rt");
  ifp2 = fopen("./data/minmax.dat","rt");

  fscanf(ifp1,"%d", &nx);
  fscanf(ifp1,"%d", &ny);
  fscanf(ifp1,"%d", &nz);
  fscanf(ifp1,"%lf", &dx);
  fscanf(ifp1,"%lf", &dy);
  fscanf(ifp1,"%lf", &dz);
  fscanf(ifp1,"%lf", &dt_ratio);
  fscanf(ifp1,"%d",  &nt);
  fscanf(ifp1,"%d",  &tmp);
  fscanf(ifp1,"%d",  &tmp);
  fscanf(ifp1,"%lf", &omega_0);
  fscanf(ifp1,"%lf", &fmax_w);
  fscanf(ifp1,"%lf", &Glim);
  fscanf(ifp1,"%d", &tmp);
  fscanf(ifp1,"%lf", &FL);
  fscanf(ifp1,"%lf", &FH);
  fscanf(ifp1,"%lf", &FS);
  fscanf(ifp1,"%lf", &AP);
  fscanf(ifp1,"%lf", &AS);

  signalX = (double *)malloc( sizeof(double)*nt );
  signalY = (double *)malloc( sizeof(double)*nt );
  signalZ = (double *)malloc( sizeof(double)*nt );
  signalTrue = (double *)malloc( sizeof(double)*nt );
  after   = (double *)malloc( sizeof(double)*nt );
  past    = (double *)malloc( sizeof(double)*nt );
  for(i=0;i<nt;i++) signalX[i]=signalY[i]=signalZ[i]=after[i]=past[i]=0.f;

  //courant = 1.0f/velocity/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  //dt = courant *dt_ratio;
  //courant_ori = 1.0f/velocity/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  fscanf(ifp2, "%lf\n", &minv);
  fscanf(ifp2, "%lf\n", &maxv);
  cmin = sqrt(2.f*omega_0/MU0/maxv);
  cmax = sqrt(2.f*omega_0/MU0/minv);
  courant = 1.0f/cmax/sqrt(1.0f/pow(dx,2) +1.0f/pow(dy,2) + 1.0f/pow(dz,2));
  dt = courant*6.f/7.f*0.999f;

  fclose(ifp1);
  fclose(ifp2);
}

/////////////////////////////////////////////////
void init_bandpass(void){
    f_min = 2.f*M_PI/nt/dt;
    f_max = M_PI/dt;
    //FL   = 2.f*f_min;
    //FH   = f_max*0.05;
    //FS   = FH*1.5f;
    //FL   = 4.f;
    //FH   = 8.f;
    //FS   = 12.f;
    //AP   = 0.5f;
    //AS   = 5.f;
    printf("fl=%14.6e   fh=%14.6e   fs=%14.6e\n",FL,FH,FS);
}

//////////////////////////////////////////////////////
void firstderiv_gauss(int nt,double *signal,double dt)
{
  int i;
  double beta = M_PI*pow(fmax_w,2);
  double t0   = M_PI/fmax_w;

  for(i=0;i<nt;i++){
      signal[i] = -2.0f*beta*(i*dt-t0)*sqrt(beta/M_PI)* \
                  exp(-beta*pow(i*dt - t0,2.f));
  }

}
//////////////////////////////////////////////////////
void zero(int tt, double *signal, double dt)
{
  int i;
  for(i=0;i<tt;i++){
    signal[i] = 0.f;
  }
}

//////////////////////////////////////////////////////
void output(double *signalX,double *signalY,double *signalZ,double *signalTrue,int tt)
{
  FILE *ofp1;
  FILE *ofp2;
  FILE *ofp3;
  FILE *ofp4;
  int i;
  ofp1 = fopen("./data/waveformX.dat","wt");
  ofp2 = fopen("./data/waveformY.dat","wt");
  ofp3 = fopen("./data/waveformZ.dat","wt");
  ofp4 = fopen("./data/true_waveformX.dat","wt");
  for(i=0;i<nt;i++){
    fprintf(ofp1, "%e  %e\n",i*dt, signalX[i]);
    fprintf(ofp2, "%e  %e\n",i*dt, signalY[i]);
    fprintf(ofp3, "%e  %e\n",i*dt, signalZ[i]);
    fprintf(ofp4, "%e  %e\n",i*dt, signalTrue[i]);
  }

  printf("cmax     = %e\n",  cmax);
  printf("cmin     = %e\n",  cmin);
  printf("tmax     = %lf\n", tmax);
  printf("dx       = %lf\n", dx);
  printf("dy       = %lf\n", dy);
  printf("dz       = %lf\n", dz);
  printf("dt_ratio = %lf\n", dt_ratio);
  printf("tt       = %d\n", nt);
  printf("delta t = %e\n",dt);
  printf("courant_ori = %e\n",courant_ori);
  printf("fmax of wave = %e\n",fmax_w);
  printf("Total calculation time = %e\n", dt*nt);
  printf("freqency band = %4e ~ %4e  \n", f_min, f_max);
  fclose(ofp1);fclose(ofp2);fclose(ofp3);fclose(ofp4);
  free(signalX);
  free(signalY);
  free(signalZ);
  free(signalTrue);
  free(after);
  free(past);
}
/////////////////////////////////////////////////
//Making Coefficients for Band-Pass Filter****
int butpas(double *h,int *m,double *gn,int *n,double fl,double fh,double fs,double ap,double as)
{
  double wl,wh,ws,clh,op,ww,ts,os,pa,sa,cc,c,dp,g,fj,rr,tt,re,ri,a,wpc,wmc;
  int k,l,j,i;
  struct {
    double r;
    double c;
  } oj,aa,cq,r[2];
  if(fabs(fl)<fabs(fh)) wl=fabs(fl)*M_PI;
  else wl=fabs(fh)*M_PI;
  if(fabs(fl)>fabs(fh)) wh=fabs(fl)*M_PI;
  else wh=fabs(fh)*M_PI;
  ws=fabs(fs)*M_PI;
  if(wl==0.0 || wl==wh || wh>=HP || ws==0.0 || ws>=HP ||
     (ws-wl)*(ws-wh)<=0.0)
    {
      fprintf(stderr,"? (butpas) invalid input : fl=%14.6e fh=%14.6e fs=%14.6e ?\n",
	      fl,fh,fs);
      *m=0;
      *gn=1.0;
      return 1;
    }
  /****  DETERMINE N & C */
  clh=1.0/(cos(wl)*cos(wh));
  op=sin(wh-wl)*clh;
  ww=tan(wl)*tan(wh);
  ts=tan(ws);
  os=fabs(ts-ww/ts);
  if(fabs(ap)<fabs(as)) pa=fabs(ap);
  else pa=fabs(as);
  if(fabs(ap)>fabs(as)) sa=fabs(ap);
  else sa=fabs(as);
  if(pa==0.0) pa=0.5;
  if(sa==0.0) sa=5.0;
  if((*n=(int)(fabs(log(pa/sa)/log(op/os))+0.5))<2) *n=2;
  cc=exp(log(pa*sa)/(double)(*n))/(op*os);
  c=sqrt(cc);
  ww=ww*cc;
  
  dp=HP/(double)(*n);
  k=(*n)/2;
  *m=k*2;
  l=0;
  g=fj=1.0;
  
  for(j=0;j<k;j++)
    {
      oj.r=cos(dp*fj)*0.5;
      oj.c=sin(dp*fj)*0.5;
      fj=fj+2.0;
      aa.r=oj.r*oj.r-oj.c*oj.c+ww;
      aa.c=2.0*oj.r*oj.c;
      rr=sqrt(aa.r*aa.r+aa.c*aa.c);
      tt=atan(aa.c/aa.r);
      cq.r=sqrt(rr)*cos(tt/2.0);
      cq.c=sqrt(rr)*sin(tt/2.0);
      r[0].r=oj.r+cq.r;
      r[0].c=oj.c+cq.c;
      r[1].r=oj.r-cq.r;
      r[1].c=oj.c-cq.c;
      g=g*cc;
      
      for(i=0;i<2;i++)
	{
	  re=r[i].r*r[i].r;
	  ri=r[i].c;
	  a=1.0/((c+ri)*(c+ri)+re);
	  g=g*a;
	  h[l  ]=0.0;
	  h[l+1]=(-1.0);
	  h[l+2]=2.0*((ri-c)*(ri+c)+re)*a;
	  h[l+3]=((ri-c)*(ri-c)+re)*a;
	  l=l+4;
	}
    }
  /****  EXIT */
  *gn=g;
  if(*n==(*m)) return 0;
  /****  FOR ODD N */
  *m=(*m)+1;
  wpc=  cc *cos(wh-wl)*clh;
  wmc=(-cc)*cos(wh+wl)*clh;
  a=1.0/(wpc+c);
  *gn=g*c*a;
  h[l  ]=0.0;
  h[l+1]=(-1.0);
  h[l+2]=2.0*wmc*a;
  h[l+3]=(wpc-c)*a;
  return 0;
}
//*****************************************************
//Recursive Filter*************************************
int recfil(double *x,double *y,int n,double *h,int nml,double *uv)
{
  /*
  +   RECURSIVE FILTERING : F(Z) = (1+A*Z+AA*Z**2)/(1+B*Z+BB*Z**2)
  +
  +   ARGUMENTS
  +   X : INPUT TIME SERIES
  +   Y : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
  +   N : LENGTH OF X & Y
  +   H : FILTER COEFFICIENTS ; H(1)=A, H(2)=AA, H(3)=B, H(4)=BB
  +   NML : >0 ; FOR NORMAL  DIRECTION FILTERING
  +       <0 ; FOR REVERSE DIRECTION FILTERING
  +   uv  : past data and results saved
  +
  +   M. SAITO  (6/XII/75)
*/
  int i,j,jd;
  double a,aa,b,bb,u1,u2,u3,v1,v2,v3;
  if(n<=0)
    {
      fprintf(stderr,"? (recfil) invalid input : n=%d ?\n",n);
      return 1;
    }
  if(nml>=0)
    {
      j=0;
      jd=1;
    }
  else
    {
      j=n-1;
      jd=(-1);
    }
  a =h[0];
  aa=h[1];
  b =h[2];
  bb=h[3];
  u1=uv[0];
  u2=uv[1];
  v1=uv[2];
  v2=uv[3];
  /****  FILTERING */
  for(i=0;i<n;i++)
    {
      u3=u2;
      u2=u1;
      u1=x[j];
      v3=v2;
      v2=v1;
      v1=u1+a*u2+aa*u3-b*v2-bb*v3;
      y[j]=v1;
      j+=jd;
    }
  uv[0]=u1;
  uv[1]=u2;
  uv[2]=v1;
  uv[3]=v2;
  return 0;
}
//********************************************************
//Calling Recursive Filter********************************
int tandem(double *x,double *y,int n,double *h,int m,int nml,double *uv)
{
  /*
  +   RECURSIVE FILTERING IN SERIES
  +
  +   ARGUMENTS
  +   X : INPUT TIME SERIES
  +   Y : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
  +   N : LENGTH OF X & Y
  +   H : COEFFICIENTS OF FILTER
  +   M : ORDER OF FILTER
  +   NML : >0 ; FOR NORMAL  DIRECTION FILTERING
  +       <0 ;   REVERSE DIRECTION FILTERING
  +   uv  : past data and results saved
  +
  +   SUBROUTINE REQUIRED : RECFIL
  +
  +   M. SAITO  (6/XII/75)
*/
  int i;
  if(n<=0 || m<=0)
    {
      fprintf(stderr,"? (tandem) invalid input : n=%d m=%d ?\n",n,m);
      return 1;
    }
  /****  1-ST CALL */
  recfil(x,y,n,h,nml,uv);
  /****  2-ND AND AFTER */
  if(m>1) for(i=1;i<m;i++) recfil(y,y,n,&h[i*4],nml,&uv[i*4]);
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
  double fp = 1.e+0;
  double Md = 1.0f;
  double phase = Md/fp;

  for(i=0;i<tt;i++){
    alpha = pow(M_PI*fp,2)*pow((double)(i-200)*dt - phase,2);
    beta  = 2.0f*alpha - 1.0f;
    signal[i] = beta * exp(-alpha);

//    printf("%e %e %e\n",alpha,beta,signal[i]);
  }

  return 0;
}

/////////////////////////////////////////////////
int sin2(int tt,double *signal,double dt)
{
    int i;
    double freq = 1.06f;
    double om   = 2.f*M_PI*freq;
    double Tw   = 4.f*M_PI/om;
    int    iTw  = (int)(Tw/dt);

    for(i=0;i<iTw;i++){
        signal[i] = 0.5f*(1.f - cos(M_PI*i*dt/Tw))*sin(om*i*dt);
    }

    for(i=iTw;i<tt;i++){
        signal[i] = sin(om*i*dt);
    }
    return 0;
}

/////////////////////////////////////////////////
int gauss_sin(int tt,double *signal,double dt)
{
    int i;
    double om   = 10.0f;
    double t0 = M_PI/om;
    double alpha = pow(2.f/t0,2.f);
    double wc   = 1.f*M_PI/2.f/t0;
    int   it0 = (int)(2.f*t0/dt);

    for(i=0;i<it0;i++){
        signal[i] = exp(-alpha*pow(i*dt-t0,2.f)) * cos(wc*(i*dt - t0));
    }

    for(i=2.f*t0;i<tt;i++){
        signal[i] = 0.f;
    }
    return 0;
}

/////////////////////////////////////////////////
int fourier1(int tt, double *signal, double dt)
{
  int n,k;
  double ReF,ImF,om;
  FILE *fourier;
  fourier=fopen("./data/fourier.data","w");


  om = 1.f/tt/dt;
  for(n=0;n<tt;n++) {
    ReF=ImF=0.0;
    for(k=0;k<tt;k++) {
      ReF+=signal[k]*cos(2*M_PI*k*n/tt);
      ImF+=signal[k]*sin(2*M_PI*k*n/tt);
    }
    fprintf(fourier,"%f %e %e %e\n",n*om,ReF,ImF,ReF*ReF+ImF*ImF);             //最後にスペクトルの強さを追加
  }

  return 0;
}

