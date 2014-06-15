// physical constant
#define eps0 8.8541878E-12
#define MU0 12.5663706e-7
#define SIG8 5.7e8
#define CC 2.99790000e8

#define RMAX 40
#define MAXITR 40
#define ELIMIT 0.005

// time condition
unsigned int step;  // step number
unsigned int it;    // total iteration step
double dt,dtime,tt;
double courant,courant_ori;

// transmitter condition
double *signalX;
double *signalY;
double *signalZ;
double *signalTrue;

// space condition
int ix,iy,iz;
int mx,my,mz;
int n;
int nx,ny,nz;
double dx,dy,dz;

// boundary reference
int nxpoint,nypoint,nzpoint;

// 3d field
double *EX, *EY, *EZ;
double *HX, *HY, *HZ;

// coefficient of fdtd calculation
double *cex, *cey, *cez;
double *cex_b, *cey_b, *cez_b;
double *cexry, *cexrz;
double *ceyrx, *ceyrz;
double *cezrx, *cezry;
double *cexry_b, *cexrz_b;
double *ceyrx_b, *ceyrz_b;
double *cezrx_b, *cezry_b;
double *chxry, *chxrz;
double *chyrx, *chyrz;
double *chzrx, *chzry;

// media id
int *id;

// media
int nmedia;
int jcent,kcent;
double *sig, *eps,*mu;
double *sig2;
double sigmax,sigmin, mu1,eps1;
double minv,maxv;
double minval, maxval;
int iseabed;
double *vel;

//// 1st mur
//double *ex1_xy,*ey1_xy;
//double *ey1_yz,*ez1_yz;
//double *ez1_zx,*ex1_zx;
//// 2nd mur
//double *ex2_xy,*ey2_xy;
//double *ey2_yz,*ez2_yz;
//double *ez2_zx,*ex2_zx;
//
//// 1st mur
//double *eyx1, *ezx1;
//double *exy1, *ezy1;
//double *exz1, *eyz1;
//// 2nd mur
//double *eyx2, *ezx2;
//double *exy2, *ezy2;
//double *exz2, *eyz2;

int  shotx, shot_num;
int  recex, rec_num;
FILE *rey_file[RMAX];

int *shot_px;
int *shot_py;
int *shot_pz;
int *rec_px;
int *rec_py;
int *rec_pz;
int isource;

int mstep,trstep;
double dt_ratio;
double vc;
double omega_0, fmax_w;
double cmax,cmin,Glim;

FILE **ofe1, **ofe2, **ofe3;
FILE **ofh1, **ofh2, **ofh3;

// pml
int NPOINT_PML;
double esig_max;
double msig_max;
double *esigx, *msigx;
double *esigy, *msigy;
double *esigz, *msigz;
double *exy, *exz;
double *eyx, *eyz;
double *ezx, *ezy;
double *hxy, *hxz;
double *hyx, *hyz;
double *hzx, *hzy;
double *cexy, *cexz;
double *ceyx, *ceyz;
double *cezx, *cezy;
double *chxy, *chxz;
double *chyx, *chyz;
double *chzx, *chzy;
// cpml
// coefficient 1
double *ca_x, *cb_x;
double *ca_y, *cb_y;
double *ca_z, *cb_z;
double *da_x, *db_x;
double *da_y, *db_y;
double *da_z, *db_z;
// coefficient 2
double *be_x, *bh_x;
double *be_y, *bh_y;
double *be_z, *bh_z;
double *ce_x, *ch_x;
double *ce_y, *ch_y;
double *ce_z, *ch_z;
double *kedx, *kedy, *kedz;
double *khdx, *khdy, *khdz;
double *psi_Eyx1, *psi_Ezx1;
double *psi_Exy1, *psi_Ezy1;
double *psi_Exz1, *psi_Eyz1;
double *psi_Hyx1, *psi_Hzx1;
double *psi_Hxy1, *psi_Hzy1;
double *psi_Hxz1, *psi_Hyz1;
double *ekappax, *mkappax;
double *ekappay, *mkappay;
double *ekappaz, *mkappaz;
double *aex, *amx;
double *aey, *amy;
double *aez, *amz;
int ncpml;

// Laplace transform
double *EX_f;
double *EY_f;
double *EZ_f;
double _Complex *EX_w, *EX_t;
double _Complex *EY_w, *EY_t;
double _Complex *EZ_w, *EZ_t;
double _Complex *JX_f, *JX_w, *JX_t;
double _Complex *JY_f, *JY_w, *JY_t;
double _Complex *JZ_f, *JZ_w, *JZ_t;
double _Complex *GX_w, *GX_t;
double _Complex *GY_w, *GY_t;
double _Complex *GZ_w, *GZ_t;

double *all_recEX;
double *all_recEY;
double *all_recEZ;
double _Complex *all_recGT;
double *err_sum;
int startF,upperF;
