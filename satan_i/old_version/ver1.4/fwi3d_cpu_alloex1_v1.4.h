// physical constant
#define EPS0 8.8541878e-12
#define MU0 12.5663706e-7
#define SIG8 5.7e8
#define CC 2.99790000e8

#define RMAX 40
#define MAXITR 20
#define ELIMIT 0.005

// time condition
extern unsigned int step;  // step number
extern unsigned int it;    // total iteration step
extern double dt,dtime,tt;
extern double courant,courant_ori;

// transmitter condition
extern double *signalX;
extern double *signalY;
extern double *signalZ;
extern double *signalTrue;

// space condition
extern int ix,iy,iz;
extern int mx,my,mz;
extern int n;
extern int nx,ny,nz;
extern double dx,dy,dz;

// boundary reference
extern int nxpoint,nypoint,nzpoint;

// 3d field
extern double *EX, *EY, *EZ;
extern double *HX, *HY, *HZ;

// coefficient of fdtd calculation
extern double *cex, *cey, *cez;
extern double *cex_b, *cey_b, *cez_b;
extern double *cexry, *cexrz;
extern double *ceyrx, *ceyrz;
extern double *cezrx, *cezry;
extern double *cexry_b, *cexrz_b;
extern double *ceyrx_b, *ceyrz_b;
extern double *cezrx_b, *cezry_b;
extern double *chxry, *chxrz;
extern double *chyrx, *chyrz;
extern double *chzrx, *chzry;

// media id
extern int *id;

// media
extern int nmedia;
extern int jcent,kcent;
extern double *sig, *eps,*mu;
extern double *sig2;
extern double sigmax,sigmin, mu1, eps1;
extern double minv,maxv;
extern double minval, maxval;
extern int iseabed;
extern double *vel;

//// 1st mur
//extern double *ex1_xy,*ey1_xy;
//extern double *ey1_yz,*ez1_yz;
//extern double *ez1_zx,*ex1_zx;
//
//// 2nd mur
//extern double *ex2_xy,*ey2_xy;
//extern double *ey2_yz,*ez2_yz;
//extern double *ez2_zx,*ex2_zx;
//
//// 1st mur
//extern double *eyx1, *ezx1;
//extern double *exy1, *ezy1;
//extern double *exz1, *eyz1;
//// 2nd mur
//extern double *eyx2, *ezx2;
//extern double *exy2, *ezy2;
//extern double *exz2, *eyz2;

extern int  shotx, shot_num;
extern int  recex, rec_num;

extern int *shot_px;
extern int *shot_py;
extern int *shot_pz;
extern int *rec_px;
extern int *rec_py;
extern int *rec_pz;
extern int isource;

extern int mstep,trstep;
extern double dt_ratio;
extern double vc;
extern double omega_0, fmax_w;
extern double cmax,cmin,Glim;

extern FILE **ofe1, **ofe2, **ofe3;
extern FILE **ofh1, **ofh2, **ofh3;

// pml
extern int NPOINT_PML;
extern double esig_max;
extern double msig_max;
extern double *esigx, *msigx;
extern double *esigy, *msigy;
extern double *esigz, *msigz;
extern double *exy, *exz;
extern double *eyx, *eyz;
extern double *ezx, *ezy;
extern double *hxy, *hxz;
extern double *hyx, *hyz;
extern double *hzx, *hzy;
extern double *cexy, *cexz;
extern double *ceyx, *ceyz;
extern double *cezx, *cezy;
extern double *chxy, *chxz;
extern double *chyx, *chyz;
extern double *chzx, *chzy;
// cpml
// coefficient 1
extern double *ca_x, *cb_x;
extern double *ca_y, *cb_y;
extern double *ca_z, *cb_z;
extern double *da_x, *db_x;
extern double *da_y, *db_y;
extern double *da_z, *db_z;
// coefficient 2
extern double *be_x, *be_x;
extern double *be_y, *be_y;
extern double *be_z, *be_z;
extern double *bh_x, *bh_x;
extern double *bh_y, *bh_y;
extern double *bh_z, *bh_z;
extern double *ce_x, *ce_x;
extern double *ce_y, *ce_y;
extern double *ce_z, *ce_z;
extern double *ch_x, *ch_x;
extern double *ch_y, *ch_y;
extern double *ch_z, *ch_z;
extern double *kedx, *kedy, *kedz;
extern double *khdx, *khdy, *khdz;
extern double *psi_Eyx1, *psi_Ezx1;
extern double *psi_Exy1, *psi_Ezy1;
extern double *psi_Exz1, *psi_Eyz1;
extern double *psi_Hyx1, *psi_Hzx1;
extern double *psi_Hxy1, *psi_Hzy1;
extern double *psi_Hxz1, *psi_Hyz1;
extern double *ekappax, *mkappax;
extern double *ekappay, *mkappay;
extern double *ekappaz, *mkappaz;
extern double *aex, *amx;
extern double *aey, *amy;
extern double *aez, *amz;
extern int ncpml;

// Laplace transform
extern double *EX_f;
extern double *EY_f;
extern double *EZ_f;
extern double _Complex *EX_w, *EX_t;
extern double _Complex *EY_w, *EY_t;
extern double _Complex *EZ_w, *EZ_t;
extern double _Complex *JX_f, *JX_w, *JX_t;
extern double _Complex *JY_f, *JY_w, *JY_t;
extern double _Complex *JZ_f, *JZ_w, *JZ_t;
extern double _Complex *GX_w, *GX_t;
extern double _Complex *GY_w, *GY_t;
extern double _Complex *GZ_w, *GZ_t;

extern double *all_recEX;
extern double *all_recEY;
extern double *all_recEZ;
extern double _Complex *all_recGT;
extern double *err_sum;
extern int startF,upperF;
