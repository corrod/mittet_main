void read_shotrec(void);
void read_trawave(void);
void init_eh_field_3d(void);
void media_coeff_3d(void);
void init_eh_mur_3d(void);
void init_FILE(int isource);
void close_FILE(void);
void set_zero(void);
void e_field(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void read_source_3d(double *EX_r,double *signalX_r,int isource,int step);
void h_field(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void e_field4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void e_field42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void h_field4(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void h_field42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void output(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ, \
        FILE **ofe1,FILE **ofe2,FILE **ofe3,FILE **ofh1,FILE **ofh2,FILE **ofh3,    \
        int step);
void checking_omp(void);
void init_gradient(void);
void init_FILE2(int isource,int iter);
void read_Eobs(void);
void set_zero_eh(void);
void copytoarray(double *EX,double *EY,double *EZ,int step,FILE **ofe1,FILE **ofe2,FILE **ofe3);
void calc_delE(int isource,int iter);
void e_field4_bp(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void h_field4_bp(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void read_backwave_3d(double *EX,double *delEx,int step);
void e_field_cpml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void h_field_cpml(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void e_field_cpml42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void h_field_cpml42(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void e_field_cpml422(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
void h_field_cpml422(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ);
