/////////////////////////////////////////////////////////////
//                                                         //
//        3D_FDTD_code_CPU            Ver 1.1              //
//                                                         //
//      * Coded by         N. Imamura  (2012/07/28)        //
//                                                         //
//     v1.0 : basic FDTD is completed  (2012/08/08)        //
//     v1.1 : fictitious FDTD is completed  (2012/09/14)   //
//                                                         //
/////////////////////////////////////////////////////////////
#include "fwi3d_cpu_header1_v1.1.h"

int main(){
  int i,iter,isource,idir,step;
  clock_t start = clock();
  setbuf(stdout,NULL);

/*** Initial Condition ***/
  ix=NXT;iy=NYT;iz=NZT;
  checking_omp();
  read_shotrec();
  lattice_time_3d();
  init_gradient();
  init_eh_field_3d();

/*** Calculating for iteratively ***/
  for(iter=0;iter<MAXITR;iter++){
//    iter = 0;
//    for(idir=0;idir<3;idir++){
    idir = 0;
    media_coeff_3d();

//////////////////////////////////////
//***  Start of FWD calculation  ***//
//////////////////////////////////////

/*** Calculating for each source ***/
    for(isource=0;isource<shot_num;isource++){
      printf("%03d / %03d of %03d th iterative\n",isource,shot_num,iter);
      read_Eobs(); // Read observed value of Ex, Ez
      init_FILE2(isource, iter);
      set_zero_eh();

/*** Calculating for step ***/
      for(step=0;step<it;step++){
        e_field4(EX,EY,EZ,HX,HY,HZ);
//        e_mur2(EX,EY,EZ,HX,HY,HZ);
        cerjan_e(EX,EY,EZ);
        if(idir==0) read_source_3d(EX,signalX,isource,step);
        if(idir==1) read_source_3d(EY,signalY,isource,step);
        if(idir==2) read_source_3d(EZ,signalZ,isource,step);
        h_field4(EX,EY,EZ,HX,HY,HZ);
        cerjan_h(HX,HY,HZ);

        copytoarray(EX,EY,EZ,step,ofe1,ofe2,ofe3);

      if(step%(it/20)==it/20-1) printf("#");
      } // END loop for step
      close_FILE();

      calc_delE(isource,iter);

/////////////////////////////
//***  Start backward   ***//
/////////////////////////////
      set_zero_eh();

/*** Calculating for step ***/
      for(step=0;step<it;step++){
        e_field4_bp(EX,EY,EZ,HX,HY,HZ);
//        e_mur2(EX,EY,EZ,HX,HY,HZ);
        cerjan_e(EX,EY,EZ);
        if(idir==0) read_backwave_3d(EX, delEx,step);
        if(idir==1) read_backwave_3d(EY, delEy,step);
        if(idir==2) read_backwave_3d(EZ, delEz,step);
        h_field4_bp(EX,EY,EZ,HX,HY,HZ);
        cerjan_h(HX,HY,HZ);

        copytoarray(EX,EY,EZ,step,ofe1,ofe2,ofe3);
      }



    } // END loop for isource
  } // END loop for iter




}
