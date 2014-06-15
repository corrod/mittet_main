////////////////////////////////////////////////////////////////////
//                                                                //
//        3D_FDTD_code_CPU_MPI        Ver 1.2                     //
//                                                                //
//      * Coded by         N. Imamura  (2012/07/28)               //
//                                                                //
//     v1.0 : basic FDTD is completed             (2012/08/08)    //
//     v1.1 : fictitious FDTD is completed        (2012/09/14)    //
//     v1.2 : convolutional PML is completed      (2012/10/04)    //
//     v1.3 : accuracy check is completed         (2012/10/17)    //
//     v1.4 : accuracy check is completed         (2012/10/19)    //
//                                                                //
////////////////////////////////////////////////////////////////////
#include "fwi3d_cpu_header1_v1.4.h"

int main(){
  int isource,step,irec,i,idir;
  clock_t start=clock();
  setbuf(stdout,NULL);

/*** Initial Condition ***/
  checking_omp();
  read_shotrec();
  init_eh_field_3d();
  true_model();
  read_trawave();
  read_waveform();
  lattice_time_3d();
  init_cpml();
  media_coeff_3d();

/*** Check conficuration ***/
  time_check();
  fmax_check();

  start_fdtd();
/*** Calculating for each isource ***/
  for(isource=0;isource<shot_num;isource++){
    init_FILE(isource);

    printf("This is %2d loop for REC2REC\n",isource);
    printf("[PROGRESS  =>]  ");

/*** Set E.H field to zero ***/
    set_zero_eh();

/*** Starting calculation ***/
    for(step=0;step<it-1;step++){
      e_field_cpml42(EX,EY,EZ,HX,HY,HZ);
      epml_abcs4(EX,EY,EZ,HX,HY,HZ);
      read_source_3d(EX,signalX,isource,step);
      h_field_cpml42(EX,EY,EZ,HX,HY,HZ);
      hpml_abcs4(EX,EY,EZ,HX,HY,HZ);
      output(EX,EY,EZ,HX,HY,HZ,ofe1,ofe2,ofe3,ofh1,ofh2,ofh3,step);

      if(step%(it/20)==it/20-1) printf("#");
    }  // END of step
    close_FILE();
    laplaceToFreq(isource);

  }   // END of "isource"

/*** Computing Elapsed time ***/
  printf(" [ You succeed to computute REC2REC ]\n");
  printf("%f [s] \n",(double)(clock()-start)/CLOCKS_PER_SEC);
}
