/////////////////////////////////////////////////////////////
//                                                         //
//        3D_FDTD_code_CPU_MPI        Ver 1.0              //
//                                                         //
//      * Coded by         N. Imamura  (2012/07/28)        //
//                                                         //
//     v1.0 : basic FDTD is completed  (2012/08/08)        //
//     v1.1 : fictitious FDTD is completed  (2012/09/14)   //
//                                                         //
/////////////////////////////////////////////////////////////
#include "fwi3d_cpu_header1_v1.1.h"

int main(){
  int i,j,k,ijk,isource,step;
  int jx,jy,jz;
  clock_t start=clock();
  setbuf(stdout,NULL);

/*** Initial Condition ***/
  checking_omp();
  read_shotrec();
  lattice_time_3d();
  read_trawave();
  init_eh_field_3d();
  init_cpml();
  media_coeff_3d();
//  init_eh_mur_3d();

/*** Calculating for each isource ***/
  for(isource=0;isource<shot_num;isource++){
    init_FILE(isource);

    printf("This is %2d loop for REC2REC\n",isource);
    printf("[PROGRESS  =>]  ");
// set eh fields to zero
    set_zero();

/*** Starting calculation ***/
    for(step=0;step<it-1;step++){
      //e_field_cpml(EX,EY,EZ,HX,HY,HZ);
      e_field_cpml42(EX,EY,EZ,HX,HY,HZ);
      epml_abcs4(EX,EY,EZ,HX,HY,HZ);
//      e_mur2(EX,EY,EZ,HX,HY,HZ);
//      cerjan_e(EX,EY,EZ);
      read_source_3d(EX,signalX,isource,step);
      //h_field_cpml(EX,EY,EZ,HX,HY,HZ);
      h_field_cpml42(EX,EY,EZ,HX,HY,HZ);
      hpml_abcs4(EX,EY,EZ,HX,HY,HZ);
//      cerjan_h(HX,HY,HZ);
      output(EX,EY,EZ,HX,HY,HZ,ofe1,ofe2,ofe3,ofh1,ofh2,ofh3,step);

      if(step%(it/20)==it/20-1) printf("#");
    }  // END of step
    close_FILE();
  }   // END of "isource"

/*** Computing Elapsed time ***/
  printf("\n");
  printf("%f [s] \n",(double)(clock()-start)/CLOCKS_PER_SEC);
}
