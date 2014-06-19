////////////////////////////////////////////////////////////////////
//                                                                //
//        3D_FDTD_code_CPU            Ver 1.1                     //
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
#include "fwi3d_cpu_header2_v1.4.h"

int main(){
  int i, j, iter,isource,idir,step,irec;
  clock_t start = clock();
  setbuf(stdout,NULL);

/*** Initial Condition ***/
  init();

  printf("Start calculation\n");
/*** Calculating for iteratively ***/
  for(iter=0;iter<MAXITR;iter++){
    printf("%03d th iterative\n",iter);
    media_coeff_3d();
    init_iterate();
    init_FILE3(iter);

//***********  Start backward calculation  **********//
    for(irec=0;irec<rec_num;irec++){
      for(idir=0;idir<NDIRECT;idir++){
        banner1(idir,irec);
        set_zero_eh();

/*** Calculating for step ***/
        for(step=0;step<it-1;step++){
          backpropagation(idir, irec, step);
          copytoEcal_b(EX,EY,EZ,step,idir, irec,ofe1,ofe2,ofe3);
          if(step%(it/20)==it/20-1) printf("#");
        }
        printf("\n");
        laplaceToFreq_back(irec,idir);
        //output_backwave(irec,idir);
      } // END loop for direction
    }// END loop for ireceiver
    close_FILE3();

//***********  Start of FWD calculation  **********//
/*** Calculating for each source ***/
    for(idir=0;idir<3;idir++){
    for(isource=0;isource<shot_num;isource++){
      printf("This is %2d loop for Transmitter \n",isource);

      read_3C_Eobs(isource,idir); // Read observed value of Ex
      //init_FILE2(isource, iter);
      printf("[PROGRESS  =>]  ");
      set_zero_eh();
      media_coeff_3d();

/*** Calculating for step ***/
      for(step=0;step<it-1;step++){
        ThreeCfwdpropagation(isource, step,idir);
        copytoEcal(EX,EY,EZ,step);
        if(step%(it/20)==it/20-1) printf("#");
      } // END loop for step
      laplaceToFreq_2(isource);
      output_fwdwave(isource);
      //close_FILE2();

      // Calculation of residual error
      residualE(iter);
      sensitivity();
      // Convolution of Greenfunc and conjugated delE
      conv_GdelE();
      // Calculation of gamma
      calc_gamma();
    } // END loop for isource

    // Calclulation of delta model
    calc_phi();
    //printf("sig2  %e\n",sig[0]);
    media_coeff_sig_tmp();

//***********  again FWD calculation  **********//
    for(isource=0;isource<shot_num;isource++){
      printf("This is %2d loop for AGAIN Transmitter \n",isource);
      printf("[PROGRESS  =>]  ");
      set_zero_eh();
/*** Calculating for step ***/
      for(step=0;step<it-1;step++){
        ThreeCfwdpropagation(isource, step,idir);
        copytoEcal(EX,EY,EZ,step);
        if(step%(it/20)==it/20-1) printf("#");
      } // END loop for step
      laplaceToFreq_3();
      calc_alpha(isource);

    } // END loop for isource
    update_para2(iter);

    show_error(iter);
      }
  } // END loop for iter
  printf("\n");
  printf("%f [s] \n",(double)(clock()-start)/CLOCKS_PER_SEC);
  fclose(ofer);
}
