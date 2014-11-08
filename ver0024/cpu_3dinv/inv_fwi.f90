!/////////////////////////////////////////////////////////////////////////////
! Full Waveform Inversion
!
! flaw chartの作成
!//////////////////////////////////////////////////////////////////////////


! program inv_fwi
!   use const_para
!   implicit none

!   write(*,*) "****************** Initial Condition **********************"
!   write(*,*) "****************** Backward calculation **********************"
!   write(*,*) "****************** Forward calculation **********************"
!   write(*,*) "****************** Start FWI **********************"
!   write(*,*) "******************  **********************"
!   write(*,*) "******************  **********************"
!   write(*,*) "******************  **********************"
!   write(*,*) "******************  **********************"

! end program inv_fwi


1.initicial Condition

  shot record位置の読み込み
  入力波形の設定
  グリッド作成
  境界条件の初期設定
  グラディエントの初期設定

2. 反復開始

  0クリア
  媒質の設定

  2.1 back propagation
    電磁場計算 a
      temporalに
        電磁場aの実部を配列EcalX_bに保存 b
        bを仮想領域から実領域へ c
        cからグリーン関数を求め、配列EcalX_bb,...に保存 d

  2.2 obserbed EMfieldの読み込み
        Eobs
        0クリア ehfield
        媒質の設定

  2.3 foward propagation
    電磁場計算 a
      sourceの読み込み
        電磁場aの実部を配列EcalX_bに保存 b
         bを仮想領域から実領域へ c
         cからグリーン関数を求め d
         dによるconvolutionに(G*J)をEcalX_ff,...に保存 e
         同時に、EcalX_ffLをもとめる f
    eの出力

  2.4 residual error の計算
    delE_f を求める delE_f = Eobs - EcalX_ff
    err_sum = err_sum + abs(delE_f)

  2.5 sensitivity の計算
    sensitiv = sensitiv + abs(EcalX_ff * EcalX_bb)
    sensitiv = sensitiv + abs(EcalY_ff * EcalY_bb)
    sensitiv = sensitiv + abs(EcalZ_ff * EcalZ_bb)

  2.6 delE を green fucntion と convolutoin する
    EcalX_back = EcalX_back + EcalX_bb * conjg(delEx_f)

  2.7 γ の計算
    grad = grad + real(EcalX_ff * EcalX_back
                     + EcalY_ff * EcalY_back
                     + EcalZ_ff * EcalZ_back) / sensitiv

  2.8 calc phi
    scaler = 0.01 * sigmax2 / gradmax
    sig2 = sig - scaler * grad

  2.9 媒質パラメータ2の作成

2回目
  2.10 foward propagatoin
    電磁場計算 a
    EcalX_b から EcalX_ff2をもとめ、

  2.11 calc alpha
    EcalX_ff2とEcalX_ffLを使用し、
    k1,k2を求める
    k1 = k1 + delEx_f * (EcalX_ff2 - EcalX_ffL)
    k1 = k1 + delEy_f * (EcalY_ff2 - EcalY_ffL)
    k1 = k1 + delEz_f * (EcalZ_ff2 - EcalZ_ffL)
    k2 = k2 + (EcalX_ff2 - EcalX_ffL)
    k2 = k2 + (EcalY_ff2 - EcalY_ffL)
    k2 = k2 + (EcalZ_ff2 - EcalZ_ffL)

  2.12 パラメータの更新
    alpha = 3 * abs(k1) / abs(k2)

    sig = sig + alpha * scaler * grad
    sig = minval (if sig < minval)
    sig = maxval (if sig >maxval)
    grad_p() = grad()


・・・・




! #include "fwi3d_cpu_header2_v1.4.h"

! int main(){
!   int i, j, iter,isource,idir,step,irec;
!   clock_t start = clock();
!   char fname[256];
!   setbuf(stdout,NULL);

! /*** Initial Condition ***/
!   init();

!   printf("Start calculation\n");
! /*** Calculating for iteratively ***/
!   for(iter=0;iter<MAXITR;iter++){
!     sprintf(fname,"./data/error%03d.dat",iter);
!     ofer = fopen(fname,"w");
!     printf("%03d th iterative\n",iter);
!     media_coeff_3d();
!     init_iterate();
!     init_FILE3(iter);

! //***********  Start backward calculation  **********//
!     for(irec=0;irec<rec_num;irec++){
!       for(idir=0;idir<NDIRECT;idir++){
!         banner1(idir,irec);
!         set_zero_eh();

! /*** Calculating for step ***/
!         for(step=0;step<it-1;step++){
!           backpropagation(idir, irec, step);
!           copytoEcal_b(EX,EY,EZ,step,idir, irec,ofe1,ofe2,ofe3);
!           if(step%(it/20)==it/20-1) printf("#");
!         }
!         printf("\n");
!         laplaceToFreq_back(irec,idir);
!         //output_backwave(irec,idir);
!       } // END loop for direction
!     }// END loop for ireceiver
!     close_FILE3();

! //***********  Start of FWD calculation  **********//
! /*** Calculating for each source ***/
!     for(isource=0;isource<shot_num;isource++){
!       printf("This is %2d loop for Transmitter \n",isource);

!       read_Eobs(isource); // Read observed value of Ex
!       //init_FILE2(isource, iter);
!       printf("[PROGRESS  =>]  ");
!       set_zero_eh();
!       media_coeff_3d();

! /*** Calculating for step ***/
!       for(step=0;step<it-1;step++){
!         fwdpropagation(isource, step);
!         copytoEcal(EX,EY,EZ,step);
!         if(step%(it/20)==it/20-1) printf("#");
!       } // END loop for step
!       laplaceToFreq_2(isource);
!       output_fwdwave(isource);
!       //close_FILE2();

!       // Calculation of residual error
!       residualE(iter);
!       sensitivity();
!       // Convolution of Greenfunc and conjugated delE
!       conv_GdelE();
!       // Calculation of gamma
!       calc_gamma();
!     } // END loop for isource

!     // Calclulation of delta model
!     calc_phi();
!     //printf("sig2  %e\n",sig[0]);
!     media_coeff_sig_tmp();

! //***********  again FWD calculation  **********//
!     for(isource=0;isource<shot_num;isource++){
!       printf("This is %2d loop for AGAIN Transmitter \n",isource);
!       printf("[PROGRESS  =>]  ");
!       set_zero_eh();
! /*** Calculating for step ***/
!       for(step=0;step<it-1;step++){
!         fwdpropagation(isource, step);
!         copytoEcal(EX,EY,EZ,step);
!         if(step%(it/20)==it/20-1) printf("#");
!       } // END loop for step
!       laplaceToFreq_3();
!       calc_alpha(isource);

!     } // END loop for isource
!     update_para2(iter);
!     //update_para3(iter);

!     show_error(iter,ofer);
!     fclose(ofer);
!   } // END loop for iter

!   printf("\n");
!   printf("%f [s] \n",(double)(clock()-start)/CLOCKS_PER_SEC);
! }
