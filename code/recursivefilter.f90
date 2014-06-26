! **********************************************************
! コンパイルの方法: f77 -o wavedata wavedata.f
! システムによっては，-lmオプションが必要な場合がある．

! GMTで波形を書かせるためのデータを作る
! RGES ascii format => GMT format (2 colunm ascii)
! chebyshev filter (bandpass) を書けることが出来る．

! ○　波形数値データについて

! 　このディレクトリ REC の下には、観測年ごとに RECxx の形のサブディレ
! クトリがあり、この中に爆破点別に RGxxS1 などの名前で波形データが入っ
! ています。各観測点ごとのフォーマットは次のとおりで、これが１爆破点分
! 切れ目なく入っています。ファイルの最後にはブランク行が入っています。

! ヘッダ　FORMAT(A8, I8, F8.0, F8.0)
! 　　　　　　　観測点名など、データ数、距離、１点目のデータの走時
! データ　FORMAT(14I5)
! 　　　　　　　100Hzサンプリングで正が地動ＵＰ

! また、１点目のデータの走時の後ろに、最大振幅が mgal 単位で入っている
! ものがあります。

! ○  注意
! １．このデータを使った研究等の発表にあたっては、データは爆破地震動研
!  究グループによることを明記し、関係する観測報告等を引用してください。
! ２．AD変換は、1981年以降のデータについて、IF-800、PC9801およびhp350H
!  を用いて、東京大学地震研究所爆破地震学部門で行なわれました。
! ３．1981年の坂出・一宇爆破のデータは、1993年11月に岐阜大学教育学部４
!  年（当時）の臼井友美さんが地震研究所でAD変換したものを、提供してい
!  ただきました。

! input parameter (from STDIN)
!  filename of wavedata(RGES format)
!   shot number, normalize factor, resion of time, resion of disance,
!   reduce velosity, filter flag, distance flag
!  low freq., high freq., cut freq.

! 各ショットの測線の端からの距離は，data文でshotという配列に入力してある．
! **********************************************************

      real start,dist,amp,x1,x2,t1,t2,fl,fh,fs,dmax,reduce
      real shot(6),h(200),data(10000)
      integer ishot,ffrag,aflag,dflag,ndata
      character input*50,obs*8

***** shot data *****
      data shot/32.7,125.2,227.5,0.0,86.6,181.8/
      read(*,'(a)') input
      read(*,*) ishot,amp,aflag,x1,x2,t1,t2,reduce,ffrag,dflag
      read(*,*) fl,fh,fs
      fl=fl/100
      fh=fh/100
      fs=fs/100
      open(10,file=input,status='old',access='sequential')
      do i=1,300
      read(10,1000,end=999)obs,ndata,dist,start
      read(10,2000,end=999)(data(j),j=1,ndata)
***** filtering (chebishev filter) *****
        if(ffrag.eq.1)then
          call chbpas(h,m,gn,n,eps,fl,fh,fs,0.5,5.0)
          call tandem(data,data,ndata,h,m,1)
          do j=1,ndata
            data(j)=data(j)*gn
          end do
        end if
***** amplitude manupilation *****
        if(aflag.eq.1)then
!          real amplitude **
          dmax=8000.0
        else if(aflag.eq.2)then
!          normalize by max amplitude **
          dmax=abs(data(1))
          do j=2,ndata
            dmax=max(dmax,abs(data(j)))
          end do
        end if
***** output wave data for GMT (piping to "psxy") *****
        do j=1,ndata
          if (dflag.eq.1)then
            x=dist+shot(ishot)+amp*data(j)/dmax
          else if (dflag.eq.2)then
            x=dist+amp*data(j)/dmax
          end if
          t=start+0.01*(j-1)-abs(dist)/reduce
          if((t1.le.t).and.(t.le.t2))then
            write(*,3000) x,t
          end if
        end do
        write(*,4000)
      end do
  999 close(10)
 1000 format(a8,i8,f8.2,f8.2)
 2000 format(14f5)
 3000 format(f8.4,f8.2)
 4000 format('>')
      stop
      end


      SUBROUTINE  CHBPAS(H,M,GN,N,EPS,FL,FH,FS,AP,AS)
        implicit none
      COMPLEX  R(2),OJ,CQ
      DIMENSION  H(4)
      DATA  PI/3.141593/,HP/1.570796/
!
!       CHEBYSHEV BAND-PASS FILTER COEFFICIENT
!
!       ARGUMENTS
!         H      : FILTER COEFFICIENT
!         M      : ORDER OF FILTER
!         GN     : GAIN FACTOR
!         N      : DEGREE OF CHEBYSHEV POLYNOMIAL
!         EPS    : RIPPLE
!         FL     : LOW FREQUENCY CUT-OFF  (NON-DIMENSIONAL)
!         FH     : HIGH FREQUENCY CUT-OFF
!         FS     : STOP BAND FREQUENCY
!         AP     : MAX. ATTENUATION IN PASS BAND
!         AS     : MIN. ATTENUATION IN STOP BAND
!
!       M. SAITO  (16/XII/75)
!
      WL = AMIN1(ABS(FL),ABS(FH))*PI
      WH = AMAX1(ABS(FL),ABS(FH))*PI
      WS = ABS(FS)*PI
      IF( WL.EQ.WH .OR. WL.EQ.0. .OR. WH.GE.HP .OR. WS.EQ.0. .OR.
     *    WS.GE.HP .OR. (WS-WL)*(WS-WH).LE.0. )  GO TO  100
!       DETERMINE N, EPS & C
      CLH= COS(WL)*COS(WH)
      C  = CLH/SIN(WH-WL)
      CC = C*C
      WW = TAN(WL)*TAN(WH)*CC
      WS = TAN(WS)*C
      OS = ABS(WS-WW/WS)
      TS = ALOG(OS+SQRT((OS-1.)*(OS+1.)))
      PA = AMIN1(ABS(AP),ABS(AS))
      SA = AMAX1(ABS(AP),ABS(AS))
      IF( PA.EQ.0. )  PA = 0.1
      IF( SA.EQ.0. )  SA = 5.
      PS = PA/SA
      N  = MAX0(2,IFIX(ABS(ALOG((1.+SQRT((1.-PS)*(1.+PS)))/PS)/TS)+0.5))
      FN = N
      EPS= SQRT(PA*SA/COSH(TS*FN))
!
      K  = N/2
      M  = K*2
      L  = 0
      FJ = 1.
      G  = 1./(EPS*FLOAT(2**(N-1)))
      PS = ALOG((1.+SQRT(1.+EPS**2))/EPS)/FN
      CH = COSH(PS)
      SH = SINH(PS)
      DP = HP/FN
!
      DO  2  J=1,K
        OJ = CMPLX(COS(DP*FJ)*CH,SIN(DP*FJ)*SH)*0.5
        FJ = FJ+2.
        CQ = CSQRT(OJ**2+WW)
        R(1) = OJ+CQ
        R(2) = OJ-CQ
        G  = G*CC
!
        DO  1  I=1,2
          RE =  REAL(R(I))**2
          RI = AIMAG(R(I))
          A  = 1./((C+RI)**2+RE)
          G  = G*A
            H(L+1) = 0.
            H(L+2) =-1.
            H(L+3) = 2.*((RI-C)*(RI+C)+RE)*A
            H(L+4) = ((RI-C)**2+RE)*A
            L = L+4
    1   CONTINUE
!
    2 CONTINUE
!       EXIT
      GN = G
      IF( N.EQ.M )  RETURN
!       FOR ODD N
      M  = M+1
      WPC= CC*COS(WH-WL)/CLH
      WMC=-CC*COS(WL+WH)/CLH
          A  = 1./(WPC+C*SH)
          GN = G*C*A
            H(L+1) = 0.
            H(L+2) =-1.
            H(L+3) = 2.*WMC*A
            H(L+4) = (WPC-C*SH)*A
      RETURN
!       ERROR
  100 WRITE(6,101)  FL,FH,FS
  101  FORMAT(/1X,5('?'),'   (CHBPAS)   INVALID INPUT   FL =',
     *        1PE14.6,3X,'FH =',E14.6,3X,'FS =',E14.6,3X,5('?')//)
      RETURN
      END


      SUBROUTINE  TANDEM(X,Y,N,H,M,NML)
      DIMENSION  X(N),Y(N),H(4)
!
!       RECURSIVE FILTERING IN SERIES
!
!       ARGUMENTS
!         X      : INPUT TIME SERIES
!         Y      : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
!         N      : LENGTH OF X & Y
!         H      : COEFFICIENTS OF FILTER
!         M      : ORDER OF FILTER
!         NML    : >0 ; FOR NORMAL  DIRECTION FILTERING
!                  <0 ;     REVERSE DIRECTION FILTERING
!
!       SUBROUTINE REQUIRED : RECFIL
!
!       M. SAITO  (6/XII/75)
!
      IF( N.LE.0 .OR. M.LE.0 )  GO TO  2
!  1-ST CALL
      CALL  RECFIL(X,Y,N,H,NML)
      IF( M.LE.1 )  RETURN
!****  2-ND AND AFTER
      DO  1  I=2,M
        CALL  RECFIL(Y,Y,N,H(I*4-3),NML)
    1 CONTINUE
!
      RETURN
!  ERROR
    2 WRITE(6,3)  N,M
    3  FORMAT(//1X,5('?'),3X,'(TANDEM)',3X,'INVALID INPUT',3X,'N =',I5,
     *        3X,'M =',I5,5('?')//)
      RETURN
      END


      SUBROUTINE  RECFIL(X,Y,N,H,NML)
      DIMENSION  X(N),Y(N),H(4)
!
!       RECURSIVE FILTERING : F(Z) = (1+A*Z+AA*Z**2)/(1+B*Z+BB*Z**2)
!
!       ARGUMENTS
!         X      : INPUT TIME SERIES
!         Y      : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
!         N      : LENGTH OF X & Y
!         H      : FILTER COEFFICIENTS ; H(1)=A, H(2)=AA, H(3)=B, H(4)=BB
!         NML    : >0 ; FOR NORMAL  DIRECTION FILTERING
!                  <0 ; FOR REVERSE DIRECTION FILTERING
!
!       M. SAITO  (6/XII/75)
!
      IF(   N.LE.0 )  GO TO  4
      IF( NML.GE.0 )  GO TO  1
!  REVERSE FILTERING
      J  = N
      JD =-1
        GO TO  2
!  NORMAL FILTERING
    1 J  = 1
      JD = 1
!
    2 A  = H(1)
      AA = H(2)
      B  = H(3)
      BB = H(4)
      U1 = 0.
      U2 = 0.
      V1 = 0.
      V2 = 0.
!  FILTERING
      DO  3  I=1,N
        U3 = U2
        U2 = U1
        U1 = X(J)
        V3 = V2
        V2 = V1
        V1 = U1+A*U2+AA*U3-B*V2-BB*V3
        Y(J) = V1
        J  = J+JD
    3 CONTINUE
!  EXIT
      RETURN
!  ERROR
    4 WRITE(6,5)  N
    5  FORMAT(//1X,5('?'),3X,'(RECFIL)',3X,'INVALID INPUT',3X,'N =',I5,
     *        3X,5('?')//)
      RETURN
      END