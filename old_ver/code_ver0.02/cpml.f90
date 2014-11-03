!///////////////////////////////////////////////////////////////////////////////
! CPML係数設定
! optToMaxをどうする？
!///////////////////////////////////////////////////////////////////////////////
subroutine init_cpml
  use const_para
  implicit none

        !!!    sig_max = -(m+1)*lnR0 / (2.0d0*(sqrt(myu/epsi))*ncpml*dx)  !ln(R(0));反射係数!!!
        !    sig_opt = (dble(m)+1.0d0) / (150.0d0*pi*sqrt(epsir)*dx)
        !    sig_max = 0.7d0*sig_opt

  !係数の設定x
open(97,file='fxcpml.d')
      write(97,"(13a)")  " # i, esig_x(i), msig_x(i), ekappa_x(i), mkappa_x(i), ae_x(i),   &
       am_x(i),    be_x(i),    bh_x(i),    ce_x(i),    ch_x(i),    kedx(i),    khdx(i)"
  delta = ncpml*dx
  sig_max = (nn+order+1.0d0)*cmax*log(1.0d0/Rcoef) / (2.0d0*delta) * optToMax  !!x方向だけoptToMaxかけるの？

! do i = 1,nx
!     if(i<=ncpml) then
do i=1,ncpml
      esig_x(i)  = sig_max * ((dble(ncpml)-dble(i)      )/(dble(ncpml)-1.0d0))**dble(nn+order)
      msig_x(i)  = sig_max * ((dble(ncpml)-dble(i)-0.5d0)/(dble(ncpml)-1.0d0))**dble(nn+order)  !!!-i-1/2の取り扱い
      ekappa_x(i)= 1.0d0 + (kappa_max-1.0d0) * ((dble(ncpml)-dble(i)      )/(dble(ncpml)-1.0d0))**dble(nn)
      mkappa_x(i)= 1.0d0 + (kappa_max-1.0d0) * ((dble(ncpml)-dble(i)-0.5d0)/(dble(ncpml)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
      ae_x(i)    = a_max * ((dble(i)      )/(dble(ncpml)-1.0d0))**dble(ma)
      am_x(i)    = a_max * ((dble(i)+0.5d0)/(dble(ncpml)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

      be_x(i)    = exp(-(esig_x(i)/ekappa_x(i)+ae_x(i)) * dt) !/epsi0)
      bh_x(i)    = exp(-(msig_x(i)/mkappa_x(i)+am_x(i)) * dt) !/epsi0)
      ce_x(i)    = esig_x(i)*(be_x(i)-1.0d0) / (esig_x(i) + ekappa_x(i)*ae_x(i)) / ekappa_x(i)
      ch_x(i)    = msig_x(i)*(bh_x(i)-1.0d0) / (msig_x(i) + mkappa_x(i)*am_x(i)) / mkappa_x(i)
      kedx(i)    = ekappa_x(i)*dx
      khdx(i)    = mkappa_x(i)*dx !!!(i-1/2)dxの取り扱い
enddo
!     else if(i>=nx-ncpml+1) then
do i=nx-ncpml+1,nx
      esig_x(i)  = sig_max * ((dble(i)-dble(nx)+1.0d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn+order)
      msig_x(i)  = sig_max * ((dble(i)-dble(nx)+0.5d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn+order)  !!!-i-1/2の取り扱い
      ekappa_x(i)= 1.0d0 + (kappa_max-1.0d0) * ((dble(i)-dble(nx)+1.0d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn)
      mkappa_x(i)= 1.0d0 + (kappa_max-1.0d0) * ((dble(i)-dble(nx)+0.5d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
      ae_x(i)    = a_max * ((dble(-i)+dble(nx)-1.0d0)/(dble(ncpml)-1.0d0))**dble(ma)
      am_x(i)    = a_max * ((dble(-i)+dble(nx)-0.5d0)/(dble(ncpml)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

      be_x(i)    = exp(-(esig_x(i)/ekappa_x(i)+ae_x(i)) * dt) !/epsi0)
      bh_x(i)    = exp(-(msig_x(i)/mkappa_x(i)+am_x(i)) * dt) !/epsi0)
      ce_x(i)    = esig_x(i)*(be_x(i)-1.0d0) / (esig_x(i) + ekappa_x(i)*ae_x(i)) / ekappa_x(i)
      ch_x(i)    = msig_x(i)*(bh_x(i)-1.0d0) / (msig_x(i) + mkappa_x(i)*am_x(i)) / mkappa_x(i)
      kedx(i)    = ekappa_x(i)*dx
      khdx(i)    = mkappa_x(i)*dx !!!(i-1/2)dxの取り扱い
enddo
!     else
do i=ncpml+1,nx-ncpml
      esig_x(i)  = 0.0d0
      msig_x(i)  = 0.0d0
      ekappa_x(i)= 1.0d0
      mkappa_x(i)= 1.0d0
      ae_x(i)    = 0.0d0
      am_x(i)    = 0.0d0
      be_x(i)    = 0.0d0
      bh_x(i)    = 0.0d0
      ce_x(i)    = 0.0d0
      ch_x(i)    = 0.0d0
      kedx(i)    = ekappa_x(i)*dx
      khdx(i)    = mkappa_x(i)*dx
!     endif
enddo
do i=1,nx
    write(97,"(I3,12e12.4)")  i,esig_x(i),msig_x(i),ekappa_x(i),mkappa_x(i),ae_x(i),am_x(i),be_x(i),bh_x(i),ce_x(i),ch_x(i), kedx(i),khdx(i)
enddo
! enddo
close(97)


!係数の設定y
open(98,file='fycpml.d')
      write(98,"(13a)")  " # j, esig_y(j), msig_y(j), ekappa_y(j), mkappa_y(j), ae_y(j),   &
       am_y(j),    be_y(j),    bh_y(j),    ce_y(j),    ch_y(j),    kedy(j),    khdy(j)"
  delta = ncpml*dy
  sig_max = (nn+order+1.0d0)*cmax*log(1.0d0/Rcoef) / (2.0d0*delta)

do j=1,ny
if(j<=ncpml) then
      esig_y(j)  = sig_max * ((dble(ncpml)-dble(j)      )/(dble(ncpml)-1.0d0))**dble(nn+order)
      msig_y(j)  = sig_max * ((dble(ncpml)-dble(j)-0.5d0)/(dble(ncpml)-1.0d0))**dble(nn+order)  !!!-i-1/2の取り扱い
      ekappa_y(j)= 1.0d0 + (kappa_max-1.0d0) * ((dble(ncpml)-dble(j)      )/(dble(ncpml)-1.0d0))**dble(nn)
      mkappa_y(j)= 1.0d0 + (kappa_max-1.0d0) * ((dble(ncpml)-dble(j)-0.5d0)/(dble(ncpml)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
      ae_y(j)    = a_max * ((dble(j)      )/(dble(ncpml)-1.0d0))**dble(ma)
      am_y(j)    = a_max * ((dble(j)+0.5d0)/(dble(ncpml)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

      be_y(j)    = exp(-(esig_y(j)/ekappa_y(j)+ae_y(j)) * dt) !/epsi0)
      bh_y(j)    = exp(-(msig_y(j)/mkappa_y(j)+am_y(j)) * dt) !/epsi0)
      ce_y(j)    = esig_y(j)*(be_y(j)-1.0d0) / (esig_y(j) + ekappa_y(j)*ae_y(j)) / ekappa_y(j)
      ch_y(j)    = msig_y(j)*(bh_y(j)-1.0d0) / (msig_y(j) + mkappa_y(j)*am_y(j)) / mkappa_y(j)
      kedy(j)    = ekappa_y(j)*dy
      khdy(j)    = mkappa_y(j)*dy !!!(i-1/2)dyの取り扱い

    else if(j>=ny-ncpml+1) then
      esig_y(j)  = sig_max * ((dble(j)-dble(ny)+1.0d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn+order)
      msig_y(j)  = sig_max * ((dble(j)-dble(ny)+0.5d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn+order)  !!!-i-1/2の取り扱い
      ekappa_y(j)= 1.0d0 + (kappa_max-1.0d0) * ((dble(j)-dble(ny)+1.0d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn)
      mkappa_y(j)= 1.0d0 + (kappa_max-1.0d0) * ((dble(j)-dble(ny)+0.5d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
      ae_y(j)    = a_max * ((dble(-j)+dble(ny)-1.0d0)/(dble(ncpml)-1.0d0))**dble(ma)
      am_y(j)    = a_max * ((dble(-j)+dble(ny)-0.5d0)/(dble(ncpml)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

      be_y(j)    = exp(-(esig_y(j)/ekappa_y(j)+ae_y(j)) * dt) !/epsi0)
      bh_y(j)    = exp(-(msig_y(j)/mkappa_y(j)+am_y(j)) * dt) !/epsi0)
      ce_y(j)    = esig_y(j)*(be_y(j)-1.0d0) / (esig_y(j) + ekappa_y(j)*ae_y(j)) / ekappa_y(j)
      ch_y(j)    = msig_y(j)*(bh_y(j)-1.0d0) / (msig_y(j) + mkappa_y(j)*am_y(j)) / mkappa_y(j)
      kedy(j)    = ekappa_y(j)*dy
      khdy(j)    = mkappa_y(j)*dy !!!(i-1/2)dyの取り扱い

    else
      esig_y(j)  = 0.0d0
      msig_y(j)  = 0.0d0
      ekappa_y(j)= 1.0d0
      mkappa_y(j)= 1.0d0
      ae_y(j)    = 0.0d0
      am_y(j)    = 0.0d0
      be_y(j)    = 0.0d0
      bh_y(j)    = 0.0d0
      ce_y(j)    = 0.0d0
      ch_y(j)    = 0.0d0
      kedy(j)    = ekappa_y(j)*dy
      khdy(j)    = mkappa_y(j)*dy
    endif
    write(98,"(I3,12e12.4)")  j,esig_y(j),msig_y(j),ekappa_y(j),mkappa_y(j),ae_y(j),am_y(j),be_y(j),bh_y(j), ce_y(j),ch_y(j), kedy(j),khdy(j)
enddo
close(98)


!係数の設定z
open(99,file='fzcpml.d')
      write(97,"(13a)")  " # k, esig_z(k), msig_z(k), ekappa_z(k), mkappa_z(k), ae_z(k),   &
       am_z(k),    be_z(k),    bh_z(k),    ce_z(k),    ch_z(k),    kedz(k),    khdz(k)"
  delta = ncpml*dz
  sig_max = (nn+order+1.0d0)*cmax*log(1.0d0/Rcoef) / (2.0d0*delta)

do k=1,nz
if(k<=ncpml) then
      esig_z(k)  = sig_max * ((dble(ncpml)-dble(k)      )/(dble(ncpml)-1.0d0))**dble(nn+order)
      msig_z(k)  = sig_max * ((dble(ncpml)-dble(k)-0.5d0)/(dble(ncpml)-1.0d0))**dble(nn+order)  !!!-i-1/2の取り扱い
      ekappa_z(k)= 1.0d0 + (kappa_max-1.0d0) * ((dble(ncpml)-dble(k)      )/(dble(ncpml)-1.0d0))**dble(nn)
      mkappa_z(k)= 1.0d0 + (kappa_max-1.0d0) * ((dble(ncpml)-dble(k)-0.5d0)/(dble(ncpml)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
      ae_z(k)    = a_max * ((dble(k)      )/(dble(ncpml)-1.0d0))**dble(ma)
      am_z(k)    = a_max * ((dble(k)+0.5d0)/(dble(ncpml)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

      be_z(k)    = exp(-(esig_z(k)/ekappa_z(k)+ae_z(k)) * dt) !/epsi0)
      bh_z(k)    = exp(-(msig_z(k)/mkappa_z(k)+am_z(k)) * dt) !/epsi0)
      ce_z(k)    = esig_z(k)*(be_z(k)-1.0d0) / (esig_z(k) + ekappa_z(k)*ae_z(k)) / ekappa_z(k)
      ch_z(k)    = msig_z(k)*(bh_z(k)-1.0d0) / (msig_z(k) + mkappa_z(k)*am_z(k)) / mkappa_z(k)
      kedz(k)    = ekappa_z(k)*dz
      khdz(k)    = mkappa_z(k)*dz !!!(i-1/2)dzの取り扱い

    else if(k>=nz-ncpml+1) then
      esig_z(k)  = sig_max * ((dble(k)-dble(nz)+1.0d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn+order)
      msig_z(k)  = sig_max * ((dble(k)-dble(nz)+0.5d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn+order)  !!!-i-1/2の取り扱い
      ekappa_z(k)= 1.0d0 + (kappa_max-1.0d0) * ((dble(k)-dble(nz)+1.0d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn)
      mkappa_z(k)= 1.0d0 + (kappa_max-1.0d0) * ((dble(k)-dble(nz)+0.5d0+dble(ncpml))/(dble(ncpml)-1.0d0))**dble(nn)  !!!-i-1/2の取扱い
      ae_z(k)    = a_max * ((dble(-k)+dble(nz)-1.0d0)/(dble(ncpml)-1.0d0))**dble(ma)
      am_z(k)    = a_max * ((dble(-k)+dble(nz)-0.5d0)/(dble(ncpml)-1.0d0))**dble(ma) !!!-i-1/2の取り扱い

      be_z(k)    = exp(-(esig_z(k)/ekappa_z(k)+ae_z(k)) * dt) !/epsi0)
      bh_z(k)    = exp(-(msig_z(k)/mkappa_z(k)+am_z(k)) * dt) !/epsi0)
      ce_z(k)    = esig_z(k)*(be_z(k)-1.0d0) / (esig_z(k) + ekappa_z(k)*ae_z(k)) / ekappa_z(k)
      ch_z(k)    = msig_z(k)*(bh_z(k)-1.0d0) / (msig_z(k) + mkappa_z(k)*am_z(k)) / mkappa_z(k)
      kedz(k)    = ekappa_z(k)*dz
      khdz(k)    = mkappa_z(k)*dz !!!(i-1/2)dzの取り扱い

    else
      esig_z(k)  = 0.0d0
      msig_z(k)  = 0.0d0
      ekappa_z(k)= 1.0d0
      mkappa_z(k)= 1.0d0
      ae_z(k)    = 0.0d0
      am_z(k)    = 0.0d0
      be_z(k)    = 0.0d0
      bh_z(k)    = 0.0d0
      ce_z(k)    = 0.0d0
      ch_z(k)    = 0.0d0
      kedz(k)    = ekappa_z(k)*dz
      khdz(k)    = mkappa_z(k)*dz
    endif
    write(99,"(I3,12e12.4)")  k,esig_z(k),msig_z(k),ekappa_z(k),mkappa_z(k),ae_z(k),am_z(k),be_z(k),bh_z(k), ce_z(k),ch_z(k), kedz(k),khdz(k)
enddo
close(99)
end subroutine init_cpml


!!!Convolutional PML_E !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sigmax amax kappamax の求め方
!
!
!
!kappa >=1 real
!sigi>0 real
!ai=alphai>0 real
!psi loop あと３枚必要？
!演算中の倍精度の表現に注意.d0,
!be_x(5)が0になってしまう問題.１になるはず,epsi0=8.854d-12
!epsi0が小さすぎる!!:q
!cb_x,cb_y,cb_zを仮想領域にあわせる
!係数が1-ncpml, nx-(nx-nxpml)でことなるので調整
!epsi0=1 from imamu
!esig_xの値をどう取るか
!sig(i,j,k) or esig_x(i) ??
!operator half length ln = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CPML_E4(Ex,Ey,Ez,Hx,Hy,Hz)!,sig)!,cmax)
    use const_para
    implicit none

    complex(kind(0d0)), intent(inout) ::Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(in) ::   Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)

!psi-update

    !xe-PML4 loop(-)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 3,ncpml+1
                psi_Ezx1(i,j,k) = be_x(i) * psi_Ezx1(i,j,k) &
                                + ce_x(i) * (c1*Hy(i,j,k)-c1*Hy(i-1,j,k) + c2*Hy(i+1,j,k)-c2*Hy(i-2,j,k)) / dx
                psi_Eyx1(i,j,k) = be_x(i) * psi_Eyx1(i,j,k) &
                                + ce_x(i) * (c1*Hz(i,j,k)-c1*Hz(i-1,j,k) + c2*Hz(i+1,j,k)-c2*Hz(i-2,j,k)) / dx
                Ez(i,j,k) = Ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
                Ey(i,j,k) = Ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
            enddo
        enddo
    enddo
     !xe-PML4 loop(+)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = nx-ncpml+1,nx-1
                psi_Ezx1(i,j,k) = be_x(i) * psi_Ezx1(i,j,k) &
                                + ce_x(i) * (c1*Hy(i,j,k)-c1*Hy(i-1,j,k) + c2*Hy(i+1,j,k)-c2*Hy(i-2,j,k)) / dx
                psi_Eyx1(i,j,k) = be_x(i) * psi_Eyx1(i,j,k) &
                                + ce_x(i) * (c1*Hz(i,j,k)-c1*Hz(i-1,j,k) + c2*Hz(i+1,j,k)-c2*Hz(i-2,j,k)) / dx
                Ez(i,j,k) = Ez(i,j,k) + cb_z(i,j,k)*psi_Ezx1(i,j,k)
                Ey(i,j,k) = Ey(i,j,k) - cb_y(i,j,k)*psi_Eyx1(i,j,k)
            enddo
        enddo
    enddo

    !ye-PML4 loop(-)
    do k = 1,nz-1
        do j = 3,ncpml+1
            do i = 1,nx-1
                psi_Exy1(i,j,k) = be_y(j) * psi_Exy1(i,j,k) &
                                + ce_y(j) * (c1*Hz(i,j,k)-c1*Hz(i,j-1,k) + c2*Hz(i,j+1,k)-c2*Hz(i,j-2,k)) / dy
                psi_Ezy1(i,j,k) = be_y(j) * psi_Ezy1(i,j,k) &
                                + ce_y(j) * (c1*Hx(i,j,k)-c1*Hx(i,j-1,k) + c2*Hx(i,j+1,k)-c2*Hx(i,j-2,k)) / dy
                Ex(i,j,k) = Ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
                Ez(i,j,k) = Ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
            enddo
        enddo
    enddo
   !ye-PML4 loop(+)
    do k = 1,nz-1
        do j = ny-ncpml+1,ny-1
            do i = 1,nx-1
                psi_Exy1(i,j,k) = be_y(j) * psi_Exy1(i,j,k) &
                                + ce_y(j) * (c1*Hz(i,j,k)-c1*Hz(i,j-1,k) + c2*Hz(i,j+1,k)-c2*Hz(i,j-2,k)) / dy
                psi_Ezy1(i,j,k) = be_y(j) * psi_Ezy1(i,j,k) &
                                + ce_y(j) * (c1*Hx(i,j,k)-c1*Hx(i,j-1,k) + c2*Hx(i,j+1,k)-c2*Hx(i,j-2,k)) / dy
                Ex(i,j,k) = Ex(i,j,k) + cb_x(i,j,k)*psi_Exy1(i,j,k)
                Ez(i,j,k) = Ez(i,j,k) - cb_z(i,j,k)*psi_Ezy1(i,j,k)
            enddo
        enddo
    enddo

    !ze-PML4 loop(-)
    do k = 3,ncpml+1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Eyz1(i,j,k) = be_z(k) * psi_Eyz1(i,j,k) &
                                + ce_z(k) * (c1*Hx(i,j,k)-c1*Hx(i,j,k-1) + c2*Hx(i,j,k+1)-c2*Hx(i,j,k-2)) / dz
                psi_Exz1(i,j,k) = be_z(k) * psi_Exz1(i,j,k) &
                                + ce_z(k) * (c1*Hy(i,j,k)-c1*Hy(i,j,k-1) + c2*Hy(i,j,k+1)-c2*Hy(i,j,k-2)) / dz
                Ey(i,j,k) = Ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
                Ex(i,j,k) = Ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
            enddo
        enddo
    enddo
    !ze-PML4 loop(+)
    do k = nz-ncpml+1,nz-1
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Eyz1(i,j,k) = be_z(k) * psi_Eyz1(i,j,k) &
                                + ce_z(k) * (c1*Hx(i,j,k)-c1*Hx(i,j,k-1) + c2*Hx(i,j,k+1)-c2*Hx(i,j,k-2)) / dz
                psi_Exz1(i,j,k) = be_z(k) * psi_Exz1(i,j,k) &
                                + ce_z(k) * (c1*Hy(i,j,k)-c1*Hy(i,j,k-1) + c2*Hy(i,j,k+1)-c2*Hy(i,j,k-2)) / dz
                Ey(i,j,k) = Ey(i,j,k) + cb_y(i,j,k)*psi_Eyz1(i,j,k)
                Ex(i,j,k) = Ex(i,j,k) - cb_x(i,j,k)*psi_Exz1(i,j,k)
            enddo
        enddo
    enddo

end subroutine CPML_E4





!Convolutional PML_H !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sigmax amax kappamax の求め方 sig*の求め方
!
!
!
!kappa >=1 real
!sigi>0 real
!ai=alphai>0 real
!psi loop あと３枚必要？
!割り算のとき倍精度d0に注意!!
!db_x,db_y,db_zを仮想領域にあわせる
!epsi0=1 from imamu
!sig* の求め方
!db_z[ijk] = dt/MU0 /(1.f+(msigz[k]*dt)/(2.f*eps2));
!operator half length ln = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CPML_H4(Ex,Ey,Ez,Hx,Hy,Hz)!,sig,myu)!,cmax)
    use const_para
    implicit none

    complex(kind(0d0)), intent(in) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)

!psi-update

!xh-PML4 loop(-)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 2,ncpml
                psi_Hzx1(i,j,k) = bh_x(i) * psi_Hzx1(i,j,k)&
                                + ch_x(i) * ( c1*Ey(i+1,j,k)-c1*Ey(i,j,k) + c2*Ey(i+2,j,k)-c2*Ey(i-1,j,k) ) / dx
                psi_Hyx1(i,j,k) = bh_x(i) * psi_Hyx1(i,j,k)&
                                + ch_x(i) * ( c1*Ez(i+1,j,k)-c2*Ez(i,j,k) + c2*Ez(i+2,j,k)-c2*Ez(i-1,j,k) ) / dx
                Hz(i,j,k) = Hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
                Hy(i,j,k) = Hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
             enddo
        enddo
    enddo
!xh-PML4 loop(+)
    do k = 1,nz-1
        do j = 1,ny-1
            do i = nx-ncpml,nx-2
                psi_Hzx1(i,j,k) = bh_x(i) * psi_Hzx1(i,j,k)&
                                + ch_x(i) * ( c1*Ey(i+1,j,k)-c1*Ey(i,j,k) + c2*Ey(i+2,j,k)-c2*Ey(i-1,j,k) ) / dx
                psi_Hyx1(i,j,k) = bh_x(i) * psi_Hyx1(i,j,k)&
                                + ch_x(i) * ( c1*Ez(i+1,j,k)-c2*Ez(i,j,k) + c2*Ez(i+2,j,k)-c2*Ez(i-1,j,k) ) / dx
                Hz(i,j,k) = Hz(i,j,k) - db_z(i,j,k)*psi_Hzx1(i,j,k)
                Hy(i,j,k) = Hy(i,j,k) + db_y(i,j,k)*psi_Hyx1(i,j,k)
           enddo
        enddo
    enddo

!yh-PML4 loop(-)
    do k = 1,nz-1
        do j = 2,ncpml
            do i = 1,nx-1
                psi_Hxy1(i,j,k) = bh_y(j) * psi_Hxy1(i,j,k)&
                                + ch_y(j) * ( c1*Ez(i,j+1,k)-c1*Ez(i,j,k) + c2*Ez(i,j+2,k)-c2*Ez(i,j-1,k) ) / dy
                psi_Hzy1(i,j,k) = bh_y(j) * psi_Hzy1(i,j,k)&
                                + ch_y(j) * ( c1*Ex(i,j+1,k)-c1*Ex(i,j,k) + c2*Ex(i,j+2,k)-c2*Ex(i,j-1,k) ) / dy
                Hx(i,j,k) = Hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
                Hz(i,j,k) = Hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
            enddo
         enddo
    enddo
!!yh-PML4 loop(+)
    do k = 1,nz-1
        do j = ny-ncpml,ny-2
            do i = 1,nx-1
                psi_Hxy1(i,j,k) = bh_y(j) * psi_Hxy1(i,j,k)&
                                + ch_y(j) * ( c1*Ez(i,j+1,k)-c1*Ez(i,j,k) + c2*Ez(i,j+2,k)-c2*Ez(i,j-1,k) ) / dy
                psi_Hzy1(i,j,k) = bh_y(j) * psi_Hzy1(i,j,k)&
                                + ch_y(j) * ( c1*Ex(i,j+1,k)-c1*Ex(i,j,k) + c2*Ex(i,j+2,k)-c2*Ex(i,j-1,k) ) / dy
                Hx(i,j,k) = Hx(i,j,k) - db_x(i,j,k)*psi_Hxy1(i,j,k)
                Hz(i,j,k) = Hz(i,j,k) + db_z(i,j,k)*psi_Hzy1(i,j,k)
            enddo
         enddo
   enddo

!zh-PML4 loop(-)
    do k = 2,ncpml
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Hyz1(i,j,k) = bh_z(k) * psi_Hyz1(i,j,k)&
                                + ch_z(k) * ( c1*Ex(i,j,k+1)-c1*Ex(i,j,k) + c2*Ex(i,j,k+2)-c2*Ex(i,j,k-1) ) / dz
                psi_Hxz1(i,j,k) = bh_y(k) * psi_Hxz1(i,j,k)&
                                + ch_z(k) * ( c1*Ey(i,j,k+1)-c1*Ey(i,j,k) + c2*Ey(i,j,k+2)-c2*Ey(i,j,k-1) ) / dz
                Hy(i,j,k) = Hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
                Hx(i,j,k) = Hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
             enddo
        enddo
    enddo
!zh-PML4 loop(+)
    do k = nz-ncpml,nz-2
        do j = 1,ny-1
            do i = 1,nx-1
                psi_Hyz1(i,j,k) = bh_z(k) * psi_Hyz1(i,j,k)&
                                + ch_z(k) * ( c1*Ex(i,j,k+1)-c1*Ex(i,j,k) + c2*Ex(i,j,k+2)-c2*Ex(i,j,k-1) ) / dz
                psi_Hxz1(i,j,k) = bh_y(k) * psi_Hxz1(i,j,k)&
                                + ch_z(k) * ( c1*Ey(i,j,k+1)-c1*Ey(i,j,k) + c2*Ey(i,j,k+2)-c2*Ey(i,j,k-1) ) / dz
                Hy(i,j,k) = Hy(i,j,k) - db_y(i,j,k)*psi_Hyz1(i,j,k)
                Hx(i,j,k) = Hx(i,j,k) + db_x(i,j,k)*psi_Hxz1(i,j,k)
             enddo
        enddo
    enddo

end subroutine CPML_H4



