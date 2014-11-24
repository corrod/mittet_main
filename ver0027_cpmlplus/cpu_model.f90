!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!モデル設定
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine model

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0


    xstart = xcenter - 2
    xend   = xcenter + 2
    ystart = ycenter - 5
    yend   = ycenter + 5
    zstart = plate - 2
    zend   = plate


    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset
    write(*,*) 'xstart,xend', xstart,xend
    write(*,*) 'ystart,yend', ystart,yend
    write(*,*) 'zstart,zend', zstart,zend
    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3


    ! model 立方体形状欠陥モデル（空気層無し）_______________

    !海水
    do k=1,nz
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    !鉄板
    do k=1,zend
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
!             !海水ー鉄板境界
!             k = zend +1
!             do j = 1,ny
!                 do i = 1,nx
!                     sig(i,j,k) = (sigwa + sigfe) * 0.5d0
!                     myu(i,j,k) = (myuwa + myufe) * 0.5d0
!                 enddo
!             enddo

    !欠陥
    do k=zstart,zend
        do j=ystart,yend
            do i=xstart,xend
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

!     do k=1,nz
!         do j=1,ny
!             do i=1,nx
!                 sig(i,j,k) = sigfe
!                 myu(i,j,k) = myufe
!                 if(k>plate) then
!                     sig(i,j,k) = sigwa
!                     myu(i,j,k) = myuwa
!                     elseif(i>=xstart .and. i<=xend &
!                      .and. j>=ystart .and. j<=yend &
!                      .and. k>=zstart .and. k<= zend) then
!                     sig(i,j,k) = sigwa
!                     myu(i,j,k) = myuwa
!                 endif
!             enddo
!         enddo
!     enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!モデル出力
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !x-z 鉛直断面
        open(100,file='model_sig_xz.dat')
            do k = 1, nz
                 j = y0
                    do i = 1, nx
                        write(100,*) i,j,k,sig(i,j,k)
                    enddo
            enddo
        close(100)

        open(101,file='model_myu_xz.dat')
            do k = 1, nz
                 j = y0
                    do i = 1, nx
                        write(101,*) i,j,k,myu(i,j,k)
                    enddo
            enddo
        close(101)

        !x-y 水平断面
        open(102,file='model_sig_xy.dat')
            k = plate-1
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate-1
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model


    ! ! model1 海水一様モデル______________________
    ! ! !海水
    ! sig(1:nx, 1:ny, 1:nz) = sigwa
    ! myu(1:nx, 1:ny, 1:nz) = myuwa


    ! model 海水-空気______________________________
    ! 海水
!     do k=1,airlayer-1
!         do j=1,ny
!             do i=1,nx
!                 sig(i,j,k) = sigwa
!                 myu(i,j,k) = myuwa
!             enddo
!         enddo
!     enddo
!     ! 空気
!     do k=airlayer,nz
!         do j=1,ny
!             do i=1,nx
!                 sig(i,j,k) = sigair
!                 myu(i,j,k) = myuair
!             enddo
!         enddo
!     enddo

    ! model2 海水-鉄板 欠陥なしモデル_____________________
    ! !海水
    ! sig(1:nx, 1:ny, 1:nz) = sigwa
    ! myu(1:nx, 1:ny, 1:nz) = myuwa

    ! !鉄板
    ! sig(1:nx, 1:ny, z0+10:nz) = sigfe
    ! myu(1:nx, 1:ny, z0+10:nz) = myufe


! ! model2-2 欠陥なしモデル（空気層あり）_____________
! ! 海水
!     do k=1,plate
!         do j=1,ny
!             do i=1,nx
!                 sig(i,j,k) = sigwa
!                 myu(i,j,k) = myuwa
!             enddo
!         enddo
!     enddo
! ! 鉄板
!     do k=plate+1,airlayer-1
!         do j=1,ny
!             do i=1,nx
!                 sig(i,j,k) = sigfe
!                 myu(i,j,k) = myufe
!             enddo
!         enddo
!     enddo
! ! 空気層
!     do k=airlayer,nz
!         do j=1,ny
!             do i=1,nx
!                 sig(i,j,k) = sigair
!                 myu(i,j,k) = myuair
!             enddo
!         enddo
!     enddo







! !model4 立方体形状欠陥モデル（空気層あり）________________
    ! !海水
    ! do k=1,nz
    !     do j=1,ny
    !         do i=1,nx
    !             sig(i,j,k) = sigwa
    !             myu(i,j,k) = myuwa
    !         enddo
    !     enddo
    ! enddo
    ! !鉄板
    ! do k=1,zend
    !     do j=1,ny
    !         do i=1,nx
    !             sig(i,j,k) = sigfe
    !             myu(i,j,k) = myufe
    !         enddo
    !     enddo
    ! enddo
    !             !海水ー鉄板境界
!             k = zend +1
!             do j = 1,ny
!                 do i = 1,nx
!                     sig(i,j,k) = (sigwa + sigfe) * 0.5d0
!                     myu(i,j,k) = (myuwa + myufe) * 0.5d0
!                 enddo
!             enddo
    ! !欠陥
    ! do k=zstart,zend
    !     do j=ystart,yend
    !         do i=xstart,xend
    !             sig(i,j,k) = sigwa
    !             myu(i,j,k) = myuwa
    !         enddo
    !     enddo
    ! enddo
    ! !空気層
    ! do k=airlayer,nz
    !     do j=1,ny
    !         do i=1,nx
    !             sig(i,j,k) = sigair
    !             myu(i,j,k) = myuair
    !         enddo
    !     enddo
    ! enddo

!             !空気層ー海水境界
!             k = airlayer
!             do j=1,ny
!                 do i =1,nx
!                     sig(i,j,k) = (sigwa + sigair) * 0.5d0
!                     myu(i,j,k) = (myuwa + myuair) * 0.5d0
!                 enddo
!             enddo

    ! ! model5 コーティング層挟む____________________
    ! !海水
    ! sig(1:nx, 1:ny, 1:nz) = sigwa
    ! myu(1:nx, 1:ny, 1:nz) = myuwa

    ! !コーティング
    ! sig(1:nx, 1:ny, x0+7:x0+9) = sig!　
    ! myu(1:nx, 1:ny, x0+7:x0+9) = myu!　

    ! !鉄板
    ! sig(1:nx, 1:ny, x0+10:nz) = sigfe
    ! myu(1:nx, 1:ny, x0+10:nz) = myufe

    ! !欠陥
    ! sig(x0-3:x0+3, y0-3:y0+3, x0+10:x0+15) = sigwa
    ! myu(x0-3:x0+3, y0-3:y0+3, x0+10:x0+15) = myuwa


