!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!モデル設定
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine model_nocrack

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0

    !欠陥1つモデル
    xstart = xcenter - 2
    xend   = xcenter + 2
    ystart = ycenter - 2
    yend   = ycenter + 2
    zstart = plate - 2
    zend   = plate

    write(*,*) '************** model_nocrack ****************'
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

!#############################################################
    ! model 立方体形状欠陥モデル（空気層無し）_______________
!#############################################################

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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_nocrack




subroutine model_simple

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0

    !欠陥1つモデル
    xstart = xcenter - 2
    xend   = xcenter + 2
    ystart = ycenter - 2
    yend   = ycenter + 2
    zstart = plate - 2
    zend   = plate

    write(*,*) '************** model_simple ****************'
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

!#############################################################
    ! model 立方体形状欠陥モデル（空気層無し）_______________
!#############################################################

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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_simple



















! H23モデル　モデル1

subroutine model_stair

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0


    !実験モデル1　階段
    xblock1 = ncpml
    zblock1 = plate
    xblock2 = xblock1 + 20
    zblock2 = plate - 1
    xblock3 = xblock2 + 20
    zblock3 = plate - 2
    xblock4 = xblock3 + 20
    zblock4 = plate - 1
    xblock5 = xblock4 + 20
    zblock5 = plate

    write(*,*) '************** model_stair ****************'

    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset

    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3



!#############################################################
    ! 実験モデル1
!     5ブロック
!     幅   200mm ずつ
!     厚み 15mm 10mm 5mm 7.5mm 15mm
!#############################################################

    !海水
    do k=1,nz
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    !鉄板 ブロックを5つに分けて
    do k=1,zblock1
        do j=1,ny
            do i=1,xblock2
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock2
        do j=1,ny
            do i=xblock2,xblock3
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock3
        do j=1,ny
            do i=xblock3,xblock4
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock4
        do j=1,ny
            do i=xblock4,xblock5
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock5
        do j=1,ny
            do i=xblock5,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo



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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_stair




! H23モデル　モデル1 large

subroutine model_stair_large

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0


    !実験モデル1　階段
    xblock1 = ncpml
    zblock1 = plate
    xblock2 = xblock1 + 20
    zblock2 = plate - 1
    xblock3 = xblock2 + 20
    zblock3 = plate - 2
    xblock4 = xblock3 + 20
    zblock4 = plate - 1
    xblock5 = xblock4 + 20
    zblock5 = plate

    write(*,*) '************** model_stair ****************'

    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset

    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3



!#############################################################
    ! 実験モデル1
!     5ブロック
!     幅   200mm ずつ
!     厚み 15mm 10mm 5mm 7.5mm 15mm
!#############################################################

    !海水
    do k=1,nz
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    !鉄板 ブロックを5つに分けて
    do k=1,zblock1
        do j=1,ny
            do i=1,(nx+1)/2-30
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock2
        do j=1,ny
            do i=(nx+1)/2-30+1,(nx+1)/2-10
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock3
        do j=1,ny
            do i=(nx+1)/2-10+1,(nx+1)/2+10
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock4
        do j=1,ny
            do i=(nx+1)/2+10+1,(nx+1)/2+30
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo
    do k=1,zblock5
        do j=1,ny
            do i=(nx+1)/2+30,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo



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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_stair_large






! モデル2　ピンホール富山　一部
!#############################################################
    ! 実験モデル2 穴モデル 深さ8mm
    !     5箇所傷
    !     厚み 15mm
    !     欠陥深さ 8mm (15mm,4mm)
!#############################################################

subroutine model_pinhole8

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0


    !実験モデル2　穴
    zs  = plate - 2 !8mm
    ze  = plate
    !円
    xs1 = 10 - 5 !100mm
    xe1 = 10 + 5
    ys1 = ycenter - 5
    ye1 = ycenter + 5

    xs2 = 25 - 1!10mm
    xe2 = 25 + 1
    ys2 = ycenter - 1
    ye2 = ycenter + 1

    xs3 = 40 - 2!20mm
    xe3 = 40 + 2
    ys3 = ycenter - 2
    ye3 = ycenter + 2

    xs4 = 55 - 3!50mm
    xe4 = 55 + 3
    ys4 = ycenter - 3
    ye4 = ycenter + 3

    !直方形状
    xs5 = 85 - 7!150mm
    xe5 = 85 + 7
    ys5 = ycenter - 1!10mm??
    ye5 = ys5 + 1

    write(*,*) '************** model_pinhole8 ****************'

    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset

    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3




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
    do k=1,plate
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo

    !欠陥
    do k=zs,ze
        do j=ys1,ye1
            do i=xs1,xe1
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys2,ye2
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys3,ye3
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys4,ye4
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys5,ye5
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo


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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_pinhole8




subroutine model_pinhole8_large

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0


    !実験モデル2　穴
    zs  = plate - 2 !8mm
    ze  = plate
    !円
    xs1 = 10 - 5  +15!100mm
    xe1 = 10 + 5  +15 
    ys1 = ycenter - 5
    ye1 = ycenter + 5

    xs2 = 25 - 1  +15!10mm
    xe2 = 25 + 1  +15
    ys2 = ycenter - 1
    ye2 = ycenter + 1

    xs3 = 40 - 2  +15!20mm
    xe3 = 40 + 2  +15
    ys3 = ycenter - 2
    ye3 = ycenter + 2

    xs4 = 55 - 3  +15!50mm
    xe4 = 55 + 3  +15
    ys4 = ycenter - 3
    ye4 = ycenter + 3

    !直方形状
    xs5 = 85 - 7  +15!150mm
    xe5 = 85 + 7  +15
    ys5 = ycenter - 1!10mm??
    ye5 = ys5 + 1

    write(*,*) '************** model_pinhole8 ****************'

    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset

    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3




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
    do k=1,plate
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo

    !欠陥
    do k=zs,ze
        do j=ys1,ye1
            do i=xs1,xe1
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys2,ye2
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys3,ye3
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys4,ye4
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys5,ye5
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo


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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_pinhole8_large





subroutine model_pinhole3

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0

    !欠陥1つモデル
    xstart = xcenter - 2
    xend   = xcenter + 2
    ystart = ycenter - 2
    yend   = ycenter + 2
    zstart = plate - 2
    zend   = plate

    !実験モデル2　穴
    zs  = plate - 2 !8mm
    ze  = plate
    !円
    xs1 = 10 - 5  +15!100mm
    xe1 = 10 + 5  +15 
    ys1 = ycenter - 5
    ye1 = ycenter + 5

    xs2 = 25 - 1  +15!10mm
    xe2 = 25 + 1  +15
    ys2 = ycenter - 1
    ye2 = ycenter + 1

    xs3 = 40 - 2  +15!20mm
    xe3 = 40 + 2  +15
    ys3 = ycenter - 2
    ye3 = ycenter + 2

    xs4 = 55 - 3  +15!50mm
    xe4 = 55 + 3  +15
    ys4 = ycenter - 3
    ye4 = ycenter + 3

    !直方形状
    xs5 = 85 - 7  +15!150mm
    xe5 = 85 + 7  +15
    ys5 = ycenter - 1!10mm??
    ye5 = ys5 + 1

    write(*,*) '************** model_pinhole8 ****************'

    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset

    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3




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
    do k=1,plate
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo


    !欠陥
    do k=zstart,zend
        do j=ystart,yend
            do i=xstart,xend
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    do k=zstart,zend
        do j=ystart,yend
            do i=xstart-25,xend-25
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    do k=zstart,zend
        do j=ystart,yend
            do i=xstart+25,xend+25
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
!     !欠陥
! !     do k=zs,ze
! !         do j=ys1,ye1
! !             do i=xs1,xe1
! !                 sig(i,j,k) = sigwa
! !                 myu(i,j,k) = myuwa
! !             enddo
! !         enddo
! !     enddo
!     do k=zs,ze
!         do j=ys2,ye2
!             do i=xs2,xe2
!                 sig(i,j,k) = sigwa
!                 myu(i,j,k) = myuwa
!             enddo
!         enddo
!     enddo
!     do k=zs,ze
!         do j=ys3,ye3
!             do i=xs3,xe3
!                 sig(i,j,k) = sigwa
!                 myu(i,j,k) = myuwa
!             enddo
!         enddo
!     enddo
!     do k=zs,ze
!         do j=ys4,ye4
!             do i=xs4,xe4
!                 sig(i,j,k) = sigwa
!                 myu(i,j,k) = myuwa
!             enddo
!         enddo
!     enddo
! !     do k=zs,ze
! !         do j=ys5,ye5
! !             do i=xs5,xe5
! !                 sig(i,j,k) = sigwa
! !                 myu(i,j,k) = myuwa
! !             enddo
! !         enddo
! !     enddo


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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_pinhole3



! モデル2　ピンホール富山　全部
!#############################################################
    ! 実験モデル2 穴モデル 深さ4, 8, 15mmmm
    !     5箇所傷
    !     厚み 15mm
!#############################################################


subroutine model_pinholeall

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0

    !実験モデル2　穴 深さ8mmver
    zs  = plate - 2 !8mm
    ze  = plate
    !円
    xs1 = 10 - 5 !100mm
    xe1 = 10 + 5
    ys1 = ycenter - 5
    ye1 = ycenter + 5

    xs2 = 25  !10mm
    xe2 = 25 + 1
    ys2 = ycenter - 1
    ye2 = ys2 + 1

    xs3 = 40 - 1!20mm
    xe3 = 40 + 1
    ys3 = ycenter - 1
    ye3 = ycenter + 1

    xs4 = 55 - 2!50mm
    xe4 = 55 + 3
    ys4 = ycenter - 2
    ye4 = ycenter + 3

    !直方形状
    xs5 = 85 - 7!150mm
    xe5 = 85 + 7
    ys5 = ycenter - 1!10mm??
    ye5 = ys5 + 1


    write(*,*) '************** model_pinholeall ****************'

    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset

    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3


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
    do k=1,plate
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo

    !欠陥8mm 中
!     do k=zs,ze
!         do j=ys1,ye1
!             do i=xs1,xe1
!                 sig(i,j,k) = sigwa
!                 myu(i,j,k) = myuwa
!             enddo
!         enddo
!     enddo
    do k=zs,ze
        do j=ys2,ye2
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys3,ye3
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys4,ye4
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys5,ye5
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    !欠陥15mm貫通　上
    do k=1,plate
        do j=ys1+15+4,ye1+15-5
            do i=xs1+4,xe1-5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys2+15,ye2+15
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys3+15,ye3+15
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys4+15,ye4+15
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys5+15,ye5+15
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    !欠陥4mm貫通　下
    do k=plate-1,plate
        do j=ys1-15,ye1-15
            do i=xs1,xe1
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys2-15,ye2-15
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys3-15,ye3-15
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys4-15,ye4-15
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys5-15,ye5-15
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_pinholeall




subroutine model_pinholeall_large

    use const_para
    implicit none

    integer :: airlayer
    integer :: xcenter, xstart, xend
    integer :: ycenter, ystart, yend
    integer :: zcenter, zstart, zend

    integer :: xblock1,xblock2,xblock3,xblock4,xblock5
    integer :: zblock1,zblock2,zblock3,zblock4,zblock5

    integer :: xs1,xs2,xs3,xs4,xs5
    integer :: xe1,xe2,xe3,xe4,xe5
    integer :: ys1,ys2,ys3,ys4,ys5
    integer :: ye1,ye2,ye3,ye4,ye5
    integer :: zs,ze

    airlayer = nz - 20

    xcenter = x0
    ycenter = y0
    zcenter = z0

    !実験モデル2　穴 深さ8mmver
    zs  = plate - 2 !8mm
    ze  = plate
    !円
    xs1 = 10 - 5   +15 !100mm
    xe1 = 10 + 5   +15
    ys1 = ycenter - 5
    ye1 = ycenter + 5

    xs2 = 25   +15  !10mm
    xe2 = 25 + 1   +15
    ys2 = ycenter - 1
    ye2 = ys2 + 1

    xs3 = 40 - 1   +15!20mm
    xe3 = 40 + 1   +15
    ys3 = ycenter - 1
    ye3 = ycenter + 1

    xs4 = 55 - 2   +15!50mm
    xe4 = 55 + 3   +15
    ys4 = ycenter - 2
    ye4 = ycenter + 3

    !直方形状
    xs5 = 85 - 7   +15!150mm
    xe5 = 85 + 7   +15
    ys5 = ycenter - 1!10mm??
    ye5 = ys5 + 1


    write(*,*) '************** model_pinholeall ****************'

    write(*,*) 'plate_z', plate
    write(*,*) 'offset', offset

    write(*,*) 'x,y,z source',x_source,y_source,z_source
    write(*,*) 'x,y,z source2',x_source2,y_source2,z_source2
    write(*,*) 'xx1,yy1,zz1',xx1,yy1,zz1
    write(*,*) 'xx2,yy2,zz2',xx2,yy2,zz2
    write(*,*) 'xx3,yy3,zz3',xx3,yy3,zz3


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
    do k=1,plate
        do j=1,ny
            do i=1,nx
                sig(i,j,k) = sigfe
                myu(i,j,k) = myufe
            enddo
        enddo
    enddo

    !欠陥8mm 中
!     do k=zs,ze
!         do j=ys1,ye1
!             do i=xs1,xe1
!                 sig(i,j,k) = sigwa
!                 myu(i,j,k) = myuwa
!             enddo
!         enddo
!     enddo
    do k=zs,ze
        do j=ys2,ye2
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys3,ye3
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys4,ye4
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=zs,ze
        do j=ys5,ye5
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    !欠陥15mm貫通　上
    do k=1,plate
        do j=ys1+15+4,ye1+15-5
            do i=xs1+4,xe1-5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys2+15,ye2+15
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys3+15,ye3+15
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys4+15,ye4+15
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=1,plate
        do j=ys5+15,ye5+15
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo

    !欠陥4mm貫通　下
    do k=plate-1,plate
        do j=ys1-15,ye1-15
            do i=xs1,xe1
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys2-15,ye2-15
            do i=xs2,xe2
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys3-15,ye3-15
            do i=xs3,xe3
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys4-15,ye4-15
            do i=xs4,xe4
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
    do k=plate-1,plate
        do j=ys5-15,ye5-15
            do i=xs5,xe5
                sig(i,j,k) = sigwa
                myu(i,j,k) = myuwa
            enddo
        enddo
    enddo
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
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(102,*) i,j,k,sig(i,j,k)
                enddo
            enddo
        close(102)

        open(103,file='model_myu_xy.dat')
            k = plate
            do j = 1, ny
                do i = 1, nx
                    write(103,*) i,j,k,myu(i,j,k)
                enddo
            enddo
        close(103)

        end subroutine model_pinholeall_large










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


