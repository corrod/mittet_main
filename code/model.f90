!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!モデル設定
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine model(sig,myu)
    use const_para
    implicit none

    real(8),intent(out) :: sig(nx,ny,nz)
    real(8),intent(out) :: myu(nx,ny,nz)

    !海水一様モデル
    sig(1:nx,1:ny,1:nz) = sigwa
    myu(1:nx,1:ny,1:nz) = myuwa

    !海水
    !sig()=sigwa
    !myu()=myuwa

!     !鉄板
!     do k= z0+10,z0+16
!         do j= y0-3,y0+3
!             do i=x0-3,x0+3
!         sig(i,j,k) = sigfe
!         myu(i,j,k) = myufe
!             enddo
!         enddo
!     enddo

            end subroutine model
