!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!モデル設定
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine model
    use const_para
    implicit none

        !     real(8),intent(out) :: sig(nx,ny,nz)
        !     real(8),intent(out) :: myu(nx,ny,nz)

            !海水一様モデル
        !     sig(1:nx,1:ny,1:nz) = sigwa
        !     myu(1:nx,1:ny,1:nz) = myuwa

!海水
sig(1:nx,1:ny,1:nz) = sigwa
myu(1:nx,1:ny,1:nz) = myuwa

!鉄板
sig(1:nx,1:ny,nz-10:nz) = sigfe
myu(1:nx,1:ny,nz-10:nz) = myufe

!欠陥
sig(x0-5:x0+5,y0-5:y0+5,nz-10:nz-5) = sigwa
myu(x0-5:x0+5,y0-5:y0+5,nz-10:nz-5) = myuwa

end subroutine model
