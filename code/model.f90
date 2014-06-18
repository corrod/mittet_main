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
    myu(1:nx,1:ny,1:nz)   = myuwa

    !海水
    !sig()=sigwa
    !myu()=myuwa

    !鉄板
    !sig()=sigfe
    !myu()=myufe
            end subroutine model
