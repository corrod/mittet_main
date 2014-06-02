!!!モデル設定********************************************************************
subroutine model(sigma,myu)
    use const_para
    implicit none

    real(8),intent(out) :: sigma(nx,ny,nz)
    real(8),intent(out) :: myu(nx,ny,nz)

    !海水一様モデル
    sigma(1:nx,1:ny,1:nz) = sigmawa
    myu(1:nx,1:ny,1:nz) = myuwa

    !海水
    !sigma()=sigmawa
    !myu()=myuwa

    !鉄板
    !sigma()=sigmafe
    !myu()=myufe
            endsubroutine model

