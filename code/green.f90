!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!グリーン関数(周波数領域)を求める
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program green
	use const_para
		implicit none

		integer :: nd,ios
		real(8), allocatable :: ii1(:), ii2(:)
		real(8), allocatable :: inp1_r(:), inp1_i(:),inp2_r(:), inp2_i(:)
		complex(kind(0d0)), allocatable :: inp1(:), inp2(:)
		complex(kind(0d0)), allocatable :: out1(:)

		
		!E inp1=E(omega)の読み込み!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		open(1,file='inp1.dat',action='read')
	    nd=0
	    do
	        read(1,'(f12.0)',iostat=ios)
	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
	        nd=nd+1
	    enddo
	    close(1)

	    allocate(ii1(1:nd),inp1_r(1:nd),inp1_i(1:nd),inp1(1:nd))

		open(1,file='inp1.dat')
			read(1,*) ii1(1:nd), inp1_r(1:nd), inp1_i(1:nd)
			inp1(1:nd) = inp1_r(1:nd) + (0.0d0,1.0d0)*inp1_i(1:nd)
		close(1)
		


		!Jn inp2=J(omega)の読み込み!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		open(2,file='inp2.dat',action='read')
	    nd=0
	    do
	        read(2,'(f12.0)',iostat=ios)
	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
	        nd=nd+1
	    enddo
	    close(2)

		allocate(ii2(1:nd), inp2_r(1:nd),inp2_i(1:nd),inp2(1:nd), out1(1:nd))

		open(2,file='inp2.dat')
			read(2,*) ii2(1:nd), inp2_r(1:nd), inp2_i(1:nd)
			inp2(1:nd) = inp2_r(1:nd) + (0.0d0,1.0d0)*inp2_i(1:nd)
		close(2)



		!グリーン関数導出 out1=Gej_in(x,omega|xs)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		open(3,file='out1.dat')
		out1(1:nd) = inp1(1:nd) / inp2(1:nd)  !!G=E/J
		
		do i=1,nd
			write(3,*) i,real(out1(i)), aimag(out1(i))
		enddo
		close(3)

		deallocate(ii1,inp1,ii2,inp2) !!!取り扱い注意

			end program green