program writefile
	implicit none
	open(1,file='1.d')
	write(1,*) "1"
	close(1)

	open(1,file="2.d")
	write(1,*) '2'
	close(1)
end program writefile