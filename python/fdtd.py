for k in xrange(nz)
for j in xrange(ny)
for i in xrange(nx)
Ex[i, j, k] = ca_x[i, j, k] + Ex[i, j, k]\
+ cb_x[i, j, k] * ((c1*Hz(i, j, k) - c1*Hz(i, j-1, k) + c2*Hz(i,j+1,k) - c2*hz(i, j-2, k)))
