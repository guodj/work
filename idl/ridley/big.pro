function big,im,x

	si	= size(im)

	if si(0) ne 2 then return,im 
	if n_e(x) ne 1 then return,im 

	n	= si(1)
	m	= si(2)

        i       = long(findgen(n*x)/(n*x-1)*n)
        j       = long(findgen(m*x)/(m*x-1)*m) * n
        i       = i # replicate(1,m*x)
        j       = replicate(1,n*x) # j

        return, im(i+j)

end
