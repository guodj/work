pro color
	If n_elements(ct) ne 1 then ct=15
	num=!d.table_size

;	***************** color 0 *********************
	red 	= 0
	green 	= 0
	blue 	= 0

;	***************** needed colors *********************

	red	= [red,	    0,  0,150, 30,  0,200,100]
	green	= [green,   0,  0,150, 30,150,  0,  0]
	blue	= [blue,  250,190,150, 30,  0,  0,  0]

	red	= [red,	    0,0]
	green	= [green,   0,0]
	blue	= [blue,    0,0]

	red	= [red,	    0,  0,  0,  0,200,150,100]
	green	= [green, 200,150,100, 50,  0,  0,  0]
	blue	= [blue,    0,  0,  0,  0,200,150,100]

	red	= [red,	  250,250,190,130]
	green	= [green, 250,  0,  0,  0]
	blue	= [blue,  250,  0,  0,  0]

	red	= [red,	  250,150,120, 180,150,230,200,170,150,110]
	green	= [green, 250,150,120,   0,  0,230,  0,170,  0,110]
	blue	= [blue,  250,255,120, 180,150,  0,  0,  0,  0,  0]

	red	= [red,	    0,  0,  0,  0,  0,  0,  0,  0]
	green	= [green,   0,  0,  0,  0,  0,  0,  0,  0]
	blue	= [blue,  255,252,244,228,200,165, 125, 70]


	red	= [red,	    0,  0,  0,255]
	green	= [green,  70,130,255,255]
	blue	= [blue,    0,  0,  0,255]

;  spider colors (43 - 45) :

	red	= [red,	    0,  0,255]
	green	= [green, 100,255,  0]
	blue	= [blue,    0,100,  0]


;	***************** remaining colors *********************

	n = num - n_elements(red) - 3
	i = findgen(n)/(n-1)

	loadct,0
	tvlct,r,g,b,/get

	r = r(i*!d.table_size)
	b = b(i*!d.table_size)
	g = g(i*!d.table_size)

	red	= [red,0,r]
	green	= [green,0,g]
	blue	= [blue,0,b]

;	***************** remaining colors *********************

	rem = replicate(255,256-n_elements(red))
	red	= [red,rem]
	green	= [green,rem]
	blue	= [blue,rem]

	tvlct,red,green,blue

        return
end
