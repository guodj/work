pro solar_flux, Temp, DayOfYear=DayOfYear, Tilt=Tilt, Time=Time, n=n, plt=plt

	common map_c,map

	if 0 then begin
	  DayOfYear 	= 90
	  Tilt		= 23.8
	  Time		= 24.0
	endif

	if n_elements(DayOfYear) ne 1 then DayOfYear = 90.
	if n_elements(Tilt) ne 1 then Tilt = 23.8
	if n_elements(Time) ne 1 then Time = 24.0
	if n_elements(n) ne 1 then n = 200
	if n_elements(plt) ne 1 then plt = 1

	m		= n/2
	Latitude	= replicate(1.,n)#(findgen(m)/(m-1)*180.-90.)
	Longitude	= (findgen(n)/(n-1)*360.)#replicate(1.,m)


	;*********************************
	;*** Convert to Usable Varibles ***
	;*********************************
	
	DOY	= float(DayOfYear-89+float(Time)/24.)/356. * 2. * !pi
	TLT	= Tilt / !radeg
	CLT	= (90. - Latitude) / !radeg
	LNG	= Longitude / !radeg
	TIM	= float(Time)/24. * 2. * !pi

	SDOY	= sin(DOY)
	CDOY	= cos(DOY)
	STLT	= sin(TLT)
	CTLT	= cos(TLT)
	SCLT	= sin(CLT)
	CCLT	= cos(CLT)
	SLNG	= sin(LNG-DOY+TIM+!pi/2)
	CLNG	= cos(LNG-DOY+TIM+!pi/2)



	CX 	= CLNG * SCLT * CDOY * CTLT $
	  + SLNG * SCLT * (-SDOY) $
	  + CCLT * CDOY * STLT

	CY 	= CLNG * SCLT * SDOY * CTLT $
	  + SLNG * SCLT * CDOY $
	  + CCLT * SDOY * STLT

	CZ	= CLNG * SCLT * (-STLT) + CCLT * CTLT

	HERE	= where(CY < 0.)
	if HERE(0) ne -1 then CY(HERE) = 0.

	FLUX	= CY * SCLT

	en	= 400/n
;	CXB	= rebin(CX,en*n,en*m)
	CYB	= rebin(CY,en*n,en*m)
;	CZB	= rebin(CZ,en*n,en*m)
	FLX	= rebin(FLUX,en*n,en*m)

	dd=DayOfYear
	caldat,dd,mo,d,y & xyouts,5,5,scr(mo)+'/'+scr(d),/dev,chars=.5

	tcyb 	= fltarr(en*m) 		;total
	tflx 	= fltarr(en*m)		;integrated
	thdl 	= fltarr(en*m)		;total number of daylight
	for i=0,en*m-1 do begin
	  tcyb(en*m-1-i) = total(cyb(*,i))/(en*n)
	  here = where(cyb(*,i) ne 0)
	  if here(0) eq -1 then hdl = 0 else $
	    hdl = 24.*n_elements(here)
	  thdl(en*m-1-i) = hdl/(en*n)
	  tflx(en*m-1-i) = total(flx(*,i))/(en*n)
	endfor


	Temp = ((1.+2.)*1.4e16 * (1. - 0.33) * tcyb / (2. * 5.67e5) )^(0.25)

	if plt then begin
	  if !d.x_size ne en*n+200 or !d.y_size ne 2*en*m then $
	    window,xsi=en*n+200,ysi=2*en*m else erase

	  if n_elements(map) lt 2 then begin
	    loadct,3
	    map_set,/cont,pos=[0,0,1,1]*([0,0,1.*en*n/!d.x_size,0] $
	      + [0,0,0,1.*en*m/!d.y_size])
	    map=tvrd(0,0,en*n,en*m)
	    map=where(map ne 0)
	  endif

;	  *** Top is (surface normal . sun earth line direction) ***

	  im = float(bytscl(CYB))/255*(!d.table_size-2)
	  im(map)=255
	  tv,im,0,en*m

	  plot_sum_or_daylight = 0
	  if plot_sum_or_daylight then begin
	    var = tcyb
	    xyout = 'SUM(norm.sunline)'
	    xra=[min(tcyb),max(tcyb)]
	  endif else begin
	    var = thdl
	    xyout = 'Hours of Daylight'
	    xra = [0,24]
	  endelse

	  plot,var,rebin(reform(clt(0,*),m),en*m),/noe, $
	    pos=[en*n,en*m,en*n+100,2*en*m], $
	    /dev,xst=5,yst=5, xra=xra, yra=mm(clt(0,*))

	  xyouts,en*n+5,en*m+5,/dev,xyout,chars=.5

	  plot,temp,rebin(reform(clt(0,*),m),en*m),/noe, $
	    pos=[en*n+100,en*m,en*n+200,2*en*m], $
	    /dev,xst=5,yst=5, yra=mm(clt(0,*)), xra=[0,350]

	  xyouts,en*n+105,en*m+5,/dev,'B.B. Temp [0,350]',chars=.5

;	  *** Bottom is Unatenuated Flux ***

	  im = float(bytscl(FLX))/255*(!d.table_size-2)
	  im(map)=255
	  tv,im

	  plot,tflx,rebin(reform(clt(0,*),m),en*m),/noe, $
	    pos=[en*n,0,en*n+200,en*m], $
	    /dev,xst=5,yst=5, yra=mm(clt(0,*))

	  xyouts,en*n+5,5,/dev,'SUM(FLUX)',chars=.5

	endif
end

	common map_c,map

	map = -1

	n = 200
	m = n/2
	en = 400/n
	nf = 15

	plt = 1

	if plt then begin
;	  xanimate, set=[n*en+200,2*m*en,nf]
	  for i=0.,nf-1 do begin
	    day=365./nf*i
	    solar_flux,day=day,time=16
;	    xanimate,frame=i,window=!d.window
	  endfor
	  color_flip
;	  xanimate
	endif else begin
	  if n_elements(temp) eq 0 then begin
	    nf = n
	    for i=0.,nf-1 do begin
	      day=365./nf*i
	      print,day
	      solar_flux,day=day,time=16,TTEMP,plt=plt
	      if i eq 0 then Temp = TTEMP else Temp=[[Temp],[TTEMP]]
	    endfor
	    temp=rotate(temp,3)
	    temp=9./5.*(temp-273.15)+32 ;F
	  endif

	  tv_image,temp,[0,365],[-90,90],Title='BB Temp'

	  f = 5
	  i=indgen(n/f)*f#replicate(1,n/f)
	  j=replicate(1,n/f)#indgen(n/f)*f
	  msurface,rotate(temp(i,j),3)

	  lat=0.
	  read,'Choose a latitude: ',lat
	  if lat gt 90 or lat lt -90 then print,'Dumb Ass!'
	  if lat gt 90 then lat = 90.
	  if lat lt -90 then lat = -90.

	  clat = ( 90 - lat ) 
	  i = (180.-clat)/180.*(n-1)
	  plot,temp(*,i),yra=mm(temp),xst=1,xticks=12,xtickn= $
	    ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug', $
	    'Sep','Oct','Nov','Dev',' '],xtitle='Month',Ytitle='BB Temp (F)'
	  mx = max(temp(*,i))
	  plots,[0,n-1],mx*[1,1],lin=1
	  xyouts,-n/40,mx,ali=1,scr(mx),chars=.5
	  mn = min(temp(*,i))
	  plots,[0,n-1],mn*[1,1],lin=1
	  xyouts,-n/40,mn,ali=1,scr(mn),chars=.5
	endelse	

end
