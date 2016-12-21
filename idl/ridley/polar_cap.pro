pro fit_ellipse, x, y, as, bs, xs, ys, percent, lun=lun, time=time

  if n_elements(lun) eq 0 then lun = 0
  if n_elements(time) eq 0 then time = '             '

  print, '  Starting Ellipse Fit'

  dp = 10.0
  dv = 10.0

  xn = float(x-100)
  yn = float(y-100)

  ymi = fix(min(yn))
  yma = 0
  dy = 5.0
  d = 0

  for y = ymi, yma do begin
    loc = where(yn gt y-dy and yn lt y+dy, count)
    if count gt d then begin
      d = count
      ys = float(abs(y))
    endif
  endfor

  xmi = fix(min(xn))
  xma = 0
  dx = 5.0
  d = 0

  for x = xmi, xma do begin
    loc = where(xn gt x-dx and xn lt x+dx, count)
    if count gt d then begin
      d = count
      xs = float(abs(x))
      yo = mean(yn(loc))
    endif
  endfor

  if (xs lt ys) then begin
    xs = 1.25*xs
    loc = where(abs(xn) gt abs(xs)-dx and abs(xn) lt abs(x)+dx, count)
    if count gt 0 then yo = mean(yn(loc))
  endif

  xo = 0.0

;  print, 'Initial guess:'
;  print, '      a   : ',xs,			$
;	 '      b   : ',ys
;  print, '      xo  : ',xo,			$
;	 '      yo  : ',yo

  npts = 50

  fac = 2.0-0.00001

  n = 0
  nmax = 500
  npzs = 0
  dpar = 10.0

  a = xs
  b = ys+yo

  as = a
  bs = b
  xs = xo
  ys = yo
  abase = a
  bbase = b
  xbase = xs
  ybase = ys

  nscale = 0

  while n le nmax do begin

    xb = (findgen(npts/2+1)^2.0/float(npts/2)^2.0-1.0)*fac*a/2.0
    xt = [xb,-reverse(xb(0:npts/2-1))]+xo
    ytp = ((1.0 - ((xt-xo)/a)^2.0)^0.5)*b + yo
    ytn = -((1.0 - ((xt-xo)/a)^2.0)^0.5)*b + yo

    npa = 0
    npz = 0

    loc = where(yn gt yo,count)
    if count gt 0 then begin

      xn2 = xn(loc)
      yn2 = yn(loc)

      for i=0,npts do begin

;        diff = abs(xt(i)-xn2)
;        loc = where(diff eq min(diff))
;        d = (xt(i)-xn2(loc(0)))^2.0+(ytp(i)-yn2(loc(0)))^2.0
;        npa = npa + 1
;        if d lt dpar then npz = npz+1

        d = (xt(i)-xn2)^2.0+(ytp(i)-yn2)^2.0

        loc = where(d eq min(d))
        if d(loc(0)) lt dpar then begin
          npz = npz + 1
          xn2(loc(0)) = 500.0
          yn2(loc(0)) = 500.0
        endif

      endfor

    endif

    loc = where(yn le yo,count)
    if count gt 0 then begin

      xn2 = xn(loc)
      yn2 = yn(loc)

      for i=0,npts do begin
        d = (xt(i)-xn2)^2.0+(ytn(i)-yn2)^2.0
        loc = where(d eq min(d))
        if d(loc(0)) lt dpar then begin
          npz = npz + 1
          xn2(loc(0)) = 500.0
          yn2(loc(0)) = 500.0
        endif
      endfor

;      for i=0,npts do begin
;        diff = abs(xt(i)-xn2)
;        loc = where(diff eq min(diff))
;        d = (xt(i)-xn2(loc(0)))^2.0+(ytn(i)-yn2(loc(0)))^2.0
;        npa = npa + 1
;        if d lt dpar then npz = npz+1
;      endfor

    endif

    if npz gt npzs then begin

      xs = xo
      ys = yo
      as = a
      bs = b
      npzs = npz

;      dp = 10.0*float(2*npts-npzs)/float(2.0*npts)
;      dv = 10.0*float(2*npts-npzs)/float(2.0*npts)

;      print, n, a,b,xo,yo,npzs

      if (bs gt as) and (abs(xs) gt abs(ys)) then begin
        as = abase
        bs = bbase
	xs = xbase
	ys = ybase
        n = n - 100
        if n lt -500 then n = 501
        dv = 5.0
      endif

      if (abs(xs) gt as/3.0) then begin
        as = 3.0*abs(xs)
	xs = 0.0
        n = n - 15
        npzs = 0
      endif

      if (bs lt 0.85*abs(ys)) then begin
        bs = abs(ys)+1.0
	ys = ys/2.0
        n = n - 10
        npzs = 0
      endif

      if npzs gt npts then begin
        dp = 0.9*dp
        dv = 0.9*dv
      endif

      xo = xs
      yo = ys
      a = 0.0
      while a le 0.0 do a  = as + 0.5*dp*(randomu(s)-0.5)
      b = 0.0
      while b le 0.0 do b  = bs + 0.5*dp*(randomu(s)-0.5)

      n = n - 1

    endif else begin

      if ((n gt nmax/2) and 						$
         ((abs(xs) gt as/3.0) or ((abs(ys) gt 1.25*bs)))) then begin

        as = 1.25*abase
        bs = 1.25*bbase
        xs = xbase
        ys = ybase
	n = n - 50
	npzs = 0

	dp = 10.0
	dv = 10.0

      endif

      xo = xs + 1.0*dv*(randomu(s)-0.5)
      yo = ys + 2.0*dv*(randomu(s)-0.5)
      a = 0.0
      while a le 0.0 do a  = as + 2.0*dp*(randomu(s)-0.5)
      b = 0.0
      while b le 0.0 do b  = bs + 2.0*dp*(randomu(s)-0.5)

    endelse

    n = n + 1

  endwhile

  print, '  Ellipse Fit completed : '
  print, '      a   : ',as,			$
	 '      b   : ',bs
  print, '      xo  : ',xs,			$
	 '      yo  : ',ys
  print, '      Percent of oval : ',100.0*float(npzs)/float(2*npts)

  percent = 100.0*float(npzs)/float(2*npts)

  return

end

print, 'Enter file list which you would like plotted : '
flist = '/l2/data/cme/970110/970110_030855.bin '
read,flist

maxran=40.0

print, ''
print, 'Enter number of times to skip in between each image : '
print, '(0 for no skipping, 1 for every other image, ...)'
nskip = 0
read, nskip
nskip = nskip + 1

print, ''
print, 'Enter color table which you would like :'
ct = '/l2/software/idl/polar/geog/polar2.ct'
read,ct

print, ''
print, 'Enter maximum reyleighs for color table (200-600) :'
maxct = 150
read,maxct

print, ''
print, 'Enter minimum reyleighs for nightime auroral oval (120) :'
bound=30
read,bound


print, ''
print, '1. Polar Cap boundary only'
print, '2. Polar Cap colored only'
print, '3. Polar Cap and boundary colored'
print, 'Enter areas you would like colored : '
boundaries = 3
read, boundaries

print, ''
print, 'Enter number of plots per page : '
ppp = 4
read, ppp

print, ''
print, 'Enter ps file name (return for screen) :'
psfile = 'test.ps'
read,psfile

if strlen(psfile) gt 0 then begin
  print, 'Would you like a different ps file for each page (y/n) ?'
  que=''
  read,que
  if (strmid(que,0,1) eq 'y') then begin
    moreppp = 1
    strl = strpos(psfile,'.ps')
    if strl gt 0 then psfile = strmid(psfile,0,strl)


    outfile = psfile+'_01'+'.ps'
  endif else begin
    moreppp = 0
    outfile = psfile
  endelse
  setdevice,outfile,'land',4,0.95
endif else moreppp = 0

print, 'Enter filename for the ellipse value output (return for none) :
fout = ''
read,fout
if strlen(fout) gt 0 then begin
  lun = 12
  openw,lun,fout
endif else lun = 0

readct, ncolors, ct

read_images, flist, nskip, uttimes, inttimes, image_save, 90-maxran

nimage = n_elements(image_save(*,0,0))

space = 0.06

pos_space, ppp, space, sizes


rows=n_elements(image_save(0,0,*))
columns=n_elements(image_save(0,*,0))
cap_pixels=fltarr(columns,rows)
edge_pixels=fltarr(columns,rows)
pixels=fltarr(columns,rows)
column=fltarr(1)
row=fltarr(1)


;-----------------------------------------------------------------------------

; Produce a median filtered array of the image to remove spurious points

;-----------------------------------------------------------------------------


cimage = fltarr(columns,rows)

for i= 0, nimage-1 do begin

  cimage(*,*) = image_save(i,*,*)

  print,'Processing image  '+strmid(uttimes(i),7,4)+':'+  $
      strmid(uttimes(i),11,2)+' UT'

  cap=0
  count=0	
  cap_pixels(*,*)=0.0
  edge_pixels(*,*)=0.0
  
  cap_color=maxct 		

  for k=0,rows-1 do begin
	pixels(*,k) =float(maxran)*cimage(*,k)/float(maxran)
  endfor

  pixels=median(pixels,6)
  count=0


;-----------------------------------------------------------------------------

;  This section scans the image to locate the poleward boundary of
;  the polar cap.

;-----------------------------------------------------------------------------

  if pixels(100,100) lt bound then begin
    ibound = bound - 0.5*(bound-pixels(100,100))
    print, 'Bounding value changed from ',bound,' to ',ibound
  endif else ibound = bound

  for j=0,99 do begin
    flag=0 
    for k=0,99 do begin

	if (pixels(j+100,k+100) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(j+100,k+100) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(j+100,k+100)=cap_color
            column=[column,j+100]
	    row=[row,k+100]
	    count=count+1
  	 endif

     endfor

   endfor

    
  for j=0,99 do begin
    flag=0 
    for k=0,99 do begin

	if (pixels(100-j,k+100) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(100-j,k+100) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(100-j,k+100)=cap_color
            column=[column,100-j]
	    row=[row,k+100]
	    count=count+1
    	endif

	

    endfor
    endfor
    

  for j=0,99 do begin
    flag=0 
    for k=0,99 do begin

	if (pixels(j+100,100-k) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(j+100,100-k) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(j+100,100-k)=cap_color
            column=[column,j+100]
	    row=[row,100-k]
	    count=count+1
 	endif

    endfor
  endfor    

  for j=0,99 do begin
    flag=0 
    for k=0,99 do begin

	if (pixels(100-j,100-k) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(100-j,100-k) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(100-j,100-k)=cap_color
            column=[column,100-j]
	    row=[row,100-k]
	    count=count+1
	endif

    endfor
  endfor


;-----------------------------------------------------------------------------

  for k=0,99 do begin
    flag=0 
    for j=0,99 do begin

	if (pixels(j+100,k+100) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(j+100,k+100) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(j+100,k+100)=cap_color
	endif

    endfor
   endfor

    
  for k=0,99 do begin
    flag=0 
    for j=0,99 do begin

	if (pixels(100-j,k+100) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(100-j,k+100) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(100-j,k+100)=cap_color
	endif

	

    endfor
    endfor
    

  for k=0,99 do begin
    flag=0 
    for j=0,99 do begin

	if (pixels(j+100,100-k) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(j+100,100-k) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(j+100,100-k)=cap_color
 	endif

    endfor
  endfor    

  for k=0,99 do begin
    flag=0 
    for j=0,99 do begin

	if (pixels(100-j,100-k) lt ibound) and (flag eq 0) then begin
	    flag=1
	endif

	if (pixels(100-j,100-k) gt ibound) and (flag eq 1) then begin
            flag=2
            edge_pixels(100-j,100-k)=cap_color
	endif

    endfor
  endfor


;----------------------------------------------------------------------------

; Overlay the polar cap image or edge on top of the auroral image and count 
; pixels within the polar cap

;----------------------------------------------------------------------------


; smooth the edge_pixels to get rid of any "lone" points:

  edge_pixels=median(edge_pixels,2)
 
 for j=0,columns-1 do begin;

    for k=0,rows-1 do begin


   if (edge_pixels(j,k) eq cap_color) then begin
 

   CASE 1 OF


	(k le 100) and (j le 100): BEGIN
	   if max(pixels(0:j,k)) lt max(pixels(j:100,k)) and $
              max(pixels(j,0:k)) lt max(pixels(j,k:100)) then begin
		 edge_pixels(j,k)=0
	    endif
            END

	(k le 100) and (j gt 100): BEGIN
	   if max(pixels(j:columns-1,k)) lt max(pixels(100:j,k)) and $
	      max(pixels(j,0:k)) lt max(pixels(j,k:100)) then begin
		 edge_pixels(j,k)=0
	   endif
	   END

	(k gt 100) and (j le 100): BEGIN
	   if max(pixels(0:j,k)) lt max(pixels(j:columns-1,k)) and $
	      max(pixels(j,k:columns-1)) lt max(pixels(j,100:k)) then begin
		 edge_pixels(j,k)=0
	   endif
	   END	

	(k gt 100) and (j gt 100): BEGIN
	   if max(pixels(j:columns-1,k)) lt max(pixels(100:j,k)) and $
	      max(pixels(j,k:columns-1)) lt max(pixels(j,100:k)) then begin
		 edge_pixels(j,k)=0
	   endif
          END

	   	
     ENDCASE		
    endif

    endfor
    

  endfor

 edge_pixels=median(edge_pixels,2)

 edge_pixels(0:30,*)=0
 edge_pixels(170:columns-1,*)=0	
 edge_pixels(*,0:30)=0
 edge_pixels(*,170:rows-1)=0

; Fit the image to an ellipsoid:

  loce = where(edge_pixels eq cap_color,counte)

  if counte gt 0 then begin

    x = loce mod 200
    y = loce/200
    fit_ellipse, x, y, a, b, xo,yo, percent

  endif else print, 'No edge found'

  column=column(1:count-1)
  row=row(1:count-1)

 
 for j=0,columns-1 do begin;

    for k=0,rows-1 do begin

	   if (edge_pixels(j,k) eq 0) then begin
               	if ((max(edge_pixels(0:j,k)) eq cap_color  and  $
                   max(edge_pixels(j:columns-1,k)) eq cap_color)  or  $
                   (max(edge_pixels(j,0:k)) eq cap_color  and  $
                   max(edge_pixels(j,k:rows-1))) eq cap_color) then begin
	   	   cap_pixels(j,k)=cap_color
	        endif
           endif
     endfor

 endfor


  cap_pixels=smooth(cap_pixels,4) 


  loc = where(cap_pixels gt 0 and pixels lt ibound, cap)
  if(boundaries eq 2 or boundaries eq 3) and (cap gt 0) then 		$
    cimage(loc) = 0.75*cap_color

  loc = where(edge_pixels eq cap_color, edge)
  if (boundaries eq 1 or boundaries eq 3) and (edge gt 0) then 		$
    cimage(loc) = cap_color

;  for j=0,columns-1 do begin;
;    for k=0,rows-1 do begin
;      	  if cap_pixels(j,k) gt 0 and pixels(j,k) lt ibound then begin
;		  cap=cap+1
;		  if boundaries eq 2 then cimage(j,k) = 0.75*cap_color
;		  if boundaries eq 3 then cimage(j,k) = 0.75*cap_color
;          endif 
;       	  if edge_pixels(j,k) eq cap_color then begin
;		  if boundaries eq 1 then cimage(j,k)  =cap_color
;		  if boundaries eq 3 then cimage(j,k)  =cap_color
;	  endif
;     endfor
;  endfor

;  print,'  Area of polar cap=',cap,' Pixels'

;-----------------------------------------------------------------------------


; convert pixels to km^2 :

  cap=cap*1970.0

  pn = i mod ppp
  pagen = i/ppp + 1

  if (pagen gt 1) and (pn eq 0) then begin
    if (moreppp) then begin
      device, /close
      outfile = psfile+'_'+chopr('0'+tostr(pagen),2)+'.ps'
      setdevice,outfile,'land',4,0.95
    endif
    if (!d.name eq 'X') then prompt_for_next 
  endif

  if (pn eq 0) then begin
    plotdumb
    itime = intarr(6)
    itime(0) = strmid(uttimes(i),0,2)
    itime(1) = strmid(uttimes(i),2,2)
    itime(2) = strmid(uttimes(i),4,2)
    c_a_to_s, itime, stime
    date = strmid(stime,0,9)
    xyouts, 0.5, 1.0, date, /norm, alignment = 0.5
    ctpos = [0.98,0.3,1.0,0.7]
    plotct, ncolors, ctpos, maxct, 'photons/(cm2 s)', /right
  endif

  get_position, ppp, space, sizes, pn, pos  

  xs = pos(2) - pos(0)
  ys = pos(3) - pos(1)

  plot, [-maxran,maxran],[-maxran,maxran],		$
	pos = pos, xstyle = 5, ystyle = 5,		$
	/nodata, /noerase

  image = float(ncolors)*cimage/float(maxct)
  tv, image, pos(0), pos(1), xsize = xs, ysize = ys, /norm  

  if lun gt 0 then 				$
    printf,lun,format="(6E12.4,2X,A13)",		$
           a,b,xo,yo,percent,cap,uttimes(i)

  ca=strarr(1)
  pix=strarr(1)
  ca(0)=cap
  pix(0)=ibound
  hms = strmid(uttimes(i),7,4)+':'+strmid(uttimes(i),11,2)+' UT'
  ca= 'Polar Cap : '+ ca(0)+ ' km^2'
  pix= 'OVAL MIN PIXEL = ' + pix(0)+ ' Rlys' 

  npts = 200

; convert pixels to degrees:
  ptod = 100.0/maxran

  a = a/ptod
  b = b/ptod
  xo = xo/ptod
  yo = yo/ptod
  r = (xo^2.0 + yo^2.0)^0.5
;  print, '  Polar cap offset by ',r,' degrees'
;  print, '  Area of ellipse = ',!pi*float(a)*float(b),' degrees'
  ae = !pi*float(a)*float(b)*111.0^2.0
  aes = 'Ellipse : '+strcompress(string(ae),/remove_all)+' km^2'

  xyouts, -maxran, 0.9*maxran, hms, color = 255, charsize=0.9
  xyouts, -maxran, -1.1*maxran,ca, charsize = 0.8
  xyouts, -maxran, -1.2*maxran,aes, charsize = 0.8

  xyouts, -maxran, -0.90*maxran,					$
	'a:'+tostrf(a*111.0),color=255,charsize=0.75
  xyouts, -maxran, -0.98*maxran,					$
	'b:'+tostrf(b*111.0),color=255,charsize=0.75

  xyouts, maxran, -0.90*maxran,						$
	'x:'+tostrf(xo*111.0),color=255,charsize=0.75, alignment=1.0
  xyouts, maxran, -0.98*maxran,						$
	'y:'+tostrf(yo*111.0),color=255,charsize=0.75, alignment=1.0

  fac = 2.0-0.00001
  xb = (findgen(npts/2+1)^2.0/float(npts/2)^2.0-1.0)*fac*a/2.0
  xt = [xb,-reverse(xb(0:npts/2-1))]+xo
  yt = ((1.0 - ((xt-xo)/a)^2.0)^0.5)*b + yo
  oplot, xt, yt, color = 255
  yt = -((1.0 - ((xt-xo)/a)^2.0)^0.5)*b + yo
  oplot, xt, yt, color = 255

  plotmlt, maxran, /white, /no00, /no12

  print, '  Area of ellipse   = ',!pi*float(a)*float(b)*111.0^2.0,' km^2'
  print, '  Area of polar cap = ',cap,' km^2'
 
endfor

if lun gt 0 then close,lun

end
