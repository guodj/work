
print, 'Enter file list which you would like plotted : '
flist = '/l2/software/idl/polar/geog/JO/970410/short/970410_21*.bin '
read,flist

maxran=40.0

print, 'Enter number of times to skip in between each image : '
print, '(0 for no skipping, 1 for every other image, ...)'
nskip = 0
read, nskip
ntimes = nskip + 1

print, 'Enter color table which you would like :'
ct = 'polar2.ct'
read,ct

print, 'Enter maximum reyleighs for color table (200-600) :'
maxct = 600
read,maxct

print, 'Enter minimum reyleighs for nightime auroral oval (120) :'
nightbound=120
read,nightbound

;print, 'Enter minimum reyleighs for daytimetime auroral oval (300) :'
daybound=300
;read,daybound

print, ''
print, '1. Auroral Oval boundaries only'
print, '2. Auroral Oval colored'
print, '3. Auroral Oval boundaries and Polar Cap colored'
print, '4. Auroral Oval and Polar Cap colored'
print, 'Enter areas you would like colored : '
boundaries = 4
read, boundaries

oval_color=maxct
cap_color=maxct
if (boundaries eq 4) then begin
	oval_color=maxct
	cap_color=round(0.5*maxct)
endif

print, 'Enter number of plots per page : '
ppp = 6
read, ppp

print, 'Enter ps file name (return for screen) :'
psfile = 'test1.ps'
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

readct, ncolors, ct

read_images, flist, nskip, uttimes, inttimes, image_save, 90-maxran

nimage = n_elements(image_save(*,0,0))

space = 0.06

pos_space, ppp, space, sizes


rows=n_elements(image_save(0,0,*))
columns=n_elements(image_save(0,*,0))
pixels=fltarr(columns,rows)


for i= 0, nimage-1 do begin


  cap=0.0
  oval=0.0


  for j=0,columns-1 do begin
	pixels(*,j)=float(maxran)*image_save(i,*,j)/float(maxran)
  endfor

  pixels=smooth(pixels,2)

;-----------------------------------------------------------------------------

;  This section scans the image by column to integrate the size of the polar
;  cap.

;-----------------------------------------------------------------------------

  for j=0,columns-1 do begin
    flag=0 
    bound=nightbound

  row=pixels(*,j)
  column=pixels(j,*)

    for k=0,rows-1 do begin

	if (column(k) gt bound) and (flag eq 0) then begin
	    flag=1
            image_save(i,j,k)=oval_color
	endif

	if (column(k) lt bound) and (flag eq 1) and (k lt 100) then begin
	    flag=2
	    cap=cap+1
            image_save(i,j,k)=cap_color	    
	endif	

	if (column(k) lt bound) and (flag eq 2) and $
                (max(column(k:rows-1)) gt bound) then begin
	    cap=cap + 1
            if (boundaries eq 3) then image_save(i,j,k)=cap_color
            if (boundaries eq 4) then image_save(i,j,k)=cap_color
	endif

	if (column(k) gt bound) and (flag eq 2) then begin
            image_save(i,j,k)=oval_color
	    flag=3
	endif

	if (column(k) lt bound) and (flag eq 3) then begin
            image_save(i,j,k)=oval_color
	    flag=4
	endif

	if (column(k) lt bound) and (flag eq 4) and  $
  	     (k lt 100) and (max(column(k:rows-1)) gt bound) then begin
	   flag=2
	   cap=cap+1 
           image_save(i,j,k)=cap_color
	endif
          
	if (column(k) lt bound) and (flag eq 4) and  $
  	     (k gt 100) and (max(column(k:rows-1)) gt bound) then begin
	   bound=daybound
	   flag=2
           image_save(i,j,k)=oval_color
	endif
          

    endfor
    

  endfor

;-----------------------------------------------------------------------------
;  This section scans the image by column to integrate the size of the auroral
;  oval.
;-----------------------------------------------------------------------------


  for j=0,columns-1 do begin

  row=pixels(*,j)
  column=pixels(j,*)

    for k=0,rows-1 do begin

	if (column(k) gt bound) then begin
	    oval=oval+1
            if (boundaries eq 2) then image_save(i,j,k)=oval_color
            if (boundaries eq 4) then image_save(i,j,k)=oval_color
	endif

    endfor
  endfor	


;-----------------------------------------------------------------------------

; This section of the code switches to scanning by row to better identify the
; edges of the auroral oval 

;-----------------------------------------------------------------------------


  for k=0,rows-1 do begin
    flag=0 
    bound=nightbound

  row=pixels(*,k)
  column=pixels(k,*)


    for j=0,columns-1 do begin


	if (row(j) gt bound) and (flag eq 0) then begin
	    flag=1
            image_save(i,j,k)=oval_color
	endif

	if (row(j) lt bound) and (flag eq 1) then begin
	    flag=2
            image_save(i,j,k)=oval_color	    
	endif	

	if (row(j) gt bound) and (flag eq 2) then begin
            image_save(i,j,k)=oval_color
	    flag=3
	endif

	if (row(j) lt bound) and (flag eq 3) then begin
            image_save(i,j,k)=oval_color
	    flag=4
	endif

	if (row(j) lt bound) and (flag eq 4) and  $
  	     (k lt 100) and (max(row(j:columns-1)) gt bound) then begin
	   flag=2
           image_save(i,j,k)=oval_color
	endif
          
	if (row(j) lt bound) and (flag eq 4) and  $
  	     (k gt 100) and (max(row(j:columns-1)) gt bound) then begin
	   bound=daybound
	   flag=2
           image_save(i,j,k)=oval_color
	endif
          

    endfor
    

  endfor

;-----------------------------------------------------------------------------


  oval=oval*1970.0
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

  image = float(ncolors)*image_save(i,*,*)/float(maxct)
  tv, image, pos(0), pos(1), xsize = xs, ysize = ys, /norm  

  ca=strarr(1)
  ov=strarr(1)
  ca(0)=cap
  ov(0)=oval
  hms = strmid(uttimes(i),7,4)+':'+strmid(uttimes(i),11,2)+' UT'
  ca= 'CAP   = '+ ca(0)+ ' km2'
  ov= 'OVAL = '+ov(0)+' km2'
  xyouts, -maxran, 0.9*maxran, hms, color = 255, charsize=0.8
  xyouts, -maxran, -1.1*maxran,ca, color = 0, charsize=0.8
  xyouts, -maxran, -1.2*maxran,ov, color = 0, charsize=0.8

  tp = sizes.nby*sizes.nbx

  no00 = 1
  if (pn mod sizes.nbx eq sizes.nbx-1) or 			$
     (i eq nimage-1) then no06 = 0 else no06 = 1
  if (pn lt sizes.nbx) then no12 = 0 else no12 = 1
  if (pn mod sizes.nbx eq 0) or 				$
     (i eq 0) then no18 = 0 else no18 = 1

  plotmlt, maxran, /white,		 			$
	no00 = no00, no06 = no06, no12 = no12, no18 = no18

;  plotmlt, maxran, /white

  print, ' '
  print,'image= '+strmid(uttimes(i),7,4)+':'+  $
      strmid(uttimes(i),11,2)+' UT','       oval=',oval,'        cap=',cap
 
endfor

end
