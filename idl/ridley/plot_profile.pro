

  line = ''
  fileout = ''
  filein = '' 

  read, 'Enter output ps file name (return for screen) : ', fileout

  ppp = 9
  print, 'Enter number of plots per page : '
  read, ppp


  if strlen(fileout) gt 0 then begin
     if ppp eq 9 then setdevice, fileout,'p',4,0.83		$
     else setdevice, fileout,'l',4,0.95
  endif else window,1,xsize=800,ysize=600

  plot, [0,1], /nodata, xstyle=5, ystyle=5
  im=1 
  flist = ''
    print, 'Enter file list for images : '
    read, flist
    print, 'Enter number of images (and patterns) to skip in between each:'
    nskip = 0
    read, nskip
    nskip = nskip
    maxran=40.0 
;------------------------------------------
    read_images, flist, nskip, uttimes, inttimes, image_save, 90.0-maxran
;------------------------------------------



nimage = n_elements(image_save(*,0,0))
rows=n_elements(image_save(0,0,*))
columns=n_elements(image_save(0,*,0))
lat=findgen(201)
lat(101:200)=lat(1:100)
lat(0:100)=-reverse(lat(0:100))
lat=lat*maxran/100.0
image=fltarr(columns,rows)
space = 0.06
print,'nimage=',nimage,'rows=',rows,'columns=',columns
pos_space, ppp, space, sizes

!P.MULTI=[0,ppp,2*ppp]


for i= 0, nimage-1 do begin

  pn = i mod ppp



  pagen = i/ppp + 1
  moreppp=1

  if (pagen gt 1) and (pn eq 0) then begin
    if (moreppp) then begin
      device, /close
      outfile = fileout+'_'+chopr('0'+tostr(pagen),2)+'.ps'
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
  endif

  get_position, ppp, space, sizes, pn, pos  

  xs = pos(2) - pos(0)
  ys = pos(3) - pos(1)

  plot, [-maxran,maxran],[-maxran,maxran],		$
	pos = pos, xstyle = 5, ystyle = 5,		$
	/nodata, /noerase

  image=fltarr(columns,rows)

  for j=0,columns-1 do begin
	image(*,j)=float(maxran)*image_save(i,*,j)/float(maxran)
  endfor


;  image=roberts(image)
  image=smooth(image,9)
  row=image(*,100)
  column=image(100,*)


!P.REGION=[pos(0),pos(1),pos(2),pos(1)+0.5*ys]
  PLOT,lat,column,XTITLE='Latitude',YTITLE='photons/(cm2 s)', $
     title='Midnight   -   Noon    ' + strmid(uttimes(i),7,4)+':' $
     +strmid(uttimes(i),11,2)+' UT'
; , yrange=[0,1200]

!P.REGION=[pos(0),pos(3)-0.5*ys,pos(2),pos(3)]
  PLOT,lat,row,XTITLE='Latitude',YTITLE='photons/(cm2 s)', $  
     title='Dusk   -   Dawn    '+strmid(uttimes(i),7,4)+':' $
     +strmid(uttimes(i),11,2)+' UT'
; , yrange=[0,1200]



endfor



;------------------------------------------
;    if im eq 1 then begin
;      xyouts, 1.0*maxran, 0.85*maxran, 				$
;	strmid(uttimes((ij-n1sub)/nskip),7,4)+			$
;	':'+strmid(uttimes((ij-n1sub)/nskip),11,2), 		$
;        alignment = 1.0, charsize=0.75, color = color
;    endif
;------------------------------------------

;------------------------------------------



  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif



end

