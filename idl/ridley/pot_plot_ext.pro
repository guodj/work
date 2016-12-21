pro newprompt_for_next
print, 'Next frame or quit (return/q)?'
que=''
read, que
if que eq 'q' then stop
end



pro pot_plot_ext

  filein = ''
  fileti = ''
  amiefile = ''
  fileout = ''
  line = ''

  read, 'Enter file name of the AMIE output : ', amiefile
;  filein = amiefile+'.dat'
;  fileti = amiefile+'.time'

  close,1
  openr,1,amiefile

  readf,1, line
  lons = fix(strmid(line,40,3))
  lats = fix(strmid(line,36,3))+1
  count = fix(strmid(line,23,4))

  ltpos = fltarr(lons,lats)
  lnpos = fltarr(lons,lats)
  lts = fltarr(lats)
  fakedata = fltarr(1,lons)
  date = strarr(count)
  time = strarr(count)
  data = fltarr(count,lons,lats)
  mon = 'JanFebMarAprMayJunJulAugSepOctNovDec'

  done = 0
  while not done do begin

    readf,1, line
    if strpos(line,"ntime iyear") gt 0 then done = 1

  endwhile

  for i=0,count-1 do begin

    year = strmid(line,5,4)
    month = fix(strmid(line,10,2))
    date(i) = strmid(mon,(month-1)*3,3)+' '+strmid(line,13,2)+', '+year
    h = fix(strmid(line,16,2))
    m = fix(strmid(line,19,2))
    time(i) = chopr('0'+tostr(h),2)+':'+chopr('0'+tostr(m),2)+' UT'

    readf,1, line
    readf,1, line

    for j=0,lats-1 do begin

      readf,1, line
      lts(j) = float(strmid(line,5,5))
      readf,1, format='(24E12.4)',fakedata
      data(i,0:lons-1,j) = fakedata(0,0:lons-1)

    endfor

    if i lt count-1 then readf,1, line

  endfor

  asub = 0
  sub = 0

  que$=''
  read,'Subtraction to be done (y/n) ? ', que$

  if (strmid(que$,0,1) eq 'y') or (strmid(que$,0,1) eq 'Y') then begin
    que$ = ''
    read,'Subtract an averaged pattern or running average rather than subtract last (a/r) ? ', que$

    if (strmid(que$,0,1) eq 'a') or (strmid(que$,0,1) eq 'A') then begin
      levfac=3.0
      n1 = 0
      n2 = 0               
      print, 'Enter averaging start n, stop n :'
      read, n1,n2

      data(0,*,*) = data(n1,*,*)
      for n=n1+1,n2 do begin
        data(0,*,*) = data(0,*,*) + data(n,*,*)
      endfor
      data(0,*,*) = data(0,*,*)/(n2-n1+1)

      for timeindex= 1,count-1 do $
          data(timeindex,*,*) = data(timeindex,*,*) - data(0,*,*)
    end else if (strmid(que$,0,1) eq 'r') or (strmid(que$,0,1) eq 'R') then begin

      levfac=3.0
      print, 'Enter halfwidth of filter (nr of frames): '
      nmin = 0
      read, nmin
      tempdata=data
      for timeindex= nmin,count-1-nmin do $
        tempdata(timeindex,*,*) = $
        total(data(timeindex-nmin:timeindex+nmin,*,*),1)/(2.0*nmin+1.0)
      data=data-tempdata
    end else begin

      for timeindex= count-1,1,-1 do $
        data(timeindex,*,*) = data(timeindex,*,*) - data(timeindex-1,*,*)
      levfac=1.0
    end
  end else levfac = 5.0


    lev = findgen(30)*levfac

  read, 'Enter ps file name (return for screen) : ', fileout

  ppp = 9
  print, 'Enter number of plots per page : '
  read, ppp

  maxran = max(lts)

  maxcolat=''
  print, ' Enter max colatitude for plot (default: ', maxran,'):'
  read, maxcolat
  if fix(maxcolat) ne 0 then maxran = maxcolat

  for i=1,lons do for j=1,lats do begin
    ltpos(i-1,j-1) = lts(j-1)
    lnpos(i-1,j-1) = 180.0-float(i)*360.0/float(lons)
  endfor

  loc = where(lnpos lt 0.0, count)

  if count gt 0 then lnpos(loc) = lnpos(loc) + 360.0

  condata = fltarr(3,lats*lons)
  psfile = fileout
  gifim = 0
  if strlen(psfile) gt 0 then begin
    if strpos(psfile,'.ps') gt -1 then begin
      setdevice, psfile, 'l', 4, 0.95
      plot, [0,1], /nodata, xstyle=5, ystyle=5
    endif else begin
      gifim = 1
      set_plot,'Z'
      device, set_res=[400,600], z=0
    endelse
  endif else begin
    window,1, xsize=800, ysize=600
    plot, [0,1], /nodata, xstyle=5, ystyle=5
  endelse

  space = 0.01
;  ppp = 9
  space = space*9.0/float(ppp)
  pos_space, ppp, space, sizes

  aaa = !pi*findgen(361.0)/180.0
  xp = sin(aaa)
  yp = cos(aaa)
  ; maxran = 40.0

  for i=0,lons-1 do for j=0,lats-1 do begin
    condata(0,i*lats+j) = lnpos(i,j)
    condata(1,i*lats+j) = ltpos(i,j)
  endfor


  print, 'Enter starting, stopping number to display (-1,-1 for individual) :'
  read,n1,n2

  if (n1 eq -1) then begin

    print, 'you have selected to enter individual numbers.'
    print, ''

    done = 0

    while (done eq 0) do begin

      que = ''
      read, 'Enter number (return to quit) : ',que
      if (strlen(que) eq 0) then done = 1			$
      else begin

        if (n_elements(list) eq 0) then list = fix(que)		$
	else list = [list,fix(que)]

      endelse

    endwhile

    n1 = 0
    n2 = n_elements(list)-1

  endif

  print, 'Enter figure heading:'
  figure=''
  read,figure

  post = fltarr(ppp,4)

  get_position,ppp,space,sizes,0,pos
  ys = 1.0-pos(3)
  for i=0,ppp-1 do begin
    get_position,ppp,space,sizes,i,pos
    pos([1,3]) = pos([1,3]) + ys
    post(i,*) = pos
  endfor

;  if strlen(fileout) eq 0 then potout = amiefile+'.pot'		$
;  else potout = strmid(fileout,0,strpos(fileout,'.ps'))+'.pot'
;  if (sub eq 1) then openw,1,potout
 

  n1sub = n1-sub
  if asub eq 0 then n1sub = n1
  if n1sub eq 0 and asub eq 0 and sub eq 1 then n1sub = 1

  filen = 0

  for ij=n1sub,n2 do begin

    if (ij eq n1-1) then n = 0 else n = ij
    if asub eq 1 then subn = 0 else subn = ij-1

    if asub eq 0 then pn = (n-n1) mod ppp			$
    else if ij ge n1 then pn = (ij-(n1-sub)+2) mod ppp else pn = 1

    pos = post(pn,*)

    if (n_elements(list) gt 0) and (ij ge 0) then n = list(ij)

    if ((pn eq 0) or (ij eq n1-1)) and (gifim eq 1) then begin
          set_plot,'Z'
          device, set_res=[400,600], z=0
          plot, [0,1], /nodata, xstyle=5, ystyle=5
    endif

    if ij eq n1-1 then begin

      plot, [0,1], /nodata, /noerase, pos = [0,0,1,1], 		$
	xstyle=5, ystyle=5
      if ppp eq 9 then begin
        xyouts, pos(0)-0.04, pos(3)-0.015, '(A)', charsize = 1.2
        xyouts, pos(0) - (pos(2)-pos(0)+space), pos(1)-space/2.0+0.005,'(B)', $
	  charsize=1.2
	        oplot, [pos(0)-(pos(2)-pos(0)+space)+0.05,		$
		pos(2)+(pos(2)-pos(0)+space)], 			$
	       [pos(1)-space/2.0,pos(1)-space/2.0]
      endif
    endif

    plot, maxran*xp, maxran*yp, xstyle=5, ystyle=5, pos=pos, /noerase

    oplot, maxran*xp, maxran*yp
    oplot, 30.0*xp, 30.0*yp,linestyle=1
    oplot, 20.0*xp, 20.0*yp,linestyle=1
    oplot, 10.0*xp, 10.0*yp,linestyle=1
    oplot, [0.0,0.0],[-maxran,maxran]
    oplot, [-maxran,maxran], [0.0,0.0]

    if pn eq 1 then 						$
    xyouts, 0.5, 1.03, date(0), alignment = 0.5, charsize=1.2, /norm

    if ij ne n1-1 then begin
      xyouts, -1.0*maxran, 0.90*maxran, time(n), 		$
        alignment = 0.0, charsize=0.9
    endif else begin
      xyouts,-0.95*maxran, 0.90*maxran, 'Steady', 		$
	alignment = 0.0, charsize=0.9	
    endelse

    if (sub eq 1) and (ij ne n1-1) then begin
      xyouts, 0.95*maxran, 0.90*maxran, 'Residual',  	$
	alignment = 1.0, charsize=0.7
    endif

    tp = sizes.nby*sizes.nbx
    if (n2-ij lt sizes.nbx) or (pn+1 gt tp-sizes.nbx) then       $
      xyouts, 0.0,maxran*(-1.1), '00', alignment=0.5, 		$
	charsize=0.8
    if (pn mod sizes.nbx eq sizes.nbx-1) or (n2-ij eq 0) then   $
      xyouts, maxran*(1.01), -0.05*maxran, '06', 		$
	charsize=0.8
    if (pn lt sizes.nbx) then                 $
      xyouts, 0.0,maxran*(1.02), '12', alignment=0.5, 		$
	charsize=0.8
    if (pn mod sizes.nbx eq 0) or (ij eq n1) then                $
      xyouts, maxran*(-1.01), -0.05*maxran, '18', 		$
	alignment=1.0, charsize=0.8

    for i=0,lons-1 do for j=0,lats-1 do begin
      if (sub eq 1) and (ij ne n1-1) then 				$
        condata(2,i*lats+j) = 						$
		(data(n,i,j)-data(subn,i,j))/1000.0			$
      else								$
        condata(2,i*lats+j) = data(n,i,j)/1000.0
    endfor

    levels = lev + float(fix(min(condata(2,*))/levfac))*levfac - levfac

    basic_contour_new, condata, maxran, pos, levels, /plotmaxmin

    if max(condata(2,*)) lt 100 then begin
      xyouts,maxran*(1.0),maxran*(-0.96), 				$
	  'Ma: '+tostr(max(condata(2,*))), alignment=1.0, 		$
		charsize=0.85
    endif else begin
      xyouts,maxran*(1.0),maxran*(-0.96), 				$
	  'Ma:'+tostr(max(condata(2,*))), alignment=1.0, 		$
		charsize=0.85
    endelse

    if min(condata(2,*)) gt -100 then begin
      xyouts,maxran*(-1.0),maxran*(-0.96),	 			$
	'Mi: '+tostr(min(condata(2,*))), alignment=0.0, 		$
		charsize=0.85
    endif else begin
      xyouts,maxran*(-1.005),maxran*(-0.96),	 			$
	'Mi:'+tostr(min(condata(2,*))), alignment=0.0, 			$
		charsize=0.85
    endelse

;    if ((sub eq 1) and (ij gt n1-sub)) then 				$
;      printf,1,min(condata(2,*)), max(condata(2,*))

    if (pn+1 eq ppp) or (ij eq n2) then begin
      xyouts, 0.5, -0.05, figure, /norm, alignment = 0.5
    endif

    if pn+1 eq ppp then begin
      if !d.name eq 'X' then newprompt_for_next
      if !d.name ne 'Z' then plot, [0,1], /nodata, xstyle=5, ystyle=5
      if gifim eq 1 then begin
        bytemap = tvrd()
        tvlct,rr,gg,bb,/get
        giffile = strmid(psfile,0,strlen(psfile)-4)+              $
                chopr('00'+tostr(filen+1),3)+'.gif'
        write_gif,giffile,bytemap,rr,gg,bb
        device, /close
        set_plot,'X'
	filen = filen+1
      endif
    endif

  endfor

;  if sub eq 1 then close,1

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

    if !d.name eq 'Z' then begin
      bytemap = tvrd()
      tvlct,rr,gg,bb,/get
      giffile = strmid(psfile,0,strlen(psfile)-4)+              $
                chopr('00'+tostr(filen+1),3)+'.gif'
      write_gif,giffile,bytemap,rr,gg,bb
      device, /close
      set_plot,'X'
    endif

  return

end

