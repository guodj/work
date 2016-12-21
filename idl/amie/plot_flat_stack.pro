
function tostr,value
  return, strcompress(string(fix(value)),/remove_all)
end

function chopr, svalue, n
  if strlen(svalue) lt n then n = strlen(svalue)
  return, strmid(svalue, strlen(svalue)-n,n)
end

function nmean, data_array

  ave = 0.0
  for i=long(0),n_elements(data_array)-1 do ave=ave+data_array(i)

  ave = ave/float(n_elements(data_array))

  return, ave

end

common ffinfo, header

dename = !d.name

dirin = ''
ext = ''
filein = ''
linein = ''
line = ''

smax = 99900.

filein = ''
read, 'Enter files to input ("filein".hed, filein can be *) : ', filein

filelist = findfile('-1 '+filein+'.hed')

nf = n_elements(filelist)
print, 'Number of files found : ',n_elements(filelist)

sdate = ''

ndays = 0.0   

filein = strmid(filelist(0),0,strpos(filelist(0),'.hed'))
setup_var_menu, 1, filein, timeint, vars, units, count, nrows, rowlen

print, 'Start time of file #001 : ',timeint(0)
print, 'End time of file #001   : ',timeint(1)
sdate = ask('Enter start date in form "dd-mon-yy" : ',timeint(0))
ndays = ask('Enter number of days to plot : ','1')

c_s_to_a, sitime, sdate
c_a_to_r, sitime, stime

etime = stime + float(ndays)*24.0*3600.0-60.0
c_r_to_a, eitime, etime

c_r_to_a, itime, stime

psfile = ''
print, 'Enter postscript file name (return to plot to screen) :'
read, psfile

ppp = 0
print, 'Enter number of magnetometers per page to plot : '
read, ppp

data = fltarr(nf,3,long(1440)*ndays*60) + 99999.0
time = dblarr(nf,3,long(1440)*ndays*60)
stats = strarr(nf,3)

for i=0,n_elements(filelist)-1 do begin
  
  filein = strmid(filelist(i),0,strpos(filelist(i),'.hed'))
  print,'Reading file # ',chopr('00'+tostr(i+1),3),': ',filein

  setup_var_menu, 1, filein, timeint, vars, units, count, nrows, rowlen
  for j=0,2 do stats(i,j) = strcompress(vars(j), /remove_all)
  
  id = 0
  header = {nr : [nrows(0)], unit : [51], rl : [rowlen(0)]}
  
  openr, header.unit(0), filein+'.dat'
  
  srow = loc_ffunix_row(id, stime)
  erow = loc_ffunix_row(id, etime)
  
  if erow-srow gt 0 then begin
  
    for row=srow,erow do begin
  
      fdget, id, row, uttime, indata
      for j=0,2 do begin
        if indata(j) eq -1.0e32 then data(i,j,row-srow) = 999999.0	$
        else data(i,j,row-srow) = indata(j)
      endfor
      time(i,0:2,row-srow) = uttime(0)
  
    endfor
  
  endif

  close, header.unit(0)

endfor

if strlen(psfile) gt 0 then setdevice, psfile, 'l', 4, 0.95	$
else window,1, xsize=800, ysize=600
  
title = ''

time = time - stime

time_axis, sitime, eitime, btr, etr,  $
	xtickname, xtitle, xtickvalue, xminor, xtickn

space = 0.02
xmin = 0.02
xmax = 0.99
ymin = 0.00-space
ymax = 1.00

if ppp eq 1 then begin

  filen = 0

  while filen le nf-1 do begin

    range = 0.0
    ii = filen

    dy = (ymax-ymin)/3.0

    for icomp = 0, 2 do begin

      loc = where(data(ii,icomp,*) lt 99999.0,count)
      if count gt 0 then begin
        d = reform(data(ii,icomp,loc))
        maxd = max(d)
        mind = min(d)
      endif else begin
        mind = 0.0
        maxd = 0.0
      endelse

      if (maxd - mind gt range) then begin
        range = maxd - mind
      endif

    endfor

    range = float(fix((1.33*range+99.0)/100.0))*100.0
    if (range lt 100.0) then range = 100.0

    for icomp = 0, 2 do begin

      if (icomp eq 0) then plotdumb

      pos = [xmin,ymax-float(icomp+1)*dy+space,xmax, ymax-float(icomp)*dy]

      if (icomp eq 2) then begin
        xtickna = xtickname
        xtitl  = xtitle
      endif else begin
        xtickna = strarr(n_elements(xtickname)) + ' '
        xtitl  = ' '
      endelse

      d = reform(data(ii,icomp,*))
      t = reform(time(ii,icomp,*))
      loc = where(d lt 99999.0,count)
      if count gt 0 then begin
        maxd = max(d(loc))
        mind = min(d(loc))
        m = nmean(d(loc))
      endif else begin
	d = fltarr(n_elements(time))+99999.0
        mind = 0.0
        maxd = 0.0
	m = 0.0
      endelse
      rs = range - (maxd-mind)
      m1 = maxd + rs/2.0
      m2 = mind - rs/2.0

      plot, [btr,etr],[m1,m2],                             $
        xstyle = 1,                                             $
        ystyle = 1,                                             $
        pos = pos,                                              $
        /nodata,                                                $
        /noerase,                                               $
        xtickname = xtickna,                                    $
        xtitle = xtitl,                                         $
        xtickv = xtickvalue,                                    $
        xminor = xminor,                                        $
        xticks = xtickn,                                        $
        title = title

      oplot, t, d, max_value = 99998.0

;      for i=long(0),etr,1440.0*60.0 do 				$
;	oplot,[i,i],[m1,m2],linestyle=2

      oplot, [btr,etr],[m,m], linestyle = 2

      xyouts, [etr+etr/100],[(m1+m2)/2.0],                    	$
          strmid(stats(ii,icomp),3,1) + ' mean : ' +            $
          tostrf(m) + ' nT',                    		$
          alignment=0.5, orientation=270

      if (icomp eq 0) then begin
          xyouts, [btr],[(m1+m2)/2.0 + range*0.55],             $
            strmid(stats(ii,icomp),0,3), charsize = 1.5
          xyouts, [1.03],[0.5],                 		$
            'Plot Scale : '+tostrf(range)+' nT',                $
            /norm,                                              $
            alignment=0.5, orientation=270, charsize = 0.9
      endif

      for i=1,ndays-1 do begin
        t = btr+float(i)*24.0*3600.0
        oplot, [t,t],[m1,m2], linestyle = 1
      endfor

    endfor

    filen = filen + 1
    if ((filen le nf-1) and (!d.name ne 'PS')) then prompt_for_next 

  endwhile

endif else begin

  range = 0

  for ii=0,nf-1 do begin

    for icomp = 0, 2 do begin

      loc = where(data(ii,icomp,*) lt 99999.0,count)
      if count gt 0 then begin
        d = data(ii,icomp,loc)
        maxd = max(d)
        mind = min(d)
      endif else begin
        mind = 0.0
        maxd = 0.0
      endelse

      if (maxd - mind gt range) then begin
        range = maxd - mind
      endif

    endfor

  endfor

  range = float(fix((1.33*range+99.0)/100.0))*100.0
  if (range lt 100.0) then range = 100.0
  if (range gt 1000.0) then range = 1000.0

  range = float(ask("offset between stack plots",tostr(fix(range))))

  offset = range

  ytickname = strarr(ppp+1)+' '
  ytickvalue = findgen(ppp+1)*offset
  yminor = 10
  ytickn = ppp
  xspa = space
  dx = (xmax-xmin)/float(3)

  filen = 0

  while filen le nf-1 do begin

    for icomp = 0, 2 do begin

      sfile = filen
      efile = filen + ppp - 1
      if efile ge nf-1 then efile = nf-1

      if (icomp eq 0) then plotdumb

      pos = [xmin+float(icomp)*dx,ymin,xmin+float(icomp+1)*(dx)-xspa, ymax]

      plot, [btr,etr], [0,offset*float(ppp+1)], $
        xstyle = 1,                                     $
        ystyle = 1,                                     $
        pos = pos,                                      $
        /nodata,                                        $
        /noerase,                                       $
        xtickname = xtickname,                          $
        xtitle = xtitle,                                $
        xtickv = xtickvalue,                            $
        xminor = xminor,                                $
        xticks = xtickn,                                $
        ytickname = ytickname,                          $
        ytickv = ytickvalue,                            $
        yminor = yminor,                                $
        title = title,                                  $
        yticks = ytickn

      for i=long(0),max(time),1440.0*60.0 do 		$
	oplot,[i,i],[0,offset*float(ppp+1)],linestyle=2

      for ii = sfile, efile do begin

        offp = offset*float(ppp-(ii-sfile))
        d = data(ii,icomp,*)
        t = time(ii,icomp,*)
	loc = where(d lt 99998.0,count)
	if count gt 0 then begin
          m = nmean(d(loc)) 
          d(loc) = d(loc) - m + offp
        endif else begin
          m = nmean(d)
	endelse
        oplot, t, d, max_value=99998.0
        xyouts, [etr+etr/100],[offp], tostrf(m),              	$
            alignment=0.5, orientation=270, charsize = 0.7
        if (icomp eq 0) then begin
          xyouts, [btr-etr/20],[offp],                          	$
            strmid(stats(ii,icomp),0,3)+strmid(stats(ii,icomp),4,1),	$
            alignment=0.5, orientation=90, charsize=0.8
        endif
        xyouts, [btr+etr/50],[offp+offset/5.0],                 $
          strmid(stats(ii,icomp),3,1)

      endfor

      for i = sfile, sfile+ppp-1 do begin
        offp = offset*float(ppp-(i-sfile))
        oplot, [btr,etr],[offp, offp],                          $
          linestyle=2
      endfor

      if icomp eq 2 then begin
        xyouts, 1.0, 0.5,                                       $
          'Division between baselines is '+tostrf(offset)+' nT',        $
          /norm, alignment = 0.5, orientation = 270.0
      endif

    endfor

    filen = efile + 1

    if (!d.name eq 'X' and filen le nf-1) then prompt_for_next

  endwhile

endelse

if strlen(psfile) gt 0 then begin
  device, /close
  set_plot, dename
endif


end
