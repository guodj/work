
list = findfile("-t i*.idl")
if strlen(list(0)) gt 0 then filein = list(0) $
else filein = 'in000000.idl'

filein = ask('filename',filein)

list = findfile(filein)
nfiles = n_elements(list)
if nfiles eq 1 and strlen(list(0)) eq 0 then begin
  print, "I can't seem to find that file."
  stop
endif

line = ''

mr = float(ask('maximum range','50'))

time = dblarr(nfiles)
iter = lonarr(nfiles)
itime = intarr(6)

for n = 0, nfiles-1 do begin

  filein = list(n)

  if (strpos(filein,'in') eq 0) and (strpos(filein,'.idl') eq 8) then begin
    iter(n) = long(strmid(filein,2,6))
  endif else iter(n) = 0

  if strpos(filein,'save') gt -1 then begin
    restore, filein
  endif else begin

    openr,1,filein

    done = 0

    while (not done) do begin

      readf,1, line

      if (strpos(mklower(line),"numerical") gt -1) then begin

        readf,1, nvars
        readf,1, nlats
        readf,1, nlons

        tmp = fltarr(nvars)
        data = fltarr(2,nvars,nlons,nlats)

        vars = strarr(nvars)

      endif

      if (strpos(mklower(line),"variable") gt -1) then begin

        for i=0,nvars-1 do begin
          readf,1,line
          vars(i) = strmid(line,6,strlen(line)-6)
        endfor

      endif

      if (strpos(mklower(line),"time") gt -1) then begin

        int_tmp = 0
        for i=0,5 do begin
          readf, 1, int_tmp
          itime(i) = int_tmp
        endfor

        c_a_to_r, itime, rtime
        time(n) = rtime

      endif

      if (strpos(mklower(line),"northern") gt -1) then begin

        for j=0,nlons-1 do for i=0,nlats-1 do begin
          readf,1,tmp
          data(0,*,j,i) = tmp
        endfor

      endif

      if (strpos(mklower(line),"southern") gt -1) then begin

        for j=0,nlons-1 do for i=0,nlats-1 do begin
          readf,1,tmp
          data(1,*,j,i) = tmp
        endfor

      endif

      if eof(1) then done = 1

    endwhile

    close,1

    if (n eq 0) then begin

      nvars_to_plot = 1
      nphi = -1
      for i=0,nvars-1 do if strpos(mklower(vars(i)),'phi') gt -1 then nphi = i
      varlist = [nphi]

      if (nvars_to_plot eq 0) then begin
	print, "Can not continue!"
	stop
      endif

      data_to_plot = fltarr(2,nfiles,nvars_to_plot,nlons,nlats)
      theta = fltarr(2,nlons,nlats)
      phi   = fltarr(2,nlons,nlats)

      nt = -1
      for i=0,nvars-1 do if strpos(mklower(vars(i)),'theta') gt -1 then nt = i
      np = -1
      for i=0,2 do if strpos(mklower(vars(i)),'psi') gt -1 then np = i
;I have it named wrong for a number of runs...
      for i=0,2 do if strpos(mklower(vars(i)),'phi') gt -1 then np = i
      if (nt eq -1) or (np eq -1) then begin
        print, "Can't file Theta or Phi variable. Please check file."
        stop
      endif

      theta(*,*,*) = reform(data(*,nt,*,*))
      phi(*,*,*)   = reform(data(*,np,*,*))

    endif

    data_to_plot(*,n,*,*,*) = reform(data(*,varlist,*,*))

  endelse

endfor


c_r_to_a, itime, time(0)
c_a_to_s, itime, mhd_time
mhd_time = strmid(mhd_time,0,15)


filelist_dmsp = findfile('a*.*')

nfiles_dmsp = n_elements(filelist)

n = 0

  openr,1,filelist_dmsp(n)

  tmp=fltarr(19)

  readf,1,tmp

  satid = tmp(0)
  itime = fltarr(6)
  itime(0) = tmp(1)
  itime(1) = 1
  itime(2) = tmp(2)
  itime(3) = 0
  itime(4) = 0
  itime(5) = 0

  c_a_to_r, itime, basetime

  nlines = fix(tmp(16))

  data = fltarr(13, nlines)

  line = ''
  for i=0,nlines-1 do begin
    readf,1,line
    data(0,i) = float(strmid(line,0,7))
    data(1,i) = float(strmid(line,7,7))
    data(2,i) = float(strmid(line,14,7))
    data(3,i) = float(strmid(line,21,7))
    data(4,i) = float(strmid(line,28,7))
    data(5,i) = float(strmid(line,35,3))
    data(6,i) = float(strmid(line,38,8))
    data(7,i) = float(strmid(line,46,7))
    data(8,i) = float(strmid(line,53,7))
    data(9,i) = float(strmid(line,60,7))
    data(10,i) = float(strmid(line,67,7))
    data(11,i) = float(strmid(line,74,8))
    data(12,i) = float(strmid(line,81,8))
  endfor

  close,1

  time = basetime + reform(data(0,*))

  vy = reform(data(1,*))
  vz = reform(data(2,*))

  pot = reform(data(6,*))
  lat = reform(data(8,*))
  lat_110  = reform(data(7,*))
  lat_save = reform(data(8,*))
  mlt_save = reform(data(9,*))

  mlt = reform(data(9,*))/24.0 * 2.0 * !pi - !pi/2.0

  if mean(lat) lt 0.0 then begin
    dmsp_hem = 1
    lat = 90.0 + lat
    yfac = -1.0
    xfac = 1.0
  endif else begin
    dmsp_hem = 0
    lat = 90.0 - lat
    yfac = 1.0
    xfac = -1.0
  endelse

  x = lat * cos(mlt)
  y = lat * sin(mlt)

  dx = x(1:nlines-1) - x(0:nlines-2)
  dy = y(1:nlines-1) - y(0:nlines-2)

  dxnew = dx
  dynew = dy

  nsmooth = 10

  for i = nsmooth/2, nlines - nsmooth/2 -2 do begin

    dxnew(i) = mean(dx(i-nsmooth/2:i+nsmooth/2-1))
    dynew(i) = mean(dy(i-nsmooth/2:i+nsmooth/2-1))

  endfor

outfile = ask('output file name',filelist_dmsp(0)+'.batsrus')

space = 0.025
ppp = 6

hems = ["Northern Hemisphere", "Southern Hemisphere"]

for hem = 0, 1 do begin

  print, "Working on ",hems(hem)

  if hem eq 0 then begin
    loc = where(reform(theta(0,0,*)) le mr)
    rang = reform(theta(0,*,loc))
    lons = reform(phi(0,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
    dran = (rang(0,2) - rang(0,0))/2.0
  endif else begin
    loc = where(reform(theta(1,0,*)) gt 180.0-mr)
    rang = 180.0-reform(theta(1,*,loc))
    lons = reform(phi(0,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
    dran = (rang(0,0) - rang(0,2))/2.0
  endelse

  dlon = (lons(2,0) - lons(0,0))/2.0

  for i = 0, nvars_to_plot-1 do begin

    if hem eq 0 and i eq 0 then begin
      pn = -1
      setdevice, list(0)+'.ps','p',4
      pos_space, ppp, space, sizes
    endif

    mini = min(data_to_plot(hem,*,i,*,loc))
    maxi = max(data_to_plot(hem,*,i,*,loc))

    if (maxi eq mini) then maxi = mini + 1.0

    if (mini ge 0.0) then begin
      levels = maxi*findgen(9)/8.0 + mini
      readct,ncolors, getenv("IDL_EXTRAS")+"white_red.ct"
    endif else begin
      maxi = max([abs(mini),maxi])
      mini = -maxi
      levels = fltarr(8)
      levels(0:3) = -maxi*(4-findgen(4))/4.0
      levels(4:7) = maxi*(findgen(4)+1)/4.0
      readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
    endelse

    c_levels = (maxi-mini)*findgen(30)/29.0 + mini
    c_colors = (ncolors-1)*findgen(30)/29.0 + 1

    pn = hem

    for n = 0, nfiles - 1 do begin

      if pn eq 0 then plotdumb
      get_position, ppp, space, sizes, pn, pos

      contour, data_to_plot(hem,n,i,*,loc), xpos, ypos,        $
        /follow, nlevels=30, /noerase, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,/cell_fill

      contour, data_to_plot(hem,n,i,*,loc), xpos, ypos,        $
        /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase

      if ((strpos(mklower(vars(varlist(i))),'phi') gt -1) and $
          (hem eq dmsp_hem)) then begin

        mhd_pot = reform(data_to_plot(hem,n,i,*,loc))
        compare_data = fltarr(6,nlines)
	index        = intarr(nlines)
	dmsptime     = dblarr(nlines)
        postemp = fltarr(2,nlines)
	npts = 0

        mhd_pot_output = mhd_pot
	xpos_output = xpos
        ypos_output = ypos

        dx = dxnew
        dy = dynew

        for ii=1, nlines - 2 do begin

          if (lat(ii) lt mr) then begin

            intx = ((mlt(ii) - !pi/2.0 + !pi*2.0) mod (2.0*!pi))/dlon
            if (hem eq 1) then begin
              inty = n_elements(rang(0,*)) - (lat(ii))/dran
            endif else begin
              inty = lat(ii)/dran
            endelse

            compare_data(0,npts) = pot(ii)
            compare_data(4,npts) = vy(ii)
            compare_data(1,npts) = interpolate(mhd_pot, intx, inty)
            compare_data(2,npts) = $
                    distance(90.0-lat(ii+1), mlt(ii+1)*180.0/!pi, $
                             90.0-lat(ii-1), mlt(ii-1)*180.0/!pi)
            compare_data(3,npts) = -2.0*31100e-9*sin(!pi/2-lat(ii)*!pi/180.0)
            index(npts) = ii
            postemp(0,npts) = interpolate(xpos, intx, inty)
            postemp(1,npts) = interpolate(ypos, intx, inty)
            dmsptime(npts) = time(ii)

            if hem eq 1 then begin
              if (inty ge n_elements(rang(0,*))-2) then $
                  compare_data(1,npts)=-999.0
            endif else begin
              if (inty le 1) then $
                  compare_data(1,npts)=-999.0
            endelse

            npts = npts + 1

            if (abs(vy(ii)) lt 10.0) then begin

              dr = sqrt(dx(ii)^2 + dy(i)^2)

              vel_y = 10.0 * vy(ii) * dx(ii) / dr * yfac
              vel_x = 10.0 * vy(ii) * dy(ii) / dr * xfac

              oplot, [x(ii), x(ii) + vel_x], [y(ii), y(ii) + vel_y]

            endif

            if npts eq 1 then begin
               if (x(ii) lt 0.0) then xyouts, x(ii)-1.0, y(ii), 'S', align=1.0
               if (x(ii) gt 0.0) then xyouts, x(ii)+1.0, y(ii), 'S'
            endif

            if (lat(ii) lt lat(ii-1)) then begin
              c_r_to_a, itime, time(ii)
              c_a_to_s, itime, hightime
            endif

          endif

        endfor

      endif

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

      no00 = 0
      no06 = 1
      no12 = 1
      no18 = 1

;      if (pn+1 gt ppp-sizes.nbx) then no00 = 0
      if (pn mod sizes.nbx eq sizes.nbx-1) then no06 = 0
      if (pn lt sizes.nbx) then no12 = 0
      if (pn mod sizes.nbx eq 0) then no18 = 0

;--------------------------------------------------------------
; Draw the MLT grid
;--------------------------------------------------------------

      plotmlt, mr, no00 = no00, no06 = no06, no12 = no12, no18 = no18

      mini_tmp = min(data_to_plot(hem,n,i,*,loc))
      maxi_tmp = max(data_to_plot(hem,n,i,*,loc))

      if (abs(mini_tmp) gt 1000.0) or (abs(maxi_tmp) gt 1000.0) or        $
         (abs(maxi_tmp) lt 0.1) then begin
        maxs = string(maxi_tmp,format="(e8.2)")
        mins = string(mini_tmp,format="(e9.2)")
      endif else begin
        maxs = string(maxi_tmp,format="(f5.1)")
        mins = string(mini_tmp,format="(f6.1)")
      endelse

      charsize = 1.0
      if ppp gt 9 then charsize = 0.75

      xyouts, pos(0),pos(1), mins, /norm, charsize = charsize
      xyouts, pos(2),pos(1), maxs, /norm, align=1.0, charsize = charsize

      if (nfiles eq 1) then begin
        xyouts, pos(2),pos(3),vars(varlist(i)), 	$
	  align=1.0, /norm
      endif

    endfor

  endfor

endfor

xyouts, 1.0-pos(2), pos(3) + space, 'Sat: F'+tostr(fix(satid)), /norm
xyouts, pos(2), pos(3) + space, strmid(hightime,0,15), /norm, align=1.0

stime = dmsptime(0)
etime = dmsptime(npts-1)

time_axis, stime, etime, s_time_range, e_time_range, 	$
	xtickname, xtitle, xtickvalue, xminor, xtickn

plottime = dmsptime(0:npts-1) - stime
plotdmsp = compare_data(0,0:npts-1)
plotmhd  = compare_data(1,0:npts-1)

vymhd = fltarr(npts)

for i=1,npts-2 do begin
  if (plotmhd(i-1) gt -999.0 and plotmhd(i+1) gt -999.0) then $
    vymhd(i) = (plotmhd(i+1) - plotmhd(i-1)) / $
      compare_data(2,i)/compare_data(3,i)/1000.0 else vymhd(i) = -999.0
endfor

loc = where(plotmhd eq -999.0,count)

if count gt 0 then begin

  istart = loc(0) - 1
  iend   = loc(count-1) + 1

  val_start = plotmhd(istart)
  val_end   = plotmhd(iend)
  slope     = (val_end - val_start) / (iend - istart)

  for i=istart+1,iend-1 do $
     plotmhd(i) = val_start + slope * (i-istart)

endif

loc = where(vymhd eq -999.0,count)

if count gt 0 then begin

  istart = loc(0) - 1
  iend   = loc(count-1) + 1

  val_start = vymhd(istart)
  val_end   = vymhd(iend)
  slope     = (val_end - val_start) / (iend - istart)

  for i=istart+1,iend-1 do $
     vymhd(i) = val_start + slope * (i-istart)

endif

if (dmsp_hem) then vymhd = -1.0*vymhd

yrange = mm([plotdmsp,plotmhd])

get_position, ppp, space, sizes, 2, pos1
get_position, ppp, space, sizes, 3, pos2

pos = pos1
pos(2) = pos2(2)

plot, plottime, plotdmsp, yrange = yrange, $
        xrange = [s_time_range, e_time_range], xtickname = strarr(30)+' ', $
        xtickv=xtickvalue, xticks= xtickn, xminor = xminor,           $
        pos = pos, /noerase, ytitle = 'Potential (kV)'

oplot, plottime, plotmhd, linestyle = 2
oplot, plottime, plottime*0.0, linestyle = 1

xpos = (pos(2) - pos(0))*0.02 + pos(0)
ypos = (pos(3) - pos(1))*0.10 + pos(1)
loc = where(abs(plotdmsp) lt 10.0,count)
if count gt 0 then begin
  pot_rms = sqrt(mean((plotdmsp(loc) - plotmhd(loc))^2.0))
  tot_rms = sqrt(mean((plotdmsp(loc))^2.0))
endif
xyouts, xpos, ypos, 'RMS : '+tostrf(pot_rms), $
        charsize = 0.8, /norm
ypos = (pos(3) - pos(1))*0.05 + pos(1)
xyouts, xpos, ypos, 'RMS : '+tostrf(tot_rms)+' (DMSP Only)', $
        charsize = 0.8, /norm

cpcp_mhd = max(plotmhd) - min(plotmhd)
cpcp_dmsp = max(plotdmsp) - min(plotdmsp)
ypos = (pos(3) - pos(1))*0.93 + pos(1)
xyouts, xpos, ypos, 'CPCP : '+tostrf(cpcp_mhd)+' (M)', $
        charsize = 0.8, /norm

xpos = (pos(2) - pos(0))*0.98 + pos(0)
ypos = (pos(3) - pos(1))*0.93 + pos(1)
xyouts, xpos, ypos, 'CPCP : '+tostrf(cpcp_dmsp)+' (D)', $
        charsize = 0.8, /norm, align = 1.0

get_position, ppp, space, sizes, 4, pos1
get_position, ppp, space, sizes, 5, pos2

pos = pos1
pos(2) = pos2(2)

pos([1,3]) = pos([1,3]) 

efield_dmsp = compare_data(4,0:npts-1)
efield_mhd  = vymhd

loc = where(abs(efield_dmsp) lt 10.0)

plot, plottime(loc), efield_dmsp(loc), /noerase,$
	xrange = [s_time_range, e_time_range], xtickname = xtickname, $
        xtickv=xtickvalue, xticks= xtickn, xminor = xminor,           $
        xtitle = xtitle, pos = pos, ytitle = 'Velocity'

oplot, plottime, efield_mhd, linestyle = 2
oplot, plottime, plottime*0.0, linestyle = 1

xpos = (pos(2) - pos(0))*0.02 + pos(0)
ypos = (pos(3) - pos(1))*0.10 + pos(1)
efield_rms = sqrt(mean((efield_dmsp(loc) - efield_mhd(loc))^2.0))
etot_rms = sqrt(mean((efield_dmsp(loc))^2.0))
xyouts, xpos, ypos, 'RMS : '+tostrf(efield_rms), $
        charsize = 0.8, /norm
ypos = (pos(3) - pos(1))*0.05 + pos(1)
xyouts, xpos, ypos, 'RMS : '+tostrf(etot_rms)+' (DMSP Only)', $
        charsize = 0.8, /norm

if (iter(0) ne 0) then $
  xyouts, 1.01, 0.0, 'iter : '+strcompress(string(iter(0))),/norm, orient = 90

closedevice

openw,1,outfile

c_r_to_a, itime, stime
jd = jday(itime(0),itime(1),itime(2))

printf,1,'  '+tostr(fix(satid(0)))+' '+tostr(itime(0))+' '+$
	chopr('00'+tostr(jd),3)+' '+chopr('0'+tostr(itime(3)),2)+$
        ' '+chopr('0'+tostr(itime(4)),2)+$
        '   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0 99 99  '+$
	chopr('00'+tostr(npts),3)+'  0.00  0.00'

plottime = dmsptime(0:npts-1) - stime

itime(3:5) = 0.0
c_a_to_r, itime, daytime
filetime = dmsptime(0:npts-1) - daytime

for i=0,npts-1 do begin
  printf,1,format = '(f7.0,2F7.3,F7.2,F7.3,I3,F8.2,6F7.2)', $
        filetime(i),vymhd(i),efield_dmsp(i),$
	plotdmsp(i),0.0,1,plotmhd(i), lat_110(index(i)), $
	lat_save(index(i)), mlt_save(index(i)),0.0,0.0,0.0
endfor

nx = n_elements(xpos_output(*,0))
ny = n_elements(xpos_output(0,*))

printf,1,nx,ny
printf,1,mhd_pot_output
printf,1,xpos_output
printf,1,ypos_output

close,1


end    
