
pro figure_out_positions, n, xxyy, pos

   if !p.charsize eq 0.0 then !p.charsize=1.0
   ppp   = 6
   space = max([float(!d.y_ch_size)/float(!d.y_size),$
                float(!d.x_ch_size)/float(!d.x_size)])*3.0*!p.charsize
   pos_space, ppp, space, sizes, ny = 3

   if (n eq 1 or n eq 2) then begin

     get_position, ppp, space, sizes, (n-1)*2, pos1, /rect
     get_position, ppp, space, sizes, (n-1)*2+1, pos2, /rect

     pos2(2) = pos2(2) - (pos2(2) - pos2(0))*0.1

     pos = [pos1(0), pos1(1), pos2(2), pos1(3)]

     width  =  max([xxyy(0),xxyy(2)]) - min([xxyy(0),xxyy(2)])
     height =  max([xxyy(1),xxyy(3)]) - min([xxyy(1),xxyy(3)])
     aspectx = width/height

  endif else begin

     get_position, ppp, space, sizes, n + 1, pos
     aspectx = 1.0

  endelse

  aspectpos = (pos(2)-pos(0))/(pos(3)-pos(1)) $
               *float(!d.x_size)/float(!d.y_size)

  aspectratio = aspectpos/aspectx

  if aspectratio gt 1 then begin
    posmid=(pos(2)+pos(0))/2.
    posdif=(pos(2)-pos(0))/2.
    pos(0)=posmid - posdif/aspectratio
    pos(2)=posmid + posdif/aspectratio
  endif else begin
    posmid=(pos(3)+pos(1))/2.
    posdif=(pos(3)-pos(1))/2.
    pos(1)=posmid - posdif*aspectratio
    pos(3)=posmid + posdif*aspectratio
  endelse

  if (n eq 2) then pos([1,3]) = pos([1,3]) ;+ space/2.0

  if (n eq 3) then pos([0,2]) = pos([0,2]) + space/3.0
  if (n eq 4) then pos([0,2]) = pos([0,2]) - space/3.0

  return

end


if (n_elements(plotgrid) eq 0) then plotgrid = 0
plotgrid = fix(ask('whether you want the grid (1=yes,0=no)',tostr(plotgrid)))

if (n_elements(plotvect) eq 0) then plotvect = 1
plotvect = fix(ask('whether you want the streamlines (1=yes,0=no)', $
	 tostr(plotvect)))

if (n_elements(plotbz0) eq 0) then plotbz0 = 1
plotbz0 = fix(ask('whether you want the bz=0 contour (1=yes,0=no)', $
	 tostr(plotbz0)))

bz0contour = 0.0
realtime = 0

maxcolat = 40.0


filelist_y = findfile('IO2/y=0*.out')
filelist_z = findfile('IO2/z=0*.out')

nfiles_y = n_elements(filelist_y)
nfiles_z = n_elements(filelist_z)

if (nfiles_y eq nfiles_z) then begin

  times = dblarr(nfiles_y)
  iterations = lonarr(nfiles_y)

  for i=0,nfiles_y-1 do begin

    filename = filelist_y(i)

    gettype,filename,filetype,npictinfile

    openfile,10,filename,filetype
    gethead,10,filetype,headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nx,eqpar,variables,pictsize=pictsize
    close,10

    iterations(i) = it
    times(i)      = time

  endfor

  ;-----------------------------------------------------------------------
  ; we now need to figure out the time offset, since the time is the
  ; code time and not the real time.  So, lets search for time = 0,
  ; get the iteration number and look in the ionosphere file to
  ; get the real time.

  loc = where(times eq 0, count)
  if count gt 0 then begin

    it = iterations(loc(count-1)) 
    its = chopr('000000'+tostr(it),6)
    filename = 'ionosphere/in'+its+'.idl'

  endif else begin

    ionofilelist = findfile('ionosphere/it*.idl')
    filename = ionofilelist(0)

  endelse

  openr, 11, filename

  line = ''

  while (strpos(line,'TIME') eq -1) do begin

    readf,11,line

  endwhile

  readf,11, iYear
  readf,11, iMonth
  readf,11, iDay
  readf,11, iHour
  readf,11, iMinute
  readf,11, iSecond
  readf,11, iMillisecond

  close,11

  itime = fix([iYear, iMonth, iDay, iHour, iMinute, iSecond])
  c_a_to_r, itime, basetime

  ; time offset is...

  basetime = basetime + double(iMillisecond)/1000.0 - times(0)

  ;-----------------------------------------------------------------------
  ; now display everything to the user.

  for i = 0, nfiles_y-1 do begin

    c_r_to_a, itime, times(i)+basetime
    print, format='(a,i3,a,f8.1,a,i5,5i3,a,i7)', 'Time # ',i,$
           ' t = ',times(i),'; Real Time = ',itime, '; it = ',iterations(i)

  endfor

  iter = fix(ask('time # of interest',tostr(nfiles_y-1)))

  if (n_elements(xmin) eq 0) then xmin = -50.0
  if (n_elements(xmax) eq 0) then xmax =  20.0
  if (n_elements(ymin) eq 0) then ymin = -20.0
  if (n_elements(ymax) eq 0) then ymax =  20.0

  xmin = float(ask('xmin',tostrf(xmin)))
  xmax = float(ask('xmax',tostrf(xmax)))
  ymin = float(ask('ymin',tostrf(ymin)))
  ymax = float(ask('ymax',tostrf(ymax)))

  xreglimits=[xmin,ymin,xmax,ymax]
  dxreg = xreglimits(2) - xreglimits(0)
  dyreg = xreglimits(3) - xreglimits(1)

  nxreg=[-dxreg,-dyreg]
  transform='regular'
  dotransform='n'
  symmtri = 1

  filename = filelist_y(iter)
  read_idl_file, filename, npict, nxreg, xreglimits, transform, $
                   nfile, physics, ndim, gencoord, it, time,      $
                   nx, x, xreg, nw, w, wreg, wnames, variables,symmtri
  yplane_wreg = wreg
  yplane_xreg = xreg
  xy = x

  ; Let's Try to figure out how big the body is:

  r = reform(sqrt(xreg(*,*,0)^2 + xreg(*,*,1)^2))
  current = reform(wreg(*,*,8)^2 + wreg(*,*,9)^2 + wreg(*,*,10)^2)
  loc_current = where(current eq 0.0, count)
  rbody = 0.0
  for i=0L,count-1 do begin
    if (r(loc_current(i)) gt rbody and r(loc_current(i)) lt 7.5) then begin
      rbody = r(loc_current(i))
    endif
  endfor

  filename = filelist_z(iter)
  read_idl_file, filename, npict, nxreg, xreglimits, transform, $
                   nfile, physics, ndim, gencoord, it, time,      $
                   nx, x, xreg, nw, w, wreg, wnames, variables,symmtri
  zplane_wreg = wreg
  zplane_xreg = xreg
  xz = x

  its = chopr('000000'+tostr(it),6)

  for i=0,n_elements(wnames)-1 do print, i,". ",wnames(i)
  if (n_elements(iVarMHD) eq 0) then iVarMHD = 1
  iVarMHD = fix(ask('variable to plot',tostr(iVarMHD)))

  ilog = 'n'
  ilog = mklower(ask('whether you want this variable to be alog10ed',ilog))

  yplane_data = yplane_wreg(*,*,iVarMHD)
  zplane_data = zplane_wreg(*,*,iVarMHD)
  VarName = wnames(iVarMHD)

  if (strpos(ilog,'y') eq 0) then begin
    yplane_data = alog10(yplane_data)
    zplane_data = alog10(zplane_data)
    VarName = "alog10("+VarName+")"
  endif

  c_r_to_a, itime, times(iter)+basetime

  iYear   = itime(0) mod 100
  sYear   = chopr('00'+tostr(iYear),2)
  sMonth  = chopr('00'+tostr(iTime(1)),2)
  sDay    = chopr('00'+tostr(iTime(2)),2)
  sHour   = chopr('00'+tostr(iTime(3)),2)
  sMinute = chopr('00'+tostr(iTime(4)),2)
  sSecond = chopr('00'+tostr(iTime(5)),2)
  sMilli  = chopr('0000'+tostr((times(iter) - long(times(iter)))*1000),4)

  psfilename = 'plot'+sYear+sMonth+sDay+'_'+sHour+sMinute+sSecond+'.ps'

  setdevice, psfilename, 'p', 4, 0.95

  plotdumb

  makect, 'mid'

  ;----------------------------------------------------------------------
  ; make plot in the X-Z plane

  figure_out_positions, 1, xreglimits, pos

  mini = min(yplane_data)
  maxi = max(yplane_data)

  if (maxi eq mini) then begin
    maxi = mini + 0.5
    mini = maxi - 1.0
  endif

  if (mini lt 0.0) then begin

    ; now this is tricky, since IDL screws up the regular grid.
    ni = n_elements(yplane_data(*,0))
    nj = n_elements(yplane_data(0,*))

    mini = min(yplane_data(1:ni-2,1:nj-2))

    ; this means that it truly is less than 0.0, instead of an artifact

    if (mini lt 0.0) then begin
      if (maxi/1.0 lt abs(mini)) then begin
        maxi = max([abs(mini),maxi])*1.01 
        mini = -maxi
      endif else max = abs(mini)
    endif else begin
      yplane_data(0,*) = yplane_data(1,*)
      yplane_data(ni-1,*) = yplane_data(ni-2,*)
      yplane_data(*,0) = yplane_data(*,1)
      yplane_data(*,nj-1) = yplane_data(*,nj-2)
      mini = 0.0
    endelse

  endif

;  mini = min(yplane_data)
;  maxi = max(yplane_data)

  c_levels = (maxi-mini)*findgen(30)/29.0 + mini
  c_colors = 253.0*findgen(29)/28.0 + 1.0
  c_colors(14) = 128.0

  contour, yplane_data, yplane_xreg(*,*,0), yplane_xreg(*,*,1), $
           /follow, nlevels = 30, pos = pos, /noerase, $
           xrange = [xmax, xmin], yrange = [ymin, ymax], $
	   xstyle = 1, ystyle = 1, /cell_fill, ytitle = 'Z (Re)', $
	   xtickname = strarr(30)+' ', levels = c_levels, c_colors = c_colors

  theta = findgen(37)*!pi*2.0/36.0
  r = rbody
  polyfill, r*cos(theta), r*sin(theta), color = 0

  ;-------------
  ; Add Traces

  if (plotvect) then begin

    if (strpos(wnames(4),'bx') lt 0 or strpos(wnames(6),'bz') lt 0) then begin
      print, "Variables out of order. I am to stupid to figure out what"
      print, "to do. Contact Aaron and yell at him!!"
      stop
    endif

    f1 = yplane_wreg(*,*,4)
    f2 = yplane_wreg(*,*,6)

    velvector = 99

    velpos = fltarr(velvector,2)

    velpos(0:velvector/3-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
    velpos(0:velvector/3-1,1) = (ymax-ymin) * findgen(velvector/3)/(velvector/3) / 2.0

    velpos(velvector/3:2*velvector/3-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
    velpos(velvector/3:2*velvector/3-1,1) = 0.0

    velpos(2*velvector/3:velvector-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
    velpos(2*velvector/3:velvector-1,1) = - (ymax-ymin) * findgen(velvector/3)/(velvector/3) / 2.0

    ; normalization
    eps=1.e-30
    v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
    ; arrows
    vector,v1,v2,yplane_xreg(*,*,0),yplane_xreg(*,*,1), $
           NVECS=velvector,MAXVAL=1.,                   $
           NSTEP=6,LENGTH=0.012,HEAD=0.5,               $
           DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = [xmax, xmin], yrange = [ymin, ymax],$
	   pos = pos, xtickname = strarr(30)+' '
    ; streamline along v1;v2
    vector,v1,v2,yplane_xreg(*,*,0),yplane_xreg(*,*,1), $
           NVECS=velvector,MAXVAL=1.,$
           NSTEP=1000,LENGTH=2.,HEAD=0.,$
           DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = [xmax, xmin], yrange = [ymin, ymax],$
	   pos = pos, xtickname = strarr(30)+' '
    ; streamline in the other direction
    v1=-v1 & v2=-v2
    vector,v1,v2,yplane_xreg(*,*,0),yplane_xreg(*,*,1), $
           NVECS=velvector,MAXVAL=1.,$
           NSTEP=1000,LENGTH=2.,HEAD=0.,$
           DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = [xmax, xmin], yrange = [ymin, ymax],$
	   pos = pos, xtickname = strarr(30)+' '

  endif

  ;-------------
  ; Add Grid

  if plotgrid then begin
    loc = where(xy(*,0,0) ge xmin and xy(*,0,0) le xmax and $
                xy(*,0,1) ge ymin and xy(*,0,1) le ymax, count)
    if (count gt 0) then $
      oplot, xy(loc,0,0), xy(loc,0,1), psym = 1, symsize=0.2
  endif

  ;----------------
  ; Add Color Table

  ncolors = 255
  minmax = mm(c_levels)
  title = ""

  if (pos(2) ge 0.92) then begin
    ct_pos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
    plotct, ncolors, ct_pos, minmax, title, /right
  endif else begin
    ct_pos = [pos(2)+0.07,pos(1),pos(2)+0.09,pos(3)]
    plotct, ncolors, ct_pos, minmax, title
  endelse

  ; Put the max and min on the plot

  mini_tmp = min(yplane_data)
  maxi_tmp = max(yplane_data)

  r = (maxi_tmp - mini_tmp)/50.0

  plots, [1.3,2.0], [mini_tmp, mini_tmp]
  plots, [1.3,1.6], [mini_tmp, mini_tmp+r]
  plots, [1.3,1.6], [mini_tmp, mini_tmp-r]

  plots, [1.3,2.0], [maxi_tmp, maxi_tmp]
  plots, [1.3,1.6], [maxi_tmp, maxi_tmp+r]
  plots, [1.3,1.6], [maxi_tmp, maxi_tmp-r]

  if (abs(mini_tmp) ge 1000.0) or (abs(maxi_tmp) ge 1000.0) or        $
     (abs(maxi_tmp) lt 0.1) then begin
    maxs = string(maxi_tmp,format="(e8.2)")
    mins = string(mini_tmp,format="(e9.2)")
  endif else begin
    maxs = string(maxi_tmp,format="(f6.2)")
    mins = string(mini_tmp,format="(f7.2)")
  endelse

  xyouts, 2.1, mini_tmp, mins, charsize = 0.8
  xyouts, 2.1, maxi_tmp, maxs, charsize = 0.8

  ;----------------------------------------------------------------------
  ; make plot in the X-Y plane

  figure_out_positions, 2, xreglimits, pos

  mini = min(zplane_data)
  maxi = max(zplane_data)

  if (maxi eq mini) then maxi = mini + 1.0

  if (mini lt 0.0) then begin

    ; now this is tricky, since IDL screws up the regular grid.
    ni = n_elements(zplane_data(*,0))
    nj = n_elements(zplane_data(0,*))

    mini = min(zplane_data(1:ni-2,1:nj-2))

    ; this means that it truly is less than 0.0, instead of an artifact

    if (mini lt 0.0) then begin
      if (maxi/1.0 lt abs(mini)) then begin
        maxi = max([abs(mini),maxi])*1.01 
        mini = -maxi
      endif
    endif else begin
      zplane_data(0,*) = zplane_data(1,*)
      zplane_data(ni-1,*) = zplane_data(ni-2,*)
      zplane_data(*,0) = zplane_data(*,1)
      zplane_data(*,nj-1) = zplane_data(*,nj-2)
      mini = 0.0
    endelse

  endif

;  mini = min(zplane_data)
;  maxi = max(zplane_data)

  c_levels = (maxi-mini)*findgen(30)/29.0 + mini
  c_colors = 253.0*findgen(29)/28.0 + 1.0
  c_colors(14) = 128.0

  contour, zplane_data, zplane_xreg(*,*,0), zplane_xreg(*,*,1), $
           /follow, nlevels = 30, pos = pos, /noerase, $
           xrange = [xmax, xmin], yrange = [ymax, ymin], $
	   xstyle = 1, ystyle = 1, /cell_fill, levels = c_levels, $
	   xtitle = 'X (Re)', ytitle = 'Y (Re)', c_colors = c_colors
  theta = findgen(37)*!pi*2.0/36.0
  r = rbody
  polyfill, r*cos(theta), r*sin(theta), color = 0

  ;-------------
  ; Add Traces

  if (plotvect) then begin

    if (strpos(wnames(1),'ux') lt 0 or strpos(wnames(2),'uy') lt 0) then begin
      print, "Variables out of order. I am to stupid to figure out what"
      print, "to do. Contact Aaron and yell at him!!"
      stop
    endif

    f1 = zplane_wreg(*,*,1)
    f2 = zplane_wreg(*,*,2)

    velpos = fltarr(velvector,2)

    ; First rake near the front boundary
    velpos(0:velvector/3-1,0) = xmax - (xmax-xmin)/20.0
    velpos(0:velvector/3-1,1) = (ymax - ymin) * findgen(velvector/3)/(velvector/3) + ymin

    ; Second rake just behind the Earth
    velpos(velvector/3:2*velvector/3-1,0) = 0.0 - (xmax-xmin)/20.0
    velpos(velvector/3:2*velvector/3-1,1) = (ymax - ymin) * findgen(velvector/3)/(velvector/3) + ymin

    ; Third rake near the back
    velpos(2*velvector/3:velvector-1,0) = -20.0
    velpos(2*velvector/3:velvector-1,1) = (ymax - ymin) * findgen(velvector/3)/(velvector/3) + ymin


    ; normalization
    eps=1.e-30
    v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
    ; arrows
    vector,v1,v2,zplane_xreg(*,*,0),zplane_xreg(*,*,1), $
           NVECS=velvector,MAXVAL=1.,                   $
           NSTEP=6,LENGTH=0.012,HEAD=0.5,               $
           DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = [xmax, xmin], yrange = [ymax, ymin],$
	   pos = pos, xtickname = strarr(30)+' '
    ; streamline along v1;v2
    vector,v1,v2,zplane_xreg(*,*,0),zplane_xreg(*,*,1), $
           NVECS=velvector,MAXVAL=1.,$
           NSTEP=1000,LENGTH=2.,HEAD=0.,$
           DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = [xmax, xmin], yrange = [ymax, ymin],$
	   pos = pos, xtickname = strarr(30)+' '
    ; streamline in the other direction
    v1=-v1 & v2=-v2
    vector,v1,v2,zplane_xreg(*,*,0),zplane_xreg(*,*,1), $
           NVECS=velvector,MAXVAL=1.,$
           NSTEP=1000,LENGTH=2.,HEAD=0.,$
           DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = [xmax, xmin], yrange = [ymax, ymin],$
	   pos = pos, xtickname = strarr(30)+' '

  endif

  ;-------------
  ; Add Traces

  if plotbz0 then begin

    if (strpos(wnames(6),'bz') lt 0) then begin
      print, "Variables out of order. I am to stupid to figure out what"
      print, "to do. Contact Aaron and yell at him!!"
      stop
    endif

    bz = zplane_wreg(*,*,6)

    contour, bz, zplane_xreg(*,*,0), zplane_xreg(*,*,1), $
             levels = [bz0contour], c_thick = 6, xstyle = 5, ystyle = 5, $
	     pos = pos, /noerase, $
             xrange = [xmax, xmin], yrange = [ymax, ymin], color = 254

    contour, bz, zplane_xreg(*,*,0), zplane_xreg(*,*,1), $
             levels = [bz0contour-2.0], c_thick = 2, xstyle = 5, ystyle = 5, $
	     pos = pos, /noerase, $
             xrange = [xmax, xmin], yrange = [ymax, ymin], color = 1

    contour, bz, zplane_xreg(*,*,0), zplane_xreg(*,*,1), $
             levels = [bz0contour+2.0], c_thick = 2, xstyle = 5, ystyle = 5, $
	     pos = pos, /noerase, $
             xrange = [xmax, xmin], yrange = [ymax, ymin], color = 1

  endif


  if plotgrid then begin
    loc = where(xz(*,0,0) ge xmin and xz(*,0,0) le xmax and $
                xz(*,0,1) ge ymin and xz(*,0,1) le ymax, count)
    if (count gt 0) then $
      oplot, xz(loc,0,0), xz(loc,0,1), psym = 1, symsize = 0.2
  endif

  ncolors = 255
  minmax = mm(c_levels)
  title = ""
  if (pos(2) ge 0.92) then begin
    ct_pos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
    plotct, ncolors, ct_pos, minmax, title, /right
  endif else begin
    ct_pos = [pos(2)+0.07,pos(1),pos(2)+0.09,pos(3)]
    plotct, ncolors, ct_pos, minmax, title
  endelse

  ; Put the max and min on the plot

  mini_tmp = min(zplane_data)
  maxi_tmp = max(zplane_data)

  r = (maxi_tmp - mini_tmp)/50.0

  plots, [1.3,2.0], [mini_tmp, mini_tmp]
  plots, [1.3,1.6], [mini_tmp, mini_tmp+r]
  plots, [1.3,1.6], [mini_tmp, mini_tmp-r]

  plots, [1.3,2.0], [maxi_tmp, maxi_tmp]
  plots, [1.3,1.6], [maxi_tmp, maxi_tmp+r]
  plots, [1.3,1.6], [maxi_tmp, maxi_tmp-r]

  plots, [1.3,2.0], [mini_tmp, mini_tmp]
  plots, [1.3,2.0], [maxi_tmp, maxi_tmp]

  if (abs(mini_tmp) ge 1000.0) or (abs(maxi_tmp) ge 1000.0) or        $
     (abs(maxi_tmp) lt 0.1) then begin
    maxs = string(maxi_tmp,format="(e8.2)")
    mins = string(mini_tmp,format="(e9.2)")
  endif else begin
    maxs = string(maxi_tmp,format="(f6.2)")
    mins = string(mini_tmp,format="(f7.2)")
  endelse

  xyouts, 2.1, mini_tmp, mins, charsize = 0.8
  xyouts, 2.1, maxi_tmp, maxs, charsize = 0.8

  ;---------------------------------------------------------------------------
  ; Next we do the ionosphere plots


  ; Figure out the file name

print, times(iter)

  c_r_to_a, itime, times(iter)+basetime

print, itime

  iYear   = itime(0) mod 100
  sYear   = chopr('00'+tostr(iYear),2)
  sMonth  = chopr('00'+tostr(iTime(1)),2)
  sDay    = chopr('00'+tostr(iTime(2)),2)
  sHour   = chopr('00'+tostr(iTime(3)),2)
  sMinute = chopr('00'+tostr(iTime(4)),2)
  sSecond = chopr('00'+tostr(iTime(5)),2)
  sMilli  = chopr('0000'+tostr((times(iter) - long(times(iter)))*1000),4)

  ionofilename = 'ionosphere/it'+sYear+sMonth+sDay+'_'+$
                 sHour+sMinute+sSecond+'_*.idl'
  ionofilename = findfile(ionofilename)

  if (strlen(ionofilename(0)) eq 0) then begin
    its = chopr('000000'+tostr(it),6)
    ionofilename = 'ionosphere/in'+its+'.idl'
    ionofilename = findfile(ionofilename)
    if (strlen(ionofilename(0)) eq 0) then begin

      if (iTime(5) ge 30) then begin
        sMinute = chopr('00'+tostr(iTime(4)+1),2)
      endif else begin
        sMinute = chopr('00'+tostr(iTime(4)-1),2)
      endelse
      ionofilename = 'ionosphere/it'+sYear+sMonth+sDay+'_'+$
                   sHour+sMinute+'*.idl'
      ionofilename = findfile(ionofilename)

      if (strlen(ionofilename(0)) eq 0) then begin
        print, "Could not find ionosphere file!!!"
        print, "should be : ",'ionosphere/it'+sYear+sMonth+sDay+'_'+$
                 sHour+sMinute+sSecond+'_*.idl'
        stop
      endif
    endif
  endif

  ionofilename = ionofilename(0)

  read_iono_file, ionofilename, ionodata, time, theta, phi, variables = [4,5]

  loc = where(reform(theta(0,0,*)) le maxcolat)
  rang = reform(theta(0,*,loc))
  lons = reform(phi(0,*,loc))*!pi/180.0 + !pi/2
  xpos = rang*cos(lons)
  ypos = rang*sin(lons)

  ocflbfilename = 'IO2/ocflb_n'+its+'.dat'
  list = findfile(ocflbfilename)
  if (strlen(list(0)) gt 0) then begin
    plot_ocflb = 1
    nmlts = 0
    openr,10,list(0)
    readf,10,nmlts
    mlts = fltarr(nmlts)
    ocflb = fltarr(nmlts)
    readf,10,mlts
    readf,10,ocflb
    close,10

    for i=0,(nmlts-1)/2 do begin
      ocflb(i) = (ocflb(i) + ocflb(nmlts-1-i))/2.0
    endfor
    for i=(nmlts-1)/2+1,nmlts-1 do begin
      ocflb(i) = ocflb(nmlts-1-i)
    endfor

    print, 'Finished Reading File ',ocflbfilename
    ocx    = (90.0-ocflb)*cos(mlts*!pi/12.0 - !pi/2.0)
    ocy    = (90.0-ocflb)*sin(mlts*!pi/12.0 - !pi/2.0)
  endif else plot_ocflb = 0

  for iVar = 0, 1 do begin

    ; (the 0 here is because we are only dealing with one file)

    mini = min(ionodata(0,0,iVar,*,loc))
    maxi = max(ionodata(0,0,iVar,*,loc))

    if (maxi eq mini) then maxi = mini + 1.0
    if (mini lt 0.0) then begin
      maxi = max([abs(mini),maxi]) 
      mini = -maxi
      levels = fltarr(8)
      levels(0:3) = -maxi*(4-findgen(4))/4.0
      levels(4:7) = maxi*(findgen(4)+1)/4.0
    endif else begin
      levels = maxi*findgen(9)/8.0 + mini
    endelse

    c_levels = (maxi-mini)*findgen(30)/29.0 + mini
    c_colors = 253.0*findgen(29)/28.0 + 1.0
    c_colors(14) = 128.0

    figure_out_positions, 3 + iVar, xreglimits, pos

    contour, ionodata(0,0,iVar,*,loc), xpos, ypos,        $
        /follow, nlevels=30, /noerase, $
	pos = pos, xstyle = 5, ystyle = 5, $
        xrange = [-maxcolat, maxcolat], yrange = [-maxcolat, maxcolat], $
	levels = c_levels, /cell_fill, c_colors = c_colors

    contour, ionodata(0,0,iVar,*,loc), xpos, ypos,        $
        /follow, /noerase, $
	pos = pos, xstyle = 5, ystyle = 5, $
        xrange = [-maxcolat, maxcolat], yrange = [-maxcolat, maxcolat], $
	levels = levels

    if plotgrid then begin
      oplot, xpos, ypos, psym = 1, symsize = 0.2
    endif

    no00 = 0
    no06 = 1
    no12 = 0
    no18 = 1

    plotmlt, maxcolat, no00 = no00, no06 = no06, no12 = no12, no18 = no18

    if (plot_ocflb) then oplot, ocx, ocy, thick = 5, color = 0

    mini_tmp = min(ionodata(0,0,iVar,*,loc))
    maxi_tmp = max(ionodata(0,0,iVar,*,loc))

    if (abs(mini_tmp) gt 1000.0) or (abs(maxi_tmp) gt 1000.0) or        $
       (abs(maxi_tmp) lt 0.1) then begin
      maxs = string(maxi_tmp,format="(e8.2)")
      mins = string(mini_tmp,format="(e9.2)")
    endif else begin
      maxs = string(maxi_tmp,format="(f6.2)")
      mins = string(mini_tmp,format="(f7.2)")
    endelse

    xyouts, pos(0),pos(1), mins, /norm
    xyouts, pos(2),pos(1), maxs, /norm, align=1.0

    if (iVar eq 0) then begin
      ct_pos = [pos(0)-0.03,pos(1),pos(0)-0.01,pos(3)]
      title = 'FAC (!Mm!XA/m!U2!N)'
      right = 0
    endif else begin
      ct_pos = [pos(2)+0.01,pos(1),pos(2)+0.03,pos(3)]
      title = 'Potential (kV)'
      right = 1
    endelse

    ncolors = 255
    minmax = mm(c_levels)
    plotct, ncolors, ct_pos, minmax, title, right=right

  endfor

  ;-------------------------------------------------------------------------
  ; Next, lets label some times and iterations

  figure_out_positions, 1, xreglimits, pos

  if realtime then begin
    c_a_to_s, itime, stime
    figure_out_positions, 1, xreglimits, pos
    xyouts, pos(0)-0.05, pos(3) + 0.03, 'Time : '+strmid(stime,0,15)+' UT', $
            /norm, charsize = 1.2
  endif else begin
    c_r_to_a, itime, times(iter)
    c_a_to_s, itime, stime
    xyouts, pos(0)-0.05, pos(3) + 0.03, $
            'Simulation time : '+strmid(stime,9,15), /norm, charsize = 1.2
  endelse

  xyouts, pos(2), pos(3) + 0.03, 'It : '+tostr(it), alignment = 1.0, $
          /norm, charsize = 1.2

  xyouts, pos(0)-0.08, pos(3) - 0.02, '(A)', /norm, charsize = 1.2

  figure_out_positions, 2, xreglimits, pos
  xyouts, pos(0)-0.08, pos(3) - 0.02, '(B)', /norm, charsize = 1.2

  xyouts, (pos(0) + pos(2))/2.0, pos(3) + 0.01, $
          "Plots of "+VarName+" in the X-Y and X-Z planes", $
          /norm, charsize=1.2, alignment = 0.5

  figure_out_positions, 3, xreglimits, pos
  xyouts, pos(0)-0.13, pos(3) - 0.03, '(C)', /norm, charsize = 1.2

  figure_out_positions, 4, xreglimits, pos
  xyouts, pos(0), pos(3) - 0.03, '(D)', /norm, charsize = 1.2

  spawn,'pwd',pwd
  xyouts, 0.0, -0.02, pwd, /norm, charsize = 0.5

  closedevice

endif else begin

  print, 'nfiles_y ne nfiles_z'

endelse


end
