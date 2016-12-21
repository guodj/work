
plotgrid = 0
plotvect = 1
plotbz0  = 0
bz0contour = 0.0
realtime = 0
bw = 0
ct = 1
rbody = 1.6

if (n_elements(plane) eq 0) then plane = 'y'
plane = mklower(ask('which plane (x/z/y)',plane))
if (strpos(plane,'y') eq 0) then Plane = 1 
if (strpos(plane,'z') eq 0) then Plane = 2 
if (strpos(plane,'x') eq 0) then Plane = 0

if (Plane eq 0) then filelist = findfile('IO2/x=0*.out') 
if (Plane eq 1) then filelist = findfile('IO2/y=0*.out') 
if (Plane eq 2) then filelist = findfile('IO2/z=0*.out') 

nfiles = n_elements(filelist)

times = dblarr(nfiles)
iterations = lonarr(nfiles)

for i=0,nfiles-1 do begin

  filename = filelist(i)

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

  print, 'Help!!  There is no time = 0. Ask Ridley to reprogram this!!"
  ionolist = findfile('ionosphere/i*.idl')
  filename = ionolist(0)

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

basetime = basetime + double(iMillisecond)/1000.0

;-----------------------------------------------------------------------
; now display everything to the user.

for i = 0, nfiles-1 do begin

  c_r_to_a, itime, times(i)+basetime
  print, format='(a,i3,a,f8.1,a,i5,5i3,a,i7)', 'Time # ',i,$
	 ' t = ',times(i),'; Real Time = ',itime, '; it = ',iterations(i)

endfor

iter = fix(ask('time # of interest',tostr(nfiles-1)))

if (n_elements(plotgrid) eq 0) then plotgrid = 0
plotgrid = fix(ask('whether you want the grid (1=yes,0=no)',tostr(plotgrid)))

if (n_elements(plotvect) eq 0) then plotvect = 1
plotvect = fix(ask('whether you want the streamlines (1=yes,0=no)', $
	 tostr(plotvect)))

if (n_elements(plotbz0) eq 0) then plotbz0 = 1
plotbz0 = fix(ask('whether you want the bz=0 contour (1=yes,0=no)', $
	 tostr(plotbz0)))

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

filename = filelist(iter)
read_idl_file, filename, npict, nxreg, xreglimits, transform, $
		 nfile, physics, ndim, gencoord, it, time,      $
		 nx, x, xreg, nw, w, wreg, wnames, variables,symmtri


r = reform(sqrt(xreg(*,*,0)^2 + xreg(*,*,1)^2))
current = reform(wreg(*,*,8)^2 + wreg(*,*,9)^2 + wreg(*,*,10)^2)
loc_current = where(current eq 0.0, count)
;rbody = 0.0
;for i=0,count-1 do begin
;  if (r(loc_current(i)) gt rbody and r(loc_current(i)) lt 7.5) then begin
;    rbody = r(loc_current(i))
;  endif
;endfor

; make some new variables.....

bmag = sqrt(wreg(*,*,4)^2.0 + wreg(*,*,5)^2.0 + wreg(*,*,6)^2.0) * 1e-9

bx = wreg(*,*,4) * 1e-9
by = wreg(*,*,5) * 1e-9
bz = wreg(*,*,6) * 1e-9

bxnorm = bx / bmag
bynorm = by / bmag
bznorm = bz / bmag

ux = wreg(*,*,1) * 1000.0
uy = wreg(*,*,2) * 1000.0
uz = wreg(*,*,3) * 1000.0

uxalong = wreg(*,*,1) * abs(bxnorm)
uyalong = wreg(*,*,2) * abs(bynorm)
uzalong = wreg(*,*,3) * abs(bznorm)

uxperp = wreg(*,*,1) - uxalong
uyperp = wreg(*,*,2) - uyalong
uzperp = wreg(*,*,3) - uzalong

ex = - (uy*bz - by*uz)
ey =   (ux*bz - bx*uz)
ez = - (ux*by - bx*uy)

;loc = where(ey lt 1.0e-3,count)
;if count gt 0 then ey(loc) = 1.0e-3

;wreg(*,*,1) = ex
;wreg(*,*,2) = abs(ey)
;wreg(*,*,3) = ez
;
;uperp = sqrt(uxperp^2+uyperp^2+uzperp^2)
;loc = where(uperp lt 1.0,count)
;if count gt 0 then uperp(loc) = 1.0
;
;uperp = sqrt((ey*bz)^2 + (bx*ey)^2)/bmag^2 / 1000.0
;loc = where(uperp lt 1.0,count)
;if count gt 0 then uperp(loc) = 1.0
;
;
;wreg(*,*,1) = uperp
;;wreg(*,*,2) = uyperp
;;wreg(*,*,3) = uzperp
;
;;wreg(*,*,1) = ux/1000.0
;wreg(*,*,2) = uy/1000.0
;wreg(*,*,3) = uz/1000.0


its = chopr('000000'+tostr(it),6)

for i=0,n_elements(wnames)-1 do print, i,". ",wnames(i)
if (n_elements(iVarMHD) eq 0) then iVarMHD = 1
iVarMHD = fix(ask('variable to plot',tostr(iVarMHD)))

ilog = 'n'
ilog = mklower(ask('whether you want this variable to be alog10ed',ilog))

data = wreg(*,*,iVarMHD)
VarName = wnames(iVarMHD)

if (strpos(ilog,'y') eq 0) then begin
  data = alog10(data)
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

if (Plane eq 0) then begin
  yzs = 'Z' 
  xys = 'Y'
endif 
if (Plane eq 1) then begin
  yzs = 'Z' 
  xys = 'X'
endif
if (Plane eq 2) then begin
  yzs = 'Y' 
  xys = 'X' 
endif

psfilename = 'plot'+sYear+sMonth+sDay+'_'+sHour+sMinute+sSecond+'_'+$
             mklower(xys)+mklower(yzs)+'.ps'

setdevice, psfilename, 'p', 4, 0.95

plotdumb

if (bw eq 0) then makect, 'mid' else makect, 'wb'

;----------------------------------------------------------------------
; make plot in the X-Z plane

figure_out_positions, 1, xreglimits, pos, ppp = 1

mini = min(data)
maxi = max(data)

print, mini, maxi

mini = -500.0
maxi = 500.0

if (maxi eq mini) then maxi = mini + 1.0
if (mini lt 0.0) then begin
  maxi = max([abs(mini),maxi])*1.01 
  mini = -maxi
endif

c_levels = (maxi-mini)*findgen(30)/29.0 + mini
c_colors = 253.0*findgen(29)/28.0 + 1.0
c_colors(14) = 128.0

if (Plane eq 0) then yrange = [ymin,ymax] 
if (Plane eq 1) then yrange = [ymin,ymax] 
if (plane eq 2) then yrange = [ymax,ymin]

xrange = [xmax, xmin]
if (plane eq 0) then xrange = [xmin, xmax]

contour, data, xreg(*,*,0), xreg(*,*,1), $
	 nlevels = 30, pos = pos, /noerase, $
	 xrange = xrange, yrange = yrange, $
	 xstyle = 1, ystyle = 1, /cell_fill, ytitle = yzs+' (Re)', $
	 levels = c_levels, c_colors = c_colors, xtitle = xys+'(Re)'

theta = findgen(181)*!pi*2.0/181.0
r = rbody
loc = where(cos(theta) ge 0 and sin(theta) ge 0.0)
;polyfill, [0.0,r*cos(theta(loc)),0.0], [0.0,r*sin(theta(loc)),0.0], color = 0
polyfill, r*cos(theta), r*sin(theta), color = 0

;-------------
; Add Traces

if (plotvect) then begin

  if (Plane eq 0) then begin

    if (strpos(wnames(5),'by') lt 0 or strpos(wnames(6),'bz') lt 0) then begin
      print, "Variables out of order. I am to stupid to figure out what"
      print, "to do. Contact Aaron and yell at him!!"
      stop
    endif

    f1 = wreg(*,*,5)
    f2 = wreg(*,*,6)

    velvector = 99
  
    velpos = fltarr(velvector,2)

;    velpos(0:velvector/2-1,0) = (xmax - xmin) * findgen(velvector/2)/(velvector/2) + xmin
;    velpos(0:velvector/2-1,1) = 0.1
;    velpos(velvector/2:velvector-1,0) = (xmax - xmin) * findgen(velvector/2)/(velvector/2) + xmin
;    velpos(velvector/2:velvector-1,1) = 4.0
  
    velpos(0:velvector/3-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
    velpos(0:velvector/3-1,1) = (ymax-ymin) * findgen(velvector/3)/(velvector/3) / 2.0
  
    velpos(velvector/3:2*velvector/3-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
    velpos(velvector/3:2*velvector/3-1,1) = 0.0
  
    velpos(2*velvector/3:velvector-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
    velpos(2*velvector/3:velvector-1,1) = - (ymax-ymin) * findgen(velvector/3)/(velvector/3) / 2.0

  endif

  if (Plane eq 1) then begin

    if (strpos(wnames(4),'bx') lt 0 or strpos(wnames(6),'bz') lt 0) then begin
      print, "Variables out of order. I am to stupid to figure out what"
      print, "to do. Contact Aaron and yell at him!!"
      stop
    endif

    f1 = wreg(*,*,4)
    f2 = wreg(*,*,6)

    velvector = 61
  
    velpos = fltarr(velvector,2)

;    velpos(0:velvector/2-1,0) = (xmax - xmin) * findgen(velvector/2)/(velvector/2) + xmin
;    velpos(0:velvector/2-1,1) = 0.1
;    velpos(velvector/2:velvector-1,0) = (xmax - xmin) * findgen(velvector/2)/(velvector/2) + xmin
;    velpos(velvector/2:velvector-1,1) = 4.0
  
;    velpos(0:velvector/3-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
;    velpos(0:velvector/3-1,1) = (ymax-ymin) * findgen(velvector/3)/(velvector/3) / 2.0
;  
;    velpos(velvector/3:2*velvector/3-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
;    velpos(velvector/3:2*velvector/3-1,1) = 0.0
;  
;    velpos(2*velvector/3:velvector-1,0) = (xmax - xmin) * findgen(velvector/3)/(velvector/3) + xmin
;    velpos(2*velvector/3:velvector-1,1) = - (ymax-ymin) * findgen(velvector/3)/(velvector/3) / 2.0

    r = 1.5
    t = findgen(velvector)/(velvector-1) * !pi * 2
    velpos(*,0) = r * cos(t)
    velpos(*,1) = r * sin(t)

  endif

  if (Plane eq 2) then begin

    if (strpos(wnames(1),'ux') lt 0 or strpos(wnames(2),'uy') lt 0) then begin
      print, "Variables out of order. I am to stupid to figure out what"
      print, "to do. Contact Aaron and yell at him!!"
      stop
    endif

    f1 = wreg(*,*,1)
    f2 = wreg(*,*,2)

    velvector = 99

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

  endif

  ; normalization
  eps=1.e-30
  v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
  ; arrows
  vector,v1,v2,xreg(*,*,0),xreg(*,*,1), $
	 NVECS=velvector,MAXVAL=1.,                   $
	 NSTEP=6,LENGTH=0.012,HEAD=0.5,               $
	 DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = xrange, yrange = yrange,$
	   pos = pos, xtickname = strarr(30)+' '
  ; streamline along v1;v2
  vector,v1,v2,xreg(*,*,0),xreg(*,*,1), $
	 NVECS=velvector,MAXVAL=1.,$
	 NSTEP=1000,LENGTH=2.,HEAD=0.,$
	 DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = xrange, yrange = yrange, $
	   pos = pos, xtickname = strarr(30)+' '
  ; streamline in the other direction
  v1=-v1 & v2=-v2
  vector,v1,v2,xreg(*,*,0),xreg(*,*,1), $
	 NVECS=velvector,MAXVAL=1.,$
	 NSTEP=1000,LENGTH=2.,HEAD=0.,$
	 DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE,   $
	   xrange = xrange, yrange = yrange, $[ymin, ymax],$
	   pos = pos, xtickname = strarr(30)+' '

endif

;-------------
; Add Grid

if plotgrid then begin
  loc = where(x(*,0,0) ge xmin and x(*,0,0) le xmax and $
	      x(*,0,1) ge ymin and x(*,0,1) le ymax, count)
  if (count gt 0) then $
    oplot, x(loc,0,0), x(loc,0,1), psym = 1, symsize=0.2
endif

;----------------
; Add Color Table

if (ct) then begin

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

;  ; Put the max and min on the plot
;
;  mini_tmp = min(data)
;  maxi_tmp = max(data)
;
;  r = (maxi_tmp - mini_tmp)/50.0
;
;  plots, [1.3,2.0], [mini_tmp, mini_tmp]
;  plots, [1.3,1.6], [mini_tmp, mini_tmp+r]
;  plots, [1.3,1.6], [mini_tmp, mini_tmp-r]
;
;  plots, [1.3,2.0], [maxi_tmp, maxi_tmp]
;  plots, [1.3,1.6], [maxi_tmp, maxi_tmp+r]
;  plots, [1.3,1.6], [maxi_tmp, maxi_tmp-r]
;
;  if (abs(mini_tmp) gt 1000.0) or (abs(maxi_tmp) gt 1000.0) or        $
;     (abs(maxi_tmp) lt 0.1) then begin
;    maxs = string(maxi_tmp,format="(e8.2)")
;    mins = string(mini_tmp,format="(e9.2)")
;  endif else begin
;    maxs = string(maxi_tmp,format="(f6.2)")
;    mins = string(mini_tmp,format="(f7.2)")
;  endelse
;
;  xyouts, 2.1, mini_tmp, mins, charsize = 0.8
;  xyouts, 2.1, maxi_tmp, maxs, charsize = 0.8

endif

;-------------------------------------------------------------------------
; Next, lets label some times and iterations

figure_out_positions, 1, xreglimits, pos, ppp = 1

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

pos(1) = pos(1) - 0.1

;VarName = 'abs(Ey x b/B)'

xyouts, (pos(0) + pos(2))/2.0, pos(1), $
        "Plot of "+VarName+" in the X-"+yzs+" plane", $
        /norm, charsize=1.2, alignment = 0.5

spawn,'pwd',pwd
xyouts, 0.0, -0.02, pwd, /norm, charsize = 0.5

closedevice

end
