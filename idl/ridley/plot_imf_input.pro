
pro plot_data, data, time, nvars, titles, overplots

  if n_elements(overplots) eq 0 then overplots = intarr(n_elements(titles))

  stime = min(time)
  c_r_to_a, itime, stime
  if (itime(3) le 1) then begin
      itime(3) = 0
      itime(4) = 0
      itime(5) = 0
      c_a_to_r, itime, stime
  endif

  etime = max(time)
  c_r_to_a, itime, etime
  if (itime(3) ge 22) then begin
      itime(3) = 24
      itime(4) = 0
      itime(5) = 0
      c_a_to_r, itime, etime
  endif

  tmptime = time - stime

  time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

  loc = where(overplots eq 0, ppp)
  if ppp eq 0 then ppp = 1
  space = 0.01
  pos_space, ppp, space, sizes, ny = ppp

  plotnum = 0

  for nv = 0,nvars-1 do begin

    if plotnum eq 0 then plotdumb

    get_position, ppp, space, sizes, plotnum, pos, /rect

    pos(0) = pos(0) + 0.1

    if nv lt nvars-1 then begin
      xtn = strarr(60)+' '
      xt  = ' '
    endif else begin
      xtn = xtickname
      xt  = xtitle
    endelse

    loc = where(data(nv,*) gt -1.0e10,count)
    if count gt 0 then range = mm(data(nv,loc)) else range = [-1.0,1.0]

    if (nv lt nvars-1) then begin
        if (overplots(nv+1)) then begin
            loc = where(data(nv+1,*) gt -1.0e10,count)
            if count gt 0 then begin
                range(0) = min([min(data(nv+1,loc)),range(0)])
                range(1) = max([max(data(nv+1,loc)),range(1)])
            endif 
        endif
    endif


    if (range(0) eq range(1)) then begin
      if range(0) eq 0.0 then begin
        range(0) = -1.0
        range(1) = 1.0
      endif else begin
        range(0) = range(0) - 0.1*range(0)
        range(1) = range(1) + 0.1*range(1)
      endelse
    endif else begin
      r = range(1) - range(0)
      if (range(0) gt 0 and range(0)-r*0.1 gt 0) then $
         range(0) = range(0) - r*0.1
      range(1) = range(1) + r*0.1
    endelse

    if (overplots(nv) eq 0) or (nv eq 0) then begin

      plot, tmptime, data(nv,*), xstyle=1,		$
	  ytitle = titles(plotnum),		$
	  xtickname = xtn,			$
	  xtitle = xt,			$
	  xtickv = xtickv,			$
	  xminor = xminor,			$
	  xticks = xtickn,   $
          pos = pos, /noerase, $
          min_val = -1.0e10, yrange = range, $
          thick = 3

      oplot, mm(tmptime),[0.0,0.0],linestyle = 1

      ndays = fix(max(tmptime)/(24.0*3600.0))

      for i=1,ndays do oplot, 24.0*3600.0*[float(i),float(i)],   $
                          [-max(data(nv,*))*10.0,max(data(nv,*))*10.0], $
                          linestyle = 1

      plotnum = plotnum + 1

    endif else begin

      oplot, tmptime, data(nv,*), linestyle = 2, thick = 3

    endelse


  endfor

  return

end

if (n_elements(file) eq 0) then begin
  filelist = findfile('IMF*.dat')
  file = filelist(0)
endif
file = ask('filename to plot',file)

read_imf_input, file, time, mag, vel, den, temp, nPts

;; file2 = "imf20031029-30.log"
;; read_imf_input, file2, time2, mag2, vel2, den2, temp2, nPts2

;;;  smooth = '60'
;;;  smooth = fix(ask('number of points for smoothing',smooth))
;;;  
;;;  mag_smooth = mag
;;;  vel_smooth = vel
;;;  den_smooth = den
;;;  temp_smooth = temp
;;;  
;;;  for i=0,npts-1 do begin
;;;  
;;;    if i-smooth/2 ge 0 then istart = i-smooth/2 else istart = 0
;;;    if i+smooth/2 le npts-1 then iend = i+smooth/2 else iend = npts-1
;;;  
;;;    for j=0,2 do begin
;;;      mag_smooth(j,i) = mean(mag_smooth(j,istart:iend))
;;;      vel_smooth(j,i) = mean(vel_smooth(j,istart:iend))
;;;    endfor
;;;    den_smooth(i) = mean(den_smooth(istart:iend))
;;;    temp_smooth(i) = mean(temp_smooth(istart:iend))
;;;  
;;;  endfor

; Plot magnetic fields on first page:

setdevice,'imf.p1.ps','p',4

titles  = ['Bx (nT)','By GSM (nT)','Bz GSM (nT)','Density (/cc)','Vx (km/s)']
titles2 = ['Vx (nT)','Vy GSM (nT)','Vz GSM (nT)','T (K)','P (nP)']
titles3 = ['P ram (nP)','P kin (nP)','P mag (nP)','ratio (Pb/Pr)','Stand Off (Re)']

nvars = 5

data = fltarr(nvars,npts)
data(0,*) = mag(0,*)
data(1,*) = mag(1,*)
data(2,*) = mag(2,*)
data(3,*) = den(*)
data(4,*) = vel(0,*)

overplots = [0,0,0,0,0]

plot_data, data, time, nvars, titles, overplots

closedevice

setdevice,'imf.p2.ps','p',4

mu0 = 4.0*!pi*1.0e-7
k   = 1.3807e-23
mp  = 1.6726e-27

b_sw = sqrt(mag(0,*)^2+mag(1,*)^2+mag(2,*)^2) / 1.0e9

p_kin = (den * 1e6) * temp * 1.3807e-23 * 1.0e9
p_ram = (den*1e6) * mp * (vel(0,*)*1000.0)^2 * 1.0e9
p_b   = (b_sw)^2/(2*mu0) * 1.0e9

nvars = 5

data = fltarr(nvars,npts)
data(0,*) = vel(0,*)
data(1,*) = vel(1,*)
data(2,*) = vel(2,*)
data(3,*) = temp(*)
data(4,*) = p_ram(*) + p_b(*) + p_kin(*)

overplots = [0,0,0,0,0]

plot_data, data, time, nvars, titles2, overplots

closedevice

nvars = 5

data = fltarr(nvars,npts)

setdevice,'imf.p3.ps','p',4

data(0,*) = p_ram(*)
data(1,*) = p_kin(*)
data(2,*) = p_b(*)
data(3,*) = p_b(*)/p_ram(*)

p = p_ram(*) + p_kin(*) + p_b(*)

beq = 33000.0e-9
r = findgen(3000)/200.0 + 1.0
b = beq/r^3
p_b_earth = (2.0*b)^2/(2*mu0) * 1.0e9

for i=0,npts-1 do begin

  dist = abs(p_b_earth -  p(i))
  loc = where(dist eq min(dist))
  data(4,i) = r(loc(0))

endfor

plot_data, data, time, nvars, titles3, overplots

closedevice

end

