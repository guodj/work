
mags = 0
nObs = 0

list = findfile("-t i*.idl")
if strlen(list(0)) gt 0 then filein = list(0) $
else filein = 'in000000.idl'

filein = ask('filename',filein)

list = findfile(filein)
nfiles = n_elements(list)
if nfiles eq 1 and strlen(list(0)) eq 0 then begin
  print, "I can't seem to find that file."
  stop
endif else begin

  if (nfiles gt 1) then begin
    print, "I found ",tostr(nfiles)," to plot."
    nskip = fix(ask('how many to skip between them','0'))+1
  endif else nskip = 1

endelse

line = ''

mr = float(ask('maximum range','40'))

if (nfiles gt 1) then begin

  relative = mklower(ask('universal time (U), or relative time (R)','R'))

  dt = 0.0
  if (strpos(relative,'r') eq 0) then begin
    dt = float(ask('offset time on first file (in seconds)','60.0'))
    ut = ''
  endif else ut = ' UT'

  if (strpos(relative,'r') ne 0) then ut = ' UT'

endif else begin

  ut = ''
  relative = 'u'
  dt = 0.0

endelse

plotsouth = mklower(ask('whether to plot southern hemisphere (Y or N)','N'))
isouth = 1
if (strpos(plotsouth,'n') eq 0) then isouth = 0

savefiles = 'n'
savefiles = mklower(ask('y to save indices',savefiles))

time = dblarr(nfiles/nskip + 1)
itime = intarr(6)

ntotal = 0

for n = 0, nfiles-1,nskip do begin

  swaptheta = 0
  filein = list(n)

  if strpos(filein,'save') gt -1 then begin

    restore, filein

    if (n eq 0) then begin
        for i=0,nvars-1 do begin
            if n eq 0 then print, chopr(' '+tostr(i),2),'. ',vars(i)
        endfor
    endif

    time(ntotal) = rtime

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
              vars = strarr(nvars)

          endif

          if (strpos(mklower(line),"variable") gt -1) then begin

              for i=0,nvars-1 do begin
                  readf,1,line
                  vars(i) = strmid(line,6,strlen(line)-6)
                  if (strpos(vars(i),'Jr [mA/m^2]') gt -1) then $
                     vars(i) = 'Jr [!Mm!XA/m!U2!N]'
                  if n eq 0 then print, chopr(' '+tostr(i),2),'. ',vars(i)
              endfor

          endif

          if (strpos(mklower(line),"time") gt -1 and $
              strpos(mklower(line),"simulation") lt 0) then begin

              int_tmp = 0
              for i=0,5 do begin
                  readf, 1, int_tmp
                  itime(i) = int_tmp
              endfor

              c_a_to_r, itime, rtime
              time(ntotal) = rtime

          endif

          if (strpos(mklower(line),"northern") gt -1) then begin

              data = fltarr(2,nvars,nlons,nlats)
              for j=0,nlons-1 do for i=0,nlats-1 do begin
                  readf,1,tmp
                  data(0,*,j,i) = tmp
              endfor

          endif

          if (strpos(mklower(line),"all") gt -1) then begin

              nlons = nlons+1
              nlats = nlats/2
              data = fltarr(2,nvars,nlons,nlats)
              for j=0,nlons-2 do begin 
                  for i=nlats-1,0,-1 do begin
                      readf,1,tmp
                      data(1,*,j,i) = tmp
                  endfor
                  for i=nlats-1,0,-1 do begin
                      readf,1,tmp
                      data(0,*,j,i) = tmp
                  endfor
              endfor
              swaptheta = 1

              data(*,*,nlons-1,*) = data(*,*,0,*)

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
          nt = -1
          for i=0,nvars-1 do $
            if strpos(mklower(vars(i)),'theta') gt -1 then nt = i
          np = -1
          for i=0,nvars-1 do if strpos(mklower(vars(i)),'psi') gt -1 then np = i
          if (nt eq -1) or (np eq -1) then begin
              print, "Can't file Theta or Psi variable. Please check file."
              stop
          endif

          theta = reform(data(*,nt,*,*))
          phi   = reform(data(*,np,*,*))

          if (swaptheta) then theta = 90.0 - theta


          ; Let's try to find the open/closed boundary:

          ocflb = 0

          ib_ = -1
          for i=0,nvars-1 do $
            if strpos(mklower(vars(i)),'1/b') gt -1 then ib_ = i
          
          if (ib_ gt -1) then begin
              ocflb = 1
              boundary = fltarr(2,nfiles/nskip+1,nlons,nlats)
          endif

      endif

  endelse

  if (n eq 0) then begin

      var = 0
      nvars_to_plot = 0
      while (var ge 0) do begin
        var = fix(ask('Variable Number to plot (-1 to exit)','-1'))
        if (var ge 0) and (var lt nvars) then begin
          if nvars_to_plot eq 0 then $
            varlist = [var]          $
          else varlist = [varlist,var]
          nvars_to_plot = nvars_to_plot + 1
        endif
      endwhile

      print, "You have selected "+tostr(nvars_to_plot)+" variables to plot." 

      if (nvars_to_plot eq 0) then begin
	print, "Can not continue!"
	stop
      endif

      data_to_plot = fltarr(2,nfiles/nskip+1,nvars_to_plot,nlons,nlats)

  endif

  data_to_plot(*,ntotal,*,*,*) = reform(data(*,varlist,*,*))

  if (ocflb) then $
    boundary(*,ntotal,*,*) = reform(data(*,ib_,*,*))

  print, 'Finished Reading File '+filein
  ntotal = ntotal + 1

endfor

toff = 0.0
if (strpos(relative,'r') eq 0) then toff = time(0) - dt - 24.0*3600.0 else toff = -dt
time = time(0:ntotal-1) - toff

if (mags) then begin

    file = 'mag_obs.txt'
    read_mag_obs, file, nObs, MagStats, MagLats, MagUts

endif

;
; We want to write an indices file, so figure out maximum & minimum for all
; Selected variables
;

indices = fltarr(2,ntotal,nvars_to_plot,2)
indices_vars = strarr(nvars_to_plot)

for i = 0, nvars_to_plot-1 do indices_vars(i) = vars(varlist(i))

for hem = 0, 1 do for i = 0, nvars_to_plot-1 do for j = 0, ntotal-1 do begin
  indices(hem,j,i,*) = mm(data_to_plot(hem,j,i,*,*))
endfor

if strpos(savefiles,'y') gt -1 then  $
  save, indices, indices_vars, time, file = list(0)+'.ind.save'

;
; Figure out plots per page...
;

space = 0.02

if ntotal eq 1 then begin
  ; if we have just 1 or two variables, just plot northern and southern
  ; hemisphere on the same page. If not, plot them on different pages:
  if nvars_to_plot le 2 then ppp = nvars_to_plot*2 else ppp = nvars_to_plot
endif else begin
  ; if we have less than 15 total images (per hemisphere), plot them all
  ; on the same page. If not, then plot each variable on a seperate page,
  ; with a limit of 9 plots per page.
  if (ntotal le 5) and (nvars_to_plot le 3) then begin
    ppp = ntotal * nvars_to_plot
    space = 0.01
  endif else begin
    if ntotal le 9 then ppp = ntotal else ppp = 9
  endelse
endelse

hems = ["Northern Hemisphere", "Southern Hemisphere"]

for hem = 0, isouth do begin

  print, "Working on ",hems(hem)

  if hem eq 0 then begin
    loc = where(reform(theta(0,0,*)) le mr)
    rang = reform(theta(0,*,loc))
    lons = reform(phi(0,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
  endif else begin
    loc = where(reform(theta(1,0,*)) gt 180.0-mr)
    rang = 180.0-reform(theta(1,*,loc))
    lons = reform(phi(0,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
  endelse

  for i = 0, nvars_to_plot-1 do begin

    if (ntotal eq 1) and (nvars_to_plot le 2) then begin
      if hem eq 0 and i eq 0 then begin
        pn = -1
        setdevice, list(0)+'.ps','p',5
	pos_space, ppp, space, sizes
      endif
    endif else begin

      if (ntotal le 5) and (nvars_to_plot le 3) then begin
        if hem eq 1 and i eq 0 then closedevice
        if i eq 0 then begin
          setdevice, list(0)+'_'+strmid(hems(hem),0,5)+'.ps','p',4
          pos_space, ppp, space, sizes, nx = nvars_to_plot
        endif
      endif else begin

        if ((i eq 0) and (ntotal eq 1)) or (ntotal gt 1) then begin
          pn = -1
          if i gt 0 then closedevice
          setdevice, list(0)+'_Var'+tostr(i)+'_'+strmid(hems(hem),0,5)+'.ps',$
                     'p',4
          pos_space, ppp, space, sizes
        endif

      endelse

    endelse

    print, "Working on variable ",vars(varlist(i))

    mini = min(data_to_plot(hem,*,i,*,loc))
    maxi = max(data_to_plot(hem,*,i,*,loc))

    if (maxi eq mini) then maxi = mini + 1.0
    if (mini lt 0.0) then begin
        maxi = max([abs(mini),maxi])
;        if (maxi gt 8.0) then maxi = (fix(maxi/8)+1)*8.0 $
;        else begin
            lo = fix(alog10(maxi))
            maxi = maxi/(10.0^(lo-1))
            maxi = (fix(maxi/8)+1)*8.0
            maxi = maxi*(10.0^(lo-1))
;        endelse
    endif
    maxi = float(ask('maximum value to plot',string(maxi)))

    if (mini ge 0.0) then begin
      levels = maxi*findgen(9)/8.0 + mini
      makect, 'all'
      ncolors = 255
    endif else begin

      mini = -maxi
      nl = 8
      nld2 = nl/2
      levels = fltarr(nl)
      levels(0:nld2-1) = -maxi*(nld2-findgen(nld2))/float(nld2)
      levels(nld2:nl-1) = maxi*(findgen(nld2)+1)/float(nld2)
      makect, 'mid'
      ncolors = 255

    endelse

    c_levels = (maxi-mini)*findgen(30)/29.0 + mini
    c_colors = (ncolors-1)*findgen(30)/29.0 + 1

    for n = 0, ntotal - 1 do begin

      if (ntotal le 5) and (nvars_to_plot le 3) then begin
        if (ntotal eq 1) and (nvars_to_plot le 2) then pn = (pn + 1) mod ppp $
        else pn = n*nvars_to_plot + i
      endif else begin 
        pn = (pn + 1) mod ppp
      endelse

      if pn eq 0 then plotdumb
      get_position, ppp, space, sizes, pn, pos
      pos(1) = pos(1) + 0.5*space*n
      pos(3) = pos(3) + 0.5*space*n

      contour, data_to_plot(hem,n,i,*,loc), xpos, ypos,        $
	/follow, nlevels=30, /noerase, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,/cell_fill

      contour, data_to_plot(hem,n,i,*,loc), xpos, ypos,        $
        /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase, $
	c_linestyle = 3.0*(levels lt 0.0), c_thick = 2

      if (ocflb) then $
        contour, boundary(hem,n,*,loc), xpos, ypos,        $
        pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase, $
	levels = [-1,0], c_thick = 4

      if (mags and nObs gt 0) then begin

          c_r_to_a, itime, time(n)
          uthr = itime(3)*60.0 + itime(4)
          for iOb = 0, nObs-1 do begin

              rm = 90.0 + (hem*2-1) * MagLats(iOb)
              tm = (uthr - MagUts(iOb)) * !pi/(12*60.0) - !pi/2 

              if (rm lt mr) then begin

                  xm = rm*cos(tm)
                  ym = rm*sin(tm)

                  oplot, [xm], [ym], psym = 4

                  xyouts, xm, ym, magstats(iOb), charsize = 0.8

              endif

          endfor
          
      endif

;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

      no00 = 1
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
        maxs = string(maxi_tmp,format="(f6.2)")
        mins = string(mini_tmp,format="(f7.2)")
      endelse

      charsize = 1.0
      if ppp gt 9 then charsize = 0.75

      xyouts, pos(0),pos(1), mins, /norm, charsize = charsize
      xyouts, pos(2),pos(1), maxs, /norm, align=1.0, charsize = charsize

      if (ntotal eq 1) then begin
        xyouts, pos(2),pos(3),vars(varlist(i)), 	$
	  align=1.0, /norm
      endif

      if (n eq ntotal-1) then begin
         ctpos = pos
         ctpos(0) = ctpos(0)+0.02
         ctpos(2) = ctpos(2)-0.02
         ctpos(1) = ctpos(1)-0.02
         ctpos(3) = ctpos(1)+0.01
         plotct, 255, ctpos, mm(levels), vars(varlist(i)), /bottom
      endif

      if (ntotal le 5) and (nvars_to_plot le 3) then begin
        if i eq 0 then begin
          c_r_to_a, itime, time(n)
          c_a_to_s, itime, stime
          xyouts, pos(0)-0.03, (pos(1)+pos(3))/2.0, strmid(stime,0,15), $
                  alignment = 0.5, /norm, orient = 90
        endif
        if n eq 0 then begin
;          xyouts, (pos(0)+pos(2))/2.0, pos(3)+2*space, vars(varlist(i)), $
;                  alignment = 0.5, /norm
          if pn eq 0 then begin
            xyouts, 0.5, pos(3)+4*space, hems(hem), alignment = 0.5, $
                    /norm, charsize = 1.25
            xyouts, -0.01, -0.01, filein, /norm, charsize = 0.5, orient = 90
          endif
        endif
      endif else begin

	p1 = pos(0) + (pos(2) - pos(0))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.95
	p2 = pos(3) - (pos(3) - pos(1))*0.50 * (1.0 - sin(45.0*!pi/180.0))*0.95

        c_r_to_a, itime, time(n)
        c_a_to_s, itime, stime
        xyouts, p1, p2, $
                strmid(stime,10,2)+strmid(stime,13,5)+ut, $
                /norm, alignment = 0.5, charsize = 0.8, orient = 45.0

        if pn eq 0 then begin
          if (ntotal gt 1) then begin
            xyouts, 0.0, 1.01, hems(hem), alignment = 0.0, /norm
            xyouts, 1.0, 1.01, vars(varlist(i)), alignment = 1.0, /norm, $
                  charsize = 1.1
            if (strpos(relative,'u') eq 0) then $
               xyouts, 0.5, 1.01, strmid(stime,0,9), /norm, alignment = 0.5
          endif else begin
            xyouts, pos(0), pos(3)+space, hems(hem), $
                    alignment = 0.0, /norm, charsize = 1.25
            c_r_to_a, itime, time(n)
            c_a_to_s, itime, stime
            xyouts, 1.0-pos(0), pos(3)+space, strmid(stime,0,15), $
                    alignment = 1.0, /norm, charsize = 1.25
          endelse
          xyouts, -0.01, -0.01, filein, /norm, charsize = 0.5
        endif

      endelse

    endfor

  endfor

endfor

closedevice

end    
