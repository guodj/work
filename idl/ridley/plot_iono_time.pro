
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

time = dblarr(nfiles/nskip + 1)
itime = intarr(6)

ntotal = 0

for n = 0, nfiles-1,nskip do begin

  swaptheta = 0
  filein = list(n)

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
        vars = strarr(nvars)

      endif

      if (strpos(mklower(line),"time") gt -1) then begin

        int_tmp = 0
        for i=0,5 do begin
          readf, 1, int_tmp
          itime(i) = int_tmp
        endfor

        c_a_to_r, itime, rtime
        time(ntotal) = rtime

      endif

      if (strpos(mklower(line),"northern") gt -1) then begin

        data = fltarr(nvars,nlons,nlats)
        for j=0,nlons-1 do for i=0,nlats-1 do begin
          readf,1,tmp
          data(*,j,i) = tmp
        endfor

	done = 1

      endif

      if (strpos(mklower(line),"variable") gt -1) then begin

        for i=0,nvars-1 do begin
          readf,1,line
          vars(i) = strmid(line,6,strlen(line)-6)
        endfor

      endif

      if (strpos(mklower(line),"all") gt -1) then begin

	nlons = nlons+1
	nlats = nlats/2
        data = fltarr(nvars,nlons,nlats)
        for j=0,nlons-2 do begin 
          for i=nlats-1,0,-1 do begin
            readf,1,tmp
            data(*,j,i) = tmp
	  endfor
	  for i=nlats-1,0,-1 do begin
	    readf,1,tmp
            data(*,j,i) = tmp
          endfor
        endfor
	swaptheta = 1

	data(*,nlons-1,*) = data(*,0,*)

      endif

    endwhile

    close,1

    if (n eq 0) then begin

      varlist = [5]
      nvars_to_plot = 1

      if (nvars_to_plot eq 0) then begin
	print, "Can not continue!"
	stop
      endif

      data_to_plot = fltarr(nfiles/nskip+1,nvars_to_plot,nlons,nlats)
      theta = fltarr(nlons,nlats)
      phi   = fltarr(nlons,nlats)

      nt = -1
      for i=0,nvars-1 do if strpos(mklower(vars(i)),'theta') gt -1 then nt = i
      np = -1
      for i=0,5 do if strpos(mklower(vars(i)),'psi') gt -1 then np = i
      if (nt eq -1) or (np eq -1) then begin
        print, "Can't file Theta or Psi variable. Please check file."
        stop
      endif

      theta(*,*) = reform(data(nt,*,*))
      phi(*,*)   = reform(data(np,*,*))

      if (swaptheta) then theta = 90.0 - theta

    endif

    data_to_plot(ntotal,*,*,*) = reform(data(varlist,*,*))

  endelse

  print, 'Finished Reading File '+filein
  ntotal = ntotal + 1

endfor

loc = where(time gt 0, ntimes)

if (ntimes gt 0) then begin

  for it = 0,ntimes-1 do begin
    if (it ne 7) then $
      data_to_plot(it,0,*,*) = data_to_plot(it,0,*,*) - data_to_plot(7,0,*,*)
  endfor
  data_to_plot(7,0,*,*) = data_to_plot(7,0,*,*) - data_to_plot(7,0,*,*)


  potential = fltarr(ntimes, nlons)
  latitudes = fltarr(ntimes, nlons)

  for it = 0,ntimes-1 do begin

    for il = 0,nlons-1 do begin

      potential(it,il) = max(data_to_plot(it,0,il,*))
      if (abs(min(data_to_plot(it,0,il,*))) gt potential(it,il)) then $
        potential(it,il) = min(data_to_plot(it,0,il,*))

      loc = where(data_to_plot(it,0,il,*) eq potential(it,il))
      latitudes(it,il) = theta(it,loc(0))

    endfor

  endfor

endif


setdevice, 'propagation.ps','p',4

maxi = max(abs(potential))
mini = -maxi

x = findgen(nlons)/(nlons-1) * 24.0 + 12.0

minx = fltarr(22)
miny = fltarr(22)
minz = fltarr(22)
maxx = fltarr(22)
maxy = fltarr(22)
maxz = fltarr(22)

ppp = 2
space = 0.01
pos_space, ppp, space, sizes, ny = 2
get_position,ppp,space,sizes,0,pos,/rect
pos(0) = pos(0) + 0.1
pos(2) = pos(2) - 0.05

plot, x, potential(8,*), $
      xstyle = 1, yrange = [mini,maxi], ystyle = 1, $
      xtickname = ['12','18','00','06','12'],xtickv=[12,18,24,30,36], $
      xticks = 4, xminor = 6, $
      xtitle = 'Magnetic Local Time', ytitle = 'Potential (kV)', $
      pos = pos

l = where(potential(8,*) eq min(potential(8,*)))
minx(0) = x(l(0))
minz(0) = latitudes(8,l(0))
miny(0) = min(potential(8,*))

l = where(potential(8,*) eq max(potential(8,*)))
maxx(0) = x(l(0))
maxz(0) = latitudes(8,l(0))
maxy(0) = max(potential(8,*))

ntimes = 30

for it = 9, ntimes-1 do begin
  oplot, x, potential(it,*)

l = where(potential(it,*) eq min(potential(it,*)))
minx(it-8) = x(l(0))
miny(it-8) = min(potential(it,*))
minz(it-8) = latitudes(it,l(0))

l = where(potential(it,*) eq max(potential(it,*)))
maxx(it-8) = x(l(0))
maxy(it-8) = max(potential(it,*))
maxz(it-8) = latitudes(it,l(0))

endfor

oplot, minx, miny
oplot, maxx, maxy

oplot, [12,36],[0.0,0.0]


closedevice

end
