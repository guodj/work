
common ffinfo, header

dir = "/Despiked/Averaged"

filelist = findfile('??'+dir)
nfiles = n_elements(filelist)

filelist_hed = strarr(nfiles)
nf = 0
for i=0,nfiles-1 do begin
  if (strpos(filelist(i),'.hed') gt 0) then begin
    filelist_hed(nf) = filelist(i)
    nf = nf + 1
  endif
endfor

filelist = filelist_hed(0:nf-1)
nfiles = nf

filelist_master = filelist

statlist = [strmid(filelist(0),0,3)]
nstats = 1

for i=1,nfiles-1 do begin

  stat = strmid(filelist(i),0,3)

  ifound = 0
  istat = 0

  while (not ifound and istat lt nstats) do begin
    if (stat eq statlist(istat)) then ifound = 1
    istat = istat + 1
  endwhile

  if (not ifound) then begin
    statlist = [statlist,stat]
    nstats = nstats+1
  endif

endfor

; We need to get the latitude of the stations - so we have to read in the
; master.psi file

openr,1,'master.psi'
ns = 0
line = ''
readf,1,ns
readf,1,line
stations = strarr(ns)
lats     = fltarr(ns)
for i=0,ns-1 do begin
  readf,1,line
  stations(i) = mklower(strmid(line,4,3))
  lats(i)     = float(strmid(line,36,6))
endfor
close,1

components = strarr(nstats)

for istat = 0,nstats-1 do begin

  print, 'working on station : ',statlist(istat),istat

  loc = where(stations eq statlist(istat),count)

  if count gt 0 then latitude = lats(loc(0)) else begin
    print, "Station not found!!!"
    latitude = 0.0
  endelse

  print, "Station latitude : ", latitude

  n = 0
  filelist = strarr(31)
  for i=0,nfiles-1 do begin
    if (strpos(filelist_master(i),statlist(istat)) ge 0) then begin
      filelist(n) = filelist_master(i)
      n = n + 1
    endif
  endfor

  filelist = filelist(0:n-1)

;  filelist = findfile('*'+statlist(istat)+'.amie')

  ndays = n_elements(filelist)

  nptsmax = 1440

  data    = fltarr(3,nptsmax,ndays+1)
  time    = dblarr(nptsmax,ndays)
  dataqdr = fltarr(3,nptsmax,ndays) - 99999.0
  
  corr = fltarr(ndays,ndays)
  cave = fltarr(ndays)
  
  itime = intarr(6)

  for n=0,ndays-1 do begin
  
    print, filelist(n)

    itime(0) = fix(strmid(filelist(n),3,2))
    itime(1) = fix(strmid(filelist(n),5,2))
    itime(2) = fix(strmid(filelist(n),7,2))
    itime(3) = 0
    itime(4) = 0
    itime(5) = 0

    c_a_to_r, itime, stime

    itime(2) = itime(2) + 1
    c_a_to_r, itime, etime

    col_scal = [0,1,2]

    phed = strpos(filelist(n),".hed")
    filename = strmid(filelist(n),7,2)+dir+'/'+ $
               strmid(filelist(n),0,phed)

    read_flat_scalor, stime, etime, col_scal, time1f, data1f, nrows,     $
        filename = filename

    if (nrows(0) ne nptsmax) then begin
      time(*,n) = stime + findgen(1440)*60.0
      for i=0,2 do data(i,*,n) = -1.0e32
      if (nrows(0) gt 0) then begin
          time(0:nrows(0)-1,n) = time1f(0,0:nrows(0)-1)
          for i=0,2 do data(i,0:nrows(i)-1,n) = data1f(i,0:nrows(i)-1)
      endif
    endif else begin
      time(*,n) = time1f(0,*)
      for i=0,2 do data(i,*,n) = data1f(i,*)
    endelse

    if (n eq 0) then $
      components(istat) = strmid(header.na(0),3,1) + $
                          strmid(header.na(1),3,1) + $
                          strmid(header.na(2),3,1)

  endfor
  
  ; First thing we have to do is search for days in which the values
  ; are all the same...  Let's search for 1 hour in which the data
  ; is exactly the same, then go back and forth from there to blank out
  ; all of the bad values.
  
  for n1 = 0,ndays-1 do begin
  
    for i=0,23 do begin
  
      for j=0,2 do begin
  
  	d = data(j,i*60:(i+1)*60-1,n1)
  	loc = where(d gt -90000.0,count)
  	if count gt 0 then sd = stddev(d) else sd = 1.0
  
  	if sd eq 0.0 then begin
  
  	  print, n1,i,j
  	  data(j,i*60:(i+1)*60-1,n1) = -99999.0
  
  	endif
  
      endfor
  
    endfor
  
  endfor
  
  
  for n1 = 0,ndays-1 do begin
  
    for n2 = 0,ndays-1 do begin
  
      if n1 eq n2 then begin
  	corr(n1,n2) = 1.0
      endif else begin
  
  	corr(n1,n2) = 0.0
  
  	for i=0,2 do begin
  
  	  d1 = data(i,*,n1)
  	  d2 = data(i,*,n2)
  
  	  loc = where(d1 gt -90000.0 and d2 gt -90000.0,count)
  
  	  if count lt 60 then begin
  	    corr(n1,n2) = 0.0
  	  endif else begin
  	    corr(n1,n2) = corr(n1,n2) + correlate(d1(loc),d2(loc))/3.0
  	  endelse
  
  	endfor
  
      endelse
  
    endfor
  
    cave(n1) = mean(corr(n1,*))
  
  endfor
  
  cutoff = median(cave)
  if ndays eq 1 then cutoff = 0.5*cutoff
  
  nave = nptsmax/8   ; this should be around 3 hours...
  
  ; For low latitude stations, we need to average a whole day together, then do
  ; 3 hour running average...
  
  for n1 = 0, ndays-1 do begin
  
    ndlook = 5
    count = 0
    nmin = 2
    if (nmin * 2 gt ndays) then nmin = ndays/2
  
    while count lt nmin do begin
  
      istart = max([0,n1-ndlook])  ; look 5 days back
      iend   = istart + ndlook*2   ; look for a total of 11 days (5 back, 5 forward)
  
      if (iend gt ndays-1) then begin
  	iend = ndays-1
  	if (iend-ndlook*2 lt 0) then istart = 0 else istart = iend-ndlook*2
      endif
  
      loc = where(cave(istart:iend) ge cutoff,count)
  
      if count lt nmin and ndays gt 2 then begin
  	ndlook = ndlook + 1
  	print, "Increasing ndlook : ", ndlook
      endif else begin
  	if ndays eq 1 then count = nmin
      endelse
  
      if (ndlook gt ndays and count lt nmin) then begin
  	  print, "Something is horribly wrong here."
  	  stop
      endif
  
    endwhile
  
    if ndays eq 1 then begin
      count = 1
      loc = [0]
      istart = 0
      iend = 0
    endif
  
    if (count eq 0) then begin      ; we didn't find any good days 
  				    ; within the 10 day window - opps.
  
      print, "No good days found!!!"
      stop
  
    endif else begin
  
      quietday  = fltarr(3,nptsmax+nave*2)
      smoothday = fltarr(3,nptsmax)
      quietnpts = intarr(3,nptsmax+nave*2)*0.0
  
      for iday = 0,count-1 do begin
  	i = istart + loc(iday)
  
  	for icomp = 0,2 do for itime = 0, nptsmax-1 do begin
  	  if (data(icomp, itime, i) gt -90000.0) then begin
  	    quietday(icomp,itime+nave) = quietday(icomp,itime+nave) + $
  					 data(icomp, itime, i)
  	    quietnpts(icomp,itime+nave) = quietnpts(icomp,itime+nave) + 1
  	  endif
  	endfor
      endfor  
  
      loc = where(quietnpts eq 0,count)
      quietday(loc) = -99999.0
  
      loc = where(quietnpts gt 0,count)
      if (count gt 0) then begin
  	quietday(loc) = quietday(loc)/quietnpts(loc)
  
  	for icomp = 0,2 do begin
  	  quietday(icomp,0:nave-1) = quietday(icomp,nptsmax : nptsmax+nave-1)
  	  quietday(icomp,nptsmax+nave:nptsmax+2*nave-1) = quietday(icomp, nave :2*nave-1)
  	endfor
  
  	for i=0,nptsmax-1 do begin
  
  	  is = nave + i - nave/2
  	  ie = is+nave
  
  	  for icomp = 0,2 do begin
  	    loc = where(quietday(icomp,is:ie) gt -90000,count)
  	    if (count gt 0) then begin
  	      smoothday(icomp,i) = mean(quietday(icomp,is+loc))
  	    endif else begin
  	      smoothday(icomp,i) = quietday(icomp,is)
  	    endelse
  	  endfor
  
  	endfor
  
  ; Fill in any missing points in the smoothed day
  
  	for icomp = 0,2 do begin
  	  loc = where(smoothday(icomp,*) lt -90000.0,count)
  	  if (loc(0) eq 0) then begin
  	      loc2 = where(smoothday(icomp,*) gt -90000.0)
  	      if (loc2(0) ne -1) then begin
  		smoothday(icomp,0) = smoothday(icomp,loc2(0))
  	      endif else smoothday(icomp,0) = 0.0
  	  endif
  	  if (count gt 0) then begin
  	    for i=1,count-1 do $
  	      smoothday(icomp,loc(i)) = smoothday(icomp,loc(i)-1) 
  	  endif
  	endfor
  
  ; If the station is above a certain latitude, then replace the smoothday
  ; with the median of the smoothday:
  
  	if (latitude gt 55.0) then begin
  	  for icomp = 0,2 do smoothday(icomp,*) = median(smoothday(icomp,*))
  	endif
  
  ; Now subtract that day from the rest of the day
  
  	for icomp = 0,2 do begin
  	  loc = where(data(icomp,*,n1) gt -90000.0)
  	  if (loc(0) gt -1) then begin
  	     dataqdr(icomp,loc,n1) = data(icomp,loc,n1) - smoothday(icomp,loc)
  	  endif
  	endfor
  
      endif else begin
  
  	print, "Yikes - no points found!!"
  	stop
  
      endelse
  
    endelse
  
  endfor
  
  
  ; Now we need to try to fix times in which there is a baseline jump.
  ; We do this by looking at the median of the day and seeing if that is
  ; close to zero.  We quantify "close" as the standard deviation over the
  ; day times 2.  We also compare the median of the day to the median of
  ; the day before or the day after.  If the median of the current day
  ; is further away from zero than the other day, and it is further
  ; away than the standard deviation time 2, then we assume that there
  ; was a baseline shift, and we remove the median.
  
  for n1 = 0, ndays-1 do begin
  
    for icomp = 0,2 do begin
      loc = where(dataqdr(icomp,*,n1) gt -90000.0,count)
      if (count gt 2) then begin
  	m = median(dataqdr(icomp,loc,n1))
  	s = stddev(dataqdr(icomp,loc,n1))
      endif else begin
  	m = -99999.0
  	s = 0.0
      endelse
      if (abs(m) gt 2.0*s) then begin
  	print, "median > 2*stddev : ",m,s,n1,icomp
  
  	m1 = -99999.0
  	if (n1 gt 0) then begin
  	  loc2 = where(dataqdr(icomp,*,n1-1) gt -90000.0,count)
  	  if count gt 0 then m1 = median(dataqdr(icomp,loc2,n1-1))
  	endif
  
  	if (n1 lt ndays-1 and m1 lt -90000.0) then begin
  	  loc2 = where(dataqdr(icomp,*,n1+1) gt -90000.0,count)
  	  if count gt 0 then m1 = median(dataqdr(icomp,loc2,n1+1))
  	endif
  
  	if (m1 gt -90000.0) then begin
  	  if (abs(m1-m) gt 2.0*s and abs(m) gt abs(m1)) then begin
  	    print, "m1-m > 2*s - This is a real bad day :",m1,m,s
  	    if (loc(0) gt -1) then begin
  	      dataqdr(icomp,loc,n1) = dataqdr(icomp,loc,n1) - m
  	    endif
  	  endif else print, "m1-m < 2*s :",m1,m,s
  
  	endif
      endif
    endfor
  
  endfor
  
  if !d.name eq 'X' then setdevice,'stats.ps','l',4,0.95
  
  ppp = 3
  space = 0.01
  pos_space, ppp, space, sizes, ny = ppp
  
  get_position,ppp,space,sizes,0,pos1,/rect
  get_position,ppp,space,sizes,1,pos2,/rect
  get_position,ppp,space,sizes,2,pos3,/rect
  pos1(0) = pos1(0) + 0.05
  pos2(0) = pos2(0) + 0.05
  pos3(0) = pos3(0) + 0.05
  
  pos = fltarr(3,4)
  pos(0,*) = pos1
  pos(1,*) = pos2
  pos(2,*) = pos3
  
  stime = min(time)
  etime = max(time)
  
  time_axis, stime, etime, s_time_range, e_time_range,        	$
  	  xtickname, xtitle, xtickvalue, xminor, xtickn
  
  comp = ['X','Y','Z']
  
  i = istat

  plotdumb

  for icomp = 0,2 do begin

    if (icomp eq 2) then begin
      xtn = xtickname
      xt  = xtitle
    endif else begin
      xtn = strarr(20)+' '
      xt  = ' '
    endelse

    title = ' '
    if (icomp eq 0) then title =statlist(i)+' '+tostr(i+1)+' of '+tostr(nstats)

    plot, time-stime,dataqdr(icomp,*,*), pos = pos(icomp,*),	 $
	xstyle = 1, /noerase,		$
	xtickname = xtn, xtickv=xtickvalue, 			$
	xticks = xtickn, xminor = xminor, xtitle = xt,		$
	xrange = [s_time_range, e_time_range],                  $
        ytitle = strmid(components(i),icomp,1), $
        yrange = yrange, min_val = -5000.0, $
        title = title

    for ifile = 1, ndays-1 do begin
       it = long(nptsmax)*long(ifile)
       oplot, [time(it), time(it)]-stime,[-5000,5000], linestyle = 1
    endfor

    oplot, [s_time_range, e_time_range], [0.0,0.0], linestyle = 1

    loc = where(dataqdr(icomp,*,*) gt -5000.0,count)
    if (count gt 0) then m = median(dataqdr(icomp,loc)) else m = -99999.

    ms = strcompress(string(m))

    xyouts, pos(icomp,2)+0.02, pos(icomp,1), ms, orient = 90, /norm

  endfor

  for n1 = 0, ndays-1 do begin
    c_r_to_a, itime, time(0,n1)
    ymd = chopr('0'+tostr(itime(0)),2)+ $
          chopr('0'+tostr(itime(1)),2)+ $
          chopr('0'+tostr(itime(2)),2)

    filename = ymd+'/'+mklower(statlist(i))+ymd+'.remca'
    openw,1,filename, error = iError

    if (iError ne 0) then begin
      spawn, "mkdir "+ymd
      openw,1,filename
    endif


    for it=0,nptsmax-1 do begin
      c_r_to_a, itime, time(it,n1)
      printf,1,format='(5i2,a,1x,a,3f10.1)', $
               itime(0) mod 100, itime(1), itime(2), itime(3), itime(4),$
               statlist(i), components(i), $
               dataqdr(0,it,n1), dataqdr(1,it,n1), dataqdr(2,it,n1)
    endfor
    close,1
  endfor

endfor

closedevice

end
