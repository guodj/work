
print, 'This program assumes that HPI data is in flatfile format.'
print, ''

dirin = '/d/ridley/d.data/d.hpi/d.flat/'
dirin = ask('directory to find flatfiles',dirin)
if (strmid(dirin,strlen(dirin)-1,1) ne '/') then dirin=dirin+'/'

sdate = '28-dec-92'
sdate = ask('starting date',sdate)

edate = '28-dec-92'
edate = ask('ending date (can be same as starting)',sdate)

fileout = 'power.ascii'
fileout = ask('output file name',fileout)

c_s_to_a, stime, sdate
c_a_to_r, stime, srtime

c_s_to_a, etime, edate
c_a_to_r, etime, ertime

if (ertime lt srtime+24.0*3600.0) then begin
  ertime = srtime+24.0*3600.0
  c_r_to_a, etime, ertime
endif

if stime(0) gt 1900 then begin
  stime(0) = stime(0) - 1900
  if stime(0) gt 100 then stime(0) = stime(0) - 100
endif

if etime(0) gt 1900 then begin
  etime(0) = etime(0) - 1900
  if etime(0) gt 100 then etime(0) = etime(0) - 100
endif

files = ''

for yr = stime(0), etime(0) do begin

  fname = dirin+'hpi'+chopr('0'+tostr(yr),2)+'*.hed'
  if (yr eq stime(0)) then files = findfile(fname)			$
  else files = [files,findfile(fname)]

endfor

if strlen(files(0)) gt 0 then begin

  nfiles = n_elements(files)

  first = 1

  for i=0,nfiles-1 do begin

    read_flat_scalor, stime, etime, [0], time, data,	$
                  rows, filename = strmid(files(i),0,strlen(files(i))-4)

    n_e = n_elements(data(0,*))

    if n_e gt 1 then begin

      if first then begin
	savedata = fltarr(n_e)
	savetime = dblarr(n_e)
        savedata(0:n_e-1) = data(0,0:n_e-1)
        savetime(0:n_e-1) = time(0,0:n_e-1)
        first = 0
      endif else begin
	tdata = fltarr(n_e)
	ttime = dblarr(n_e)
        tdata(0:n_e-1) = data(0,0:n_e-1)
        ttime(0:n_e-1) = time(0,0:n_e-1)
        savedata = [savedata,tdata]
        savetime = [savetime,ttime]
      endelse

    endif

  endfor

  sort_a,savetime,savedata

  plot,savetime,savedata

  range = [srtime, ertime]
  nbin = n_elements(savetime)

  bin_y, savedata, savetime, nbin, output, smooth=5, range=range,	$
       missing = 0.0

  outtime = double(findgen(nbin))*(range(1)-range(0))/double(nbin) + range(0)

  oplot, outtime, output, linestyle = 1

  openw,1,fileout

  n = 0

  for i=0, nbin-1 do begin

    if i eq 0 then line = 'POWER = '

    if output(i) gt 0.0 then begin

      if n ne 0 then line = line + ', '

      c_r_to_a, itime, outtime(i)
      for j=2,4 do line = line + string(itime(j),format='(i2)')+','

      line = line + string(output(i),format='(f5.1)')

      n = n + 1

      if (n eq 4) or (i eq nbin-1) then begin
        printf,1,line
        line = ' '
        n = 0
      endif

    endif

  endfor

  close,1

endif else print, 'No files found!'

end