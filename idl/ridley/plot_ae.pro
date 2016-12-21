
filein = findfile('ae07.final')

nfiles = long(n_elements(filein))

dst = fltarr(1440*nfiles)
ae = fltarr(1440*nfiles)
t = findgen(1440*nfiles)*60.0
dst_p = fltarr(1440*nfiles)
stime = intarr(6)
etime = intarr(6)

for j=0L,nfiles-1 do for i=0L,1439 do begin

  if (i eq 0) then begin

    close,1
    openr,1,filein(j)

    line = ''
    readf,1,line
    readf,1,line
    readf,1,line
    readf,1,line

  endif

  readf,1,line

  if (j eq 0 and i eq 0) then 			$
    for k=0,4 do stime(k) = fix(strmid(line,k*3,3))

  if (j eq nfiles-1 and i eq 1439) then 			$
    for k=0,4 do etime(k) = fix(strmid(line,k*3,3))

  ae(j*1440+i) = float(strmid(line,22,6))
  dst(j*1440+i) = float(strmid(line,62,6))
  dst_p(j*1440+i) = float(strmid(line,72,6))

endfor

etime(4) = etime(4) + 1

close,1

mmdst = [-150.0,50.0]

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

plot, t, dst, xstyle=1, ystyle=1,		$
	yrange = mmdst, ytitle = 'DST',		$
	xtickname = xtickname,			$
	xtitle = xtitle,			$
	xtickv = xtickv,			$
	xminor = xminor,			$
	xticks = xtickn

oplot, t,dst_p, linestyle=2
oplot, [btr,etr],[0.0,0.0], linestyle = 1

plot, t, ae, xstyle=1,				$
	ytitle = 'AE',				$
	xtickname = xtickname,			$
	xtitle = xtitle,			$
	xtickv = xtickv,			$
	xminor = xminor,			$
	xticks = xtickn

end

