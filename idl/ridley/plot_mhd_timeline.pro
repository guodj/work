
xpoints = 15.0 - findgen(31)
ypoints = fltarr(31) + 0.0
;xpoints = [10.0, 5.0,  0.0,  -5.0, -10.0]
;ypoints = [ 5.0, 5.0,  5.0,   5.0,   5.0]

npoints = n_elements(xpoints)

if (n_elements(plane) eq 0) then plane = 'y'
plane = mklower(ask('which plane (z/y)',plane))
if (strpos(plane,'y') eq 0) then IsYPlane = 1 else IsYPlane = 0

if (IsYPlane) then filelist = findfile('IO2/y=0*.out') $
else filelist = findfile('IO2/z=0*.out')

nfiles = n_elements(filelist)

times = dblarr(nfiles)
iterations = lonarr(nfiles)

for ifile=0,nfiles-1 do begin

  filename = filelist(ifile)

  print, 'Reading File ',filename

  ; Read data from file

  gettype,filename,filetype,npict

  openfile,10,filename,filetype
  gethead,10,filetype,headline,physics,it,time,gencoord, $
          ndim,neqpar,nw,nx,eqpar,variables,pictsize=pictsize
  close,10

  phys=physics
  times(ifile) = time
  iterations(ifile) = it

  openfile,10,filename,filetype

  getpict,10,filetype,npict,x,w,headline,phys,it,time,$
          gencoord,ndim,neqpar,nw,nx,eqpar,variables,error

  if (ifile eq 0) then data = fltarr(nfiles, npoints, nw)

  for ipoint = 0, npoints-1 do begin

    r2 = (reform(x(*,0,0))-xpoints(ipoint))^2 + $
         (reform(x(*,0,1))-ypoints(ipoint))^2

    loc = where(r2 eq min(r2))
    loc = loc(0)

    if (r2(loc) ne 0.0) then begin
      xpoints(ipoint) = x(loc,0,0)
      ypoints(ipoint) = x(loc,0,1)
    endif

    data(ifile, ipoint, *) = w(loc, 0, *)

  endfor

endfor

ivar = 1

v = reform(sqrt(data(*,*,1)^2+data(*,*,2)^2+data(*,*,3)^2))
range = fltarr(npoints)

for ipoint = 0, npoints-1 do begin

  range(ipoint) = max(v(*, ipoint)) - min(v(*, ipoint))

  v(*,ipoint) = (v(*,ipoint)-min(v(*,ipoint)))/range(ipoint)

endfor

times = times/60.0

last = max(times)

mstart = 

end
