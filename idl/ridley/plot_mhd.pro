
transform = 'regular'
gencoord = 1
doask = 0
plotgrid = 1

files = "z=0*.out"
files = ask("filenames to plot",files)

filelist = findfile(files)
nfiles = n_elements(filelist)

if (nfiles eq 1) and (strlen(filelist(0)) eq 0) then begin
  print, "No files found"
  stop
endif

psfile = "idl.ps"
psfile = ask("ps file name ('' for screen)",psfile)

ppp = tostr(nfiles)
ppp = fix(ask("plots per page",ppp))

;titles = ['Speed of Light = 1.0 * C','Speed of Light = 0.005 * C']

for ifile=0,nfiles-1 do begin

  gettype,filelist(ifile),filetype,npict

  openfile,1,filelist(ifile),filetype

  gethead,1,filetype,phys, $
      headline,it,filetime,gencoord,ndim,neqpar,nw,nx,eqpar,variables

  close,1

  openfile,1,filelist(ifile),filetype

  getpict,1,filetype,npict(0),x,w,headline,it,filetime

  close,1

  readtransform,ndim,nx,gencoord,transform,nxreg,wregpad,$
		physics,nvector,vectors,grid,doask, rid_grid

  regulargrid,x_old,nxreg_old,x,xreg,nxreg,dxreg,$
                     w,wreg,nw,wregpad,triangles,rid_grid

  if (ifile eq 0) then begin
    if strpos(variables(3),'mx') gt -1 then begin
      data = fltarr(nfiles, nxreg(0), nxreg(1), nw+4)
    endif else begin
      data = fltarr(nfiles, nxreg(0), nxreg(1), nw)
    endelse
    time = fltarr(nfiles)
  endif

  time(ifile) = filetime
  data(ifile,*,*,0:nw-1) = wreg
  data(ifile,*,*,0) = alog10(data(ifile,*,*,0))

  if strpos(variables(3),'mx') gt -1 then begin
    data(ifile,*,*,nw) = 50.0*wreg(*,*,1)/wreg(*,*,0)
    data(ifile,*,*,nw+1) = 50.0*wreg(*,*,2)/wreg(*,*,0)
    data(ifile,*,*,nw+2) = 50.0*wreg(*,*,3)/wreg(*,*,0)
    data(ifile,*,*,nw+3) = sqrt(data(ifile,*,*,nw)^2   + $
                                data(ifile,*,*,nw+1)^2 + $
                                data(ifile,*,*,nw+2)^2)
    variables_new = [variables(0:nw+1),'Vx','Vy','Vz','|V|']
    nw = nw + 4
    variables = variables_new
  endif

endfor

rBody = rid_grid(5)

for i=0,nw-1 do print, i,' ',variables(i+2)
if n_elements(plotvar) eq 0 then plotvar = 0
plotvar = fix(ask('plotvar',tostr(plotvar)))

if (strlen(psfile) gt 0) then begin
  if ppp gt 1 then setdevice, psfile, 'p',4,0.95    $
  else setdevice, psfile, 'l',4,0.95
endif

maxi = max(data(*,*,*,plotvar))
mini = min(data(*,*,*,plotvar))

if (mini lt 0.0) then begin
  readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"
  maxi = max([abs(maxi),abs(mini)])
  mini = -maxi
endif else begin
readct,ncolors, getenv("IDL_EXTRAS")+"white_red.ct"
;  readct,ncolors, getenv("IDL_EXTRAS")+"jim6.ct"
endelse

clevels = (ncolors)*findgen(30)/29.0 + 1

mini = float(ask('minimum value',string(mini)))
maxi = float(ask('maximum value',string(maxi)))

levels = findgen(30)/29.0*(maxi-mini) + mini
ctlevels = mm(levels)

space = 0.02
pos_space, ppp, space, sizes, ny=ppp

aspect = (rid_grid(2) - rid_grid(0))/(rid_grid(3) - rid_grid(1))

for ifile = 0, nfiles-1 do begin

  pn = ifile mod ppp

  if pn eq 0 then begin
    plotdumb
    polyfill, [0,0,1,1,0],[0,1,1,0,0], color = ncolors
  endif

  if pn eq ppp-1 then begin
    xtn = strarr(20) 
    xt = variables(0)+' (Re)'
  endif else begin
    xtn = strarr(20)+' '
    xt = ' '
  endelse

  get_position, ppp, space, sizes, pn, pos
  pos(2) = pos(0) + aspect*2*(pos(2)-pos(0))/2.0

  diff = 0.5 - mean(pos([0,2]))
  pos([0,2]) = pos([0,2]) + diff

  if (pos(2) gt 0.97) then begin
    xr = pos(2) - pos(0)
    if (xr gt 0.97) then begin
      ppp = ppp*2
      pos_space, ppp, space, sizes, nx = 2
      get_position, ppp, space, sizes, pn, pos
      pos(2) = pos(0) + aspect*2*(pos(2)-pos(0))/2.0
      diff = 0.5 - mean(pos([0,2]))
      pos([0,2]) = pos([0,2]) + diff
    endif else begin
      pos(2) = 0.97
      pos(0) = pos(2) - xr
    endelse
  endif

  data_to_plot = reform(data(ifile,*,*,plotvar))

  loc = where(data_to_plot lt mini,count)
  if count gt 0 then data_to_plot(loc) = mini*1.001

  loc = where(data_to_plot gt maxi,count)
  if count gt 0 then data_to_plot(loc) = maxi*0.99

  contour, data_to_plot, xreg(*,*,0), xreg(*,*,1), 		$
    levels = levels, c_colors = clevels,  /fill, pos = pos,		$
    color = 0, /noerase, xtitle = xt, ytitle = variables(1)+' (Re)',	$
    ystyle = 1, xstyle = 1, xtickname = xtn

  if (max(xreg(*,*,0)) gt 0) and (min(xreg(*,*,0)) lt 0) then begin
    theta = findgen(37)/36.0*2.0*!pi
    polyfill, rBody*cos(theta), rBody*sin(theta), color = 0
  endif

;  xyouts, pos(0)+0.02,pos(1)+0.02,titles(ifile),/norm,color = 200

  if (plotgrid) then oplot, x(*,0,0),x(*,0,1),psym=3, color = 0

  if (strpos(variables(1),'z') gt -1) then begin

    b1 = reform(data(ifile,*,*,4))
    b2 = reform(data(ifile,*,*,6))
    btotal = sqrt(b1^2 + b2^2 + 0.00001)
    v1 = b1/btotal
    v2 = b2/btotal

    velvector = 60
    velpos = fltarr(velvector,2)
    velpos(0:59,0) = findgen(velvector)* n_elements(xreg(*,0)) / velvector
    velpos(0:59,1) = n_elements(xreg(0,*))*0.5

;    velpos(30:59,0) = 2*findgen(velvector/2)* n_elements(xreg(*,0)) / velvector
;    velpos(30:59,1) = n_elements(xreg(0,*))*0.25

;    velvector = 100

  endif else begin

    v1 = reform(data(ifile,*,*,1))
    v2 = reform(data(ifile,*,*,2))
    vtotal = sqrt(v1^2 + v2^2 + 0.00001)
    v1 = v1/vtotal
    v2 = v2/vtotal

    velvector = 60
;    velpos = fltarr(velvector,2)
;    velpos(0:29,1) = 2*findgen(velvector/2)* n_elements(xreg(0,*)) / velvector
;    velpos(0:29,0) = n_elements(xreg(*,0))*0.9

;    velpos(30:59,1) = 2*findgen(velvector/2)* n_elements(xreg(0,*)) / velvector
;    velpos(30:59,0) = n_elements(xreg(*,0))*0.5

  endelse

  !p.position = pos

;  vel,v1,v2,NVECS=velvector,MAXVAL=1.,$
;      NSTEP=6,LENGTH=0.06,HEAD=0.0,$
;      DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE

  vel,v1,v2,NVECS=velvector,MAXVAL=1.,$
      NSTEP=300,LENGTH=1.0,HEAD=0.0,$
      DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE

  v1=-v1 & v2=-v2
  vel,v1,v2,NVECS=velvector,MAXVAL=1.,$
      NSTEP=300,LENGTH=1.,HEAD=0.,$
      DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE

  if (time(ifile) gt 0) then begin
    xyouts, pos(0)+0.02,pos(1)+0.02,tostr(time(ifile))+' s',/norm,color = 0
  endif else begin
    xyouts, pos(0)+0.02,pos(1)+0.02,tostr(time(ifile))+' s',/norm,color = 0
  endelse

  ctpos = [pos(2)+0.01,pos(1), pos(2)+0.03,pos(3)]
  plotct, ncolors, ctpos, ctlevels, variables(plotvar+2), /right, 	$
	color = 0

endfor

end

