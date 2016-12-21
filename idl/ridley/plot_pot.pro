list = findfile('ionosphere_n00000*.dat.save')
nfiles = n_elements(list)

cpcp = fltarr(2,nfiles)
ave_hall = fltarr(2,nfiles)
ave_ped = fltarr(2,nfiles)
min_fac = fltarr(2,nfiles)
max_fac = fltarr(2,nfiles)
iter = lonarr(nfiles)

for l=0,nfiles-1 do begin

  file = list(l)

  iter(l) = long(strmid(file,strpos(file,'_n')+2,6))
  print, 'Reading from file : ',file
  restore,file

  if l eq 0 then begin
    all_pot = fltarr(nfiles,2,nlons,nlats)
    all_hall = fltarr(nfiles,2,nlons,nlats)
    all_ped = fltarr(nfiles,2,nlons,nlats)
    all_fac = fltarr(nfiles,2,nlons,nlats)
  endif

  all_pot(l,0,*,*) = data(0,8,*,*)
  all_hall(l,0,*,*) = data(0,5,*,*)
  all_ped(l,0,*,*) = data(0,6,*,*)
  all_fac(l,0,*,*) = data(0,7,*,*)

  all_pot(l,1,*,*) = data(1,8,*,*)
  all_hall(l,1,*,*) = data(1,5,*,*)
  all_ped(l,1,*,*) = data(1,6,*,*)
  all_fac(l,1,*,*) = data(1,7,*,*)

  cpcp(0,l) = max(all_pot(l,0,*,*)) - min(all_pot(l,0,*,*))
  ave_hall(0,l) = mean(all_hall(l,0,*,*))
  ave_ped(0,l) = mean(all_ped(l,0,*,*))
  min_fac(0,l) = min(all_fac(l,0,*,*))
  max_fac(0,l) = max(all_fac(l,0,*,*))

  cpcp(1,l) = max(all_pot(l,1,*,*)) - min(all_pot(l,1,*,*))
  ave_hall(1,l) = mean(all_hall(l,1,*,*))
  ave_ped(1,l) = mean(all_ped(l,1,*,*))
  min_fac(1,l) = min(all_fac(l,1,*,*))
  max_fac(1,l) = max(all_fac(l,1,*,*))

endfor

setdevice

fac = reform(all_fac(nfiles-1,0,*,*))
fac = reform(all_pot(nfiles-1,0,*,*))
mr = 40.0
loc = where(reform(data(0,3,0,*)) le mr)
rang = reform(data(0,3,*,loc))
lons = reform(data(0,4,*,loc))*!pi/180.0 + !pi/2
xpos = rang*cos(lons)
ypos = rang*sin(lons)

ppp = 1
space = 0.01
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos

readct,ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"

maxi = max(abs(fac(*,loc)))*1.05
mini = -maxi
levels = (maxi-mini)*findgen(9)/8.0 + mini
levels = [mini*(5.0-findgen(5))/5.0,maxi*(findgen(5)+1)/5.0]
c_levels = (maxi-mini)*findgen(30)/29.0 + mini
c_colors = (ncolors-1)*findgen(30)/29.0 + 1

get_position, ppp, space, sizes, 0, pos
contour, fac(*,loc), xpos, ypos, /follow, nlevels=30, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,	$
	/cell_fill, /noerase
contour, fac(*,loc), xpos, ypos, /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase, c_linestyle = 3.0*(levels lt 0.0)
plotmlt, mr, /no00
mini = min(fac(*,loc))
maxi = max(fac(*,loc))
maxs = "Max:"+string(maxi,format="(e9.2)")
mins = "Min:"+string(mini,format="(e9.2)")
xyouts, pos(0),pos(1)-0.02, mins, /norm
xyouts, pos(2),pos(1)-0.02, maxs, /norm, align=1.0
xyouts, (pos(0)+pos(2))/2.0,pos(3)+space/2.0,"Field Aligned Current", $
	align=0.5, /norm, charsize = 1.25

plot, iter,cpcp(0,*)
plot, iter,max_fac(0,*), yrange = mm([max_fac,min_fac])
oplot, iter, min_fac(0,*)
oplot, [0.0,max(iter)*10.0], [0.0,0.0], linestyle = 1
oplot, iter, max_fac(0,*)+min_fac(0,*),linestyle = 2

closedevice

end    