
nLts = 49
nLats = 91

dLts = 360.0/(nLts-1)
dLat = 180.0/(nLats-1)

tec = fltarr(nLts,nLats,3)
nTec = intarr(nLts,nLats)

lats = findgen(nLats)*dLat - 90.0
lts  = findgen(nLts)*dLts

lats2d = fltarr(nLts, nLats)
lts2d  = fltarr(nLts, nLats)

for i = 0, nLts-1 do lats2d(i,*) = lats
for i = 0, nLats-1 do lts2d(*,i) = lts

filelist = findfile('tec_gitm*.txt')
nFiles = n_elements(filelist)

itime = intarr(6)
nPts = 0
tmp = fltarr(5)

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   close,1
   openr,1,file

   readf,1,itime
   ut = (float(itime(3))+float(itime(4))/60.0+float(itime(5))/3600.0)*15.0

   print, ut

   readf,1,nPts

   for iPt = 0, nPts-1 do begin

      readf,1,tmp

      iLat = round((tmp(0)+90.0)/dlat)
      iLt = round((tmp(1)+360.0*2+ut)/dLts) mod nLts
      tec( iLt,iLat, *) = tec(iLt, iLat, *) + tmp(2:4)
      nTec(iLt,iLat) = nTec(iLt,iLat) + 1

   endfor

   close,1

endfor

print,'Smoothing...'

l = where(nTec gt 1,c)
if (c gt 0) then begin
   for i=0,2 do begin
      vs = reform(tec(*,*,i))
      vs(l) = vs(l)/nTec(l)
      smooth_image, vs, 0.0
      tec(*,*,i) = vs
      tec(0,*,i) = (tec(0,*,i) + tec(nLts-1,*,i))/2.0
      tec(nLts-1,*,i) = tec(0,*,i)
   endfor
endif

models = ['GPS','GITM','TIEGCM']

c_a_to_ymd,itime,ymd

for iModel = 0,2 do begin

   model = models(iModel)
   setdevice, 'tec_'+ymd+'_'+model+'.ps','p',5

   makect,'wyrb'

   ppp = 4
   space = 0.01
   pos_space, ppp, space, sizes

   plotdumb

   get_position, ppp, space, sizes, 0, pos

   no00 = 1
   no06 = 1
   no12 = 0
   no18 = 0

nLevels = 31

mini = 0.0
maxi = float(fix(max(tec))+1.0)
maxi = fix(maxi/10.0+1)*10.0

maxrange = 30.0
contour_circle, reform(tec(*,*,iModel)), lts, lats, $
                mini = mini, maxi = maxi/2, $
                nLevels = nLevels, $
                no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                pos = pos, $
                maxrange = maxrange

   xc = (pos(2)+pos(0))/2
   xr = (pos(2)-pos(0))/2 * 1.01
   yc = (pos(3)+pos(1))/2
   yr = (pos(3)-pos(1))/2 * 1.01
   xp = xc - xr*sin(!pi/4)
   yp = yc + yr*sin(!pi/4)
   xyouts, xp, yp, 'North', $
           /norm, charsize = 0.9, align = 0.5, orient = 45

   xm = 0.5-xr/2
   xp = 0.5+xr/2
   ctpos = [xm, pos(3)+0.01, xp, pos(3)+0.03]
   range = [mini,maxi/2]
   units = 'TEC (10!U16!N/m!U2!N)'
   ncolors = 255
   plotct, ncolors, ctpos, range, units, /bottom

get_position, ppp, space, sizes, 1, pos

no00 = 1
no06 = 0
no12 = 0
no18 = 1

contour_circle, reform(tec(*,*,iModel)), lts, -lats, $
                mini = mini, maxi = maxi/2, $
                nLevels = nLevels, $
                no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                pos = pos, $
                maxrange = maxrange

xc = (pos(2)+pos(0))/2
xr = (pos(2)-pos(0))/2 * 1.01
yc = (pos(3)+pos(1))/2
yr = (pos(3)-pos(1))/2 * 1.01
xp = xc + xr*sin(!pi/4)
yp = yc + yr*sin(!pi/4)
xyouts, xp, yp, 'South', $
        /norm, charsize = 0.9, align = 0.5, orient = -45

get_position, ppp, space, sizes, 2, pos1
get_position, ppp, space, sizes, 3, pos2

pos = pos1
pos(2) = pos2(2)-0.04

nLevels = 31

!p.position = pos

map_set, 0.0, 180.0, /noerase

levels = findgen(nlevels)*(maxi-mini)/(nlevels-1) + mini

tec1 = reform(tec(*,*,iModel))

l = where(tec1 lt levels(1),c)
if (c gt 0) then tec1(l) = levels(1)

contour, tec1, lts2d, lats2d, $
         /follow, /cell_fill, /over, $
         levels = levels

nLevels = 11
levels = findgen(nlevels)*(maxi-mini)/(nlevels-1) + mini
contour, tec1, lts2d, lats2d, $
         /follow, /over, $
         levels = levels, thick = 2

map_grid, lats = findgen(19)*10-90, glinethick=1

contour, lats2d, lts2d, lats2d, $
         /over, thick = 5, c_linestyle = [2,2,2,2], $
         levels = [-60,-30,30,60]

dy = (pos(3)-pos(1))/6

x = pos(0)-0.01
y = pos(3)-dy/2
xyouts, x, y, 'N.P.', orient = 90, align = 0.5, /norm
y = pos(3)-dy*1.5
xyouts, x, y, 'N.M.L.', orient = 90, align = 0.5, /norm
y = pos(3)-dy*3
xyouts, x, y, 'EQ.', orient = 90, align = 0.5, /norm
y = pos(1)+dy*1.5
xyouts, x, y, 'S.M.L.', orient = 90, align = 0.5, /norm
y = pos(1)+dy*0.5
xyouts, x, y, 'S.P.', orient = 90, align = 0.5, /norm
dx = (pos(2)-pos(0))/8.0
for i=0,8 do begin
   x = pos(0) + dx*i
   y = pos(1)-0.02
   xyouts, x, y, tostr(i*3,2)+' LT', /norm, align = 0.5, charsize=0.9
endfor

maxs = 'Max : '+string(max(tec1),format="(f5.1)")+'(10!U16!N/m!U2!N)'

xp = pos(2)
yp = pos(1)-yr/10.0*2
xyouts, xp, yp, maxs, charsize = 0.9, align = 1.0, /norm

xp = pos(0)
yp = pos(3)+(pos(3)-pos(1))*1.1 ; pos(3)+yr/20.0
c_a_to_s, itime, stime
xyouts, xp, yp, strmid(stime,0,9), $
        /norm, charsize = 1.1, align = 0.0

xp = pos(2)
xyouts, xp, yp, model, /norm, charsize=1.1, align=1.0

  
ctpos = [pos(2)+0.005, pos(1), pos(2)+0.02, pos(3)]
range = [mini,maxi]
units = 'TEC (10!U16!N/m!U2!N)'
ncolors = 255
plotct, ncolors, ctpos, range, units, /right

closedevice

endfor

end

