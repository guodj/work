
filelist = findfile('3D*0000.bin')
nFiles = n_elements(filelist)

;nFiles = 6

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)
   print, 'Reading File : ',filename
   read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                           vars, data, rb, cb, bl_cnt, iTime, Version

   alt = reform(data(2,0,0,*)) / 1000.0
   lat = reform(data(1,*,*,*)) / !dtor
   lon = reform(data(0,2:nLons-3,0,0)) / !dtor

   t = reform(data(15,*,*,*))

   if (iFile eq 0) then begin
      tSave = fltarr(nFiles,4,nLons-4)
      l = where(alt gt 240.0)
      iAlt = l(0)
      iLon = 3
      iLat = 1
   endif

   tSave(iFile,0,*)=shift(reform(t(2:nLons-3,(nLats-4)/2 + 1 - iLat*5,iAlt)),iLon*iFile)
   tSave(iFile,1,*)=shift(reform(t(2:nLons-3,(nLats-4)/2 + 1,iAlt)),iLon*iFile)
   tSave(iFile,2,*)=shift(reform(t(2:nLons-3,(nLats-4)/2 + 2,iAlt)),iLon*iFile)
   tSave(iFile,3,*)=shift(reform(t(2:nLons-3,(nLats-4)/2 + 2 + iLat*5,iAlt)),iLon*iFile)

endfor

setdevice, 'mtm.ps', 'p', 5

mini = min(tSave)
maxi = max(tSave)
r = (maxi-mini)
mini = mini - 0.1*r
maxi = maxi + 0.1*r

ppp = 3
space = 0.05
pos_space, ppp, space, sizes, ny = ppp

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0)+0.1
plot, lon/15,tsave(0,0,*), yrange = [mini,maxi], ystyle = 1, $
      /noerase, xstyle = 1, pos = pos, $
      ytitle = 'Temperature (K, -5 Lat)', $
      xtickv = findgen(9)*3, xticks = 9, xminor = 3
for iFile = 1, nFiles-1 do oplot, lon/15,tsave(iFile,0,*)

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0)+0.1
plot, lon/15,(tsave(0,1,*)+tsave(0,2,*))/2, yrange = [mini,maxi], ystyle = 1, $
      /noerase, xstyle = 1, pos = pos, $
      ytitle = 'Temperature (K, Equator)', $
      xtickv = findgen(9)*3, xticks = 9, xminor = 3
for iFile = 1, nFiles-1 do oplot, lon/15,(tsave(iFile,1,*)+tsave(iFile,2,*))/2

get_position, ppp, space, sizes, 2, pos, /rect
pos(0) = pos(0)+0.1
plot, lon/15,tsave(0,3,*), yrange = [mini,maxi], ystyle = 1, $
      /noerase, xstyle = 1, pos = pos, $
      ytitle = 'Temperature (K, +5 Lat)', $
      xtickv = findgen(9)*3, xticks = 9, xminor = 3, xtitle = 'Local Time (hr)'
for iFile = 1, nFiles-1 do oplot, lon/15, tsave(iFile,3,*)

closedevice

end

