
nAE = 3
MaxAE = 1200.0
dAE = MaxAE/nAE

nMLTs = 8
nMAGs = 4
MaxMAG = 6.0
dMAG = MaxMAG/nMAGs

HallBinned     = fltarr(nMLTs, nAE)
CCBinned       = fltarr(nMLTs, nAE)
CCBinned_amie  = fltarr(nMLTs, nAE)
CCLatBinned    = fltarr(nMLTs, nAE)
LatDiffBinned  = fltarr(nMLTs, nAE)
absCCLatBinned = fltarr(nMLTs, nAE)
nMagBinned     = fltarr(nMLTs, nAE)
nPts           = fltarr(nMLTs, nAE)

HallBinned_mag     = fltarr(nMAGs, nAE)
CCBinned_mag       = fltarr(nMAGs, nAE)
CCBinned_amie_mag  = fltarr(nMAGs, nAE)
CCLatBinned_mag    = fltarr(nMAGs, nAE)
LatDiffBinned_mag  = fltarr(nMAGs, nAE)
absCCLatBinned_mag = fltarr(nMAGs, nAE)
nMagBinned_mag     = fltarr(nMAGs, nAE)
nPts_mag           = fltarr(nMAGs, nAE)


psfile = ask('ps file','metrics_hall.ps')

setdevice,psfile,'p',4,0.9

threshold = float(ask('threshold for cc w/AVE cut off','0.5'))
nMagThreshold = float(ask('threshold for nMags','0'))

filelist = findfile('*_metrics')
nFiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

  amie_file = filelist(iFile)

  print, amie_file

  openr, 2, amie_file

  readf,2,nOrbits
  readf,2,nMetrics

  Metrics = fltarr(nMetrics,2,nOrbits)
  MetricsTimes = dblarr(nOrbits)
  MetricTitle = strarr(nMetrics)
  line = ''
  for i=0,nMetrics-1 do begin
    readf,2,line
    MetricTitle(i) = line
  endfor
  readf,2,metrics
  readf,2,MetricsTimes

  close,2

  for i=0,nOrbits-1 do for j=0,1 do begin
    mini = min(metrics(*,j,i))
    if metrics(3,j,i) ge threshold and $
       metrics(5,j,i) ge nMagThreshold and $
       mini gt -999.0 and $
       metrics(2,j,i) gt 0.0 then begin

      mlt = metrics(7,j,i)/24.0 * nMLTs
      if mlt gt nMLTs-1 then mlt = 0

      ae = min([fix(metrics(6,j,i)/dAE),nAE-1])

      if (metrics(5,j,i) eq 0) then mag = 0
      if (metrics(5,j,i) eq 1) then mag = 1
      if (metrics(5,j,i) eq 2) then mag = 1
      if (metrics(5,j,i) eq 3) then mag = 2
      if (metrics(5,j,i) eq 4) then mag = 2
      if (metrics(5,j,i) ge 5) then mag = 3

;      mag = min([fix(metrics(5,j,i)/dMAG),nMAGs-1])

      nPts(mlt,ae) = nPts(mlt,ae) + 1.0

      HallBinned(mlt,ae)  = HallBinned(mlt,ae)  + $
	metrics(8,j,i)*metrics(3,j,i)
      CCBinned_amie(mlt,ae) = CCBinned_amie(mlt,ae)    + $
        metrics(10,j,i)*metrics(3,j,i)
      CCBinned(mlt,ae)    = CCBinned(mlt,ae)    + metrics(3,j,i)
      LatDiffBinned(mlt,ae) = LatDiffBinned(mlt,ae) + $
	metrics(9,j,i)*metrics(3,j,i)
      CCLatBinned(mlt,ae) = CCLatBinned(mlt,ae) + $
	metrics(11,j,i)*metrics(3,j,i)
      absCCLatBinned(mlt,ae) = absCCLatBinned(mlt,ae) + $
	abs(metrics(11,j,i))*metrics(3,j,i)
      nMagBinned(mlt,ae)  = nMagBinned(mlt,ae)  + $
	metrics(5,j,i)*metrics(3,j,i)

;      if (metrics(7,j,i) gt 18.5 or metrics(7,j,i) lt 2.5) then begin
      if (mlt eq 0 or mlt ge nMLTs-2) then begin

        nPts_mag(mag,ae) = nPts_mag(mag,ae) + 1.0

        HallBinned_mag(mag,ae)  = HallBinned_mag(mag,ae)  + $
	  metrics(8,j,i)*metrics(2,j,i)
        CCBinned_amie_mag(mag,ae)  = CCBinned_amie_mag(mag,ae)  + $
	  metrics(10,j,i)*metrics(2,j,i)
        CCBinned_mag(mag,ae)    = CCBinned_mag(mag,ae)    + metrics(2,j,i)
        CCLatBinned_mag(mag,ae) = CCLatBinned_mag(mag,ae) + $
	  metrics(11,j,i)*metrics(2,j,i)
        LatDiffBinned_mag(mag,ae) = LatDiffBinned_mag(mag,ae) + $
	  metrics(9,j,i)*metrics(2,j,i)
        absCCLatBinned_mag(mag,ae) = absCCLatBinned_mag(mag,ae) + $
	  abs(metrics(9,j,i))*metrics(2,j,i)
        nMagBinned_mag(mag,ae)  = nMagBinned_mag(mag,ae)  + $
	  metrics(5,j,i)*metrics(2,j,i)

      endif

    endif
 
  endfor

endfor

loc = where(CCBinned gt 0)

HallBinned(loc)  = HallBinned(loc) / CCBinned(loc)
CCLatBinned(loc) = CCLatBinned(loc) / CCBinned(loc)
CCBinned_amie(loc) = CCBinned_amie(loc) / CCBinned(loc)
LatDiffBinned(loc) = LatDiffBinned(loc) / CCBinned(loc)
absCCLatBinned(loc) = absCCLatBinned(loc) / CCBinned(loc)
nMagBinned(loc)  = nMagBinned(loc) / CCBinned(loc)

CCBinned(loc)    = CCBinned(loc) / nPts(loc)

loc = where(CCBinned_mag gt 0)

HallBinned_mag(loc)  = HallBinned_mag(loc) / CCBinned_mag(loc)
CCBinned_amie_mag(loc) = CCBinned_amie_mag(loc) / CCBinned_mag(loc)
CCLatBinned_mag(loc) = CCLatBinned_mag(loc) / CCBinned_mag(loc)
LatDiffBinned_mag(loc) = LatDiffBinned_mag(loc) / CCBinned_mag(loc)
absCCLatBinned_mag(loc) = absCCLatBinned_mag(loc) / CCBinned_mag(loc)
nMagBinned_mag(loc)  = nMagBinned_mag(loc) / CCBinned_mag(loc)

CCBinned_mag(loc)    = CCBinned_mag(loc) / nPts_mag(loc)

ae  = fltarr(nMLTs, nAE)
mlt = fltarr(nMLTs, nAE)

for i=0,nMLTs-1 do ae(i,*) = findgen(nAE)*dAE
for i=0,nAE-1 do mlt(*,i) = findgen(nMLTs)/nMLTs*24.0

locnpts = where(npts gt 1.0)

plotdumb

readct, ncolors, getenv("IDL_EXTRAS")+"blue_white_red.ct"

ppp = 6
space = 0.01
pos_space, ppp, space, sizes, nx = 1

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

maxi = alog10(fix(float(max(npts))/5.0+1)*5)

levels = findgen(21)*maxi/20.0
clevels = findgen(21)*ncolors/20.0
contour, alog10(npts(locnpts)), mlt(locnpts), ae(locnpts), $
	/follow, levels = levels, c_color = clevels, $
	pos = pos, /noerase, ystyle = 1, xstyle = 1, /fill, $
	ytitle = 'AE (nT)', xtickname=strarr(10)+' ', /irr, xrange=[0,24], yrange = [0,max(ae)]

contour, alog10(npts(locnpts)), mlt(locnpts), ae(locnpts), /follow, levels = levels, $
	pos = pos, /noerase, ystyle = 1, xstyle = 1, $
	xtickname=strarr(10)+' ', /irr, xrange=[0,24], yrange = [0,max(ae)]

oplot, [6,6],[0,2000]
oplot, [12,12],[0,2000]
oplot, [18,18],[0,2000]

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.05,pos(3)]
plotct, ncolors, ctpos, mm(levels), 'alog10(nPts)', /right

;-----------------------------------------------------------------------

get_position, ppp, space, sizes, 1, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

levels = findgen(21)*30.0/20.0
clevels = findgen(21)*ncolors/20.0

if max(nMagBinned) eq 0 then nMagBinned = levels(10)

contour, nMagBinned(locnpts), mlt(locnpts), ae(locnpts), $
        /follow, levels = levels, c_color = clevels, $
	pos = pos, /noerase, ystyle = 1, xstyle = 1, /fill, $
	ytitle = 'AE (nT)', xtickname=strarr(10)+' ',/irr, xrange=[0,24], yrange = [0,max(ae)]

contour, nMagBinned(locnpts), mlt(locnpts), ae(locnpts), $
        /follow, levels = levels, $
	pos = pos, /noerase, ystyle = 5, xstyle = 5,/irr, xrange=[0,24], yrange = [0,max(ae)]

oplot, [6,6],[0,2000]
oplot, [12,12],[0,2000]
oplot, [18,18],[0,2000]

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.05,pos(3)]
plotct, ncolors, ctpos, mm(levels), 'nMag', /right

;-----------------------------------------------------------------------

get_position, ppp, space, sizes, 2, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

levels = findgen(21)*1.0/20.0
clevels = findgen(21)*ncolors/20.0
contour, CCBinned_amie(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, c_color = clevels, $
	pos = pos, /noerase, ystyle = 1, xstyle = 1, /fill, $
	ytitle = 'AE (nT)', xtickname=strarr(10)+' ',/irr, xrange=[0,24], yrange = [0,max(ae)]

contour, CCBinned_amie(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, $
	pos = pos, /noerase, ystyle = 5, xstyle = 5,/irr, xrange=[0,24], yrange = [0,max(ae)]

oplot, [6,6],[0,2000]
oplot, [12,12],[0,2000]
oplot, [18,18],[0,2000]

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.05,pos(3)]
plotct, ncolors, ctpos, mm(levels), 'Correlation', /right

;-----------------------------------------------------------------------

get_position, ppp, space, sizes, 3, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

levels = findgen(21)*max(abs(CCLatBinned))/10.0 - max(abs(CCLatBinned))
clevels = findgen(21)*ncolors/20.0
contour, CCLatBinned(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, c_color = clevels, $
	pos = pos, /noerase, ystyle = 1, xstyle = 1, /fill, $
	ytitle = 'AE (nT)', xtickname=strarr(10)+' ',/irr, xrange=[0,24], yrange = [0,max(ae)] 

contour, CCLatBinned(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, $
	pos = pos, /noerase, ystyle = 5, xstyle = 5,/irr, xrange=[0,24], yrange = [0,max(ae)] 

oplot, [6,6],[0,2000]
oplot, [12,12],[0,2000]
oplot, [18,18],[0,2000]

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.05,pos(3)]
plotct, ncolors, ctpos, mm(levels), 'Lat. Diff.', /right

;-----------------------------------------------------------------------

get_position, ppp, space, sizes, 4, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

levels = findgen(21)*max(absCCLatBinned)/20.0
clevels = findgen(21)*ncolors/20.0
contour, absCCLatBinned(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, $
	c_color = clevels, $
	pos = pos, /noerase, ystyle = 1, xstyle = 1, /fill, $
	ytitle = 'AE (nT)', xtickname=strarr(10)+' ',/irr, xrange=[0,24], yrange = [0,max(ae)]

contour, absCCLatBinned(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, $
	pos = pos, /noerase, ystyle = 5, xstyle = 5,/irr, xrange=[0,24], yrange = [0,max(ae)]

oplot, [6,6],[0,2000]
oplot, [12,12],[0,2000]
oplot, [18,18],[0,2000]

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.05,pos(3)]
plotct, ncolors, ctpos, mm(levels), 'abs(Lat. Diff.)', /right

;-----------------------------------------------------------------------

get_position, ppp, space, sizes, 5, pos, /rect
pos(0) = pos(0) + 0.05
pos(2) = pos(2) - 0.05

mini = min(HallBinned(locnpts))
maxi = max([max(HallBinned(locnpts)), abs(mini)])
levels = findgen(21)*2.0*maxi/20.0 - maxi
clevels = findgen(21)*ncolors/20.0
contour, HallBinned(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, c_color = clevels, $
	pos = pos, /noerase, ystyle = 1, xstyle = 1, /fill, $
	ytitle = 'AE (nT)', xtitle = 'MLT',/irr, xrange=[0,24], yrange = [0,max(ae)]

contour, HallBinned(locnpts), mlt(locnpts), ae(locnpts), /follow, levels = levels, $
	pos = pos, /noerase, ystyle = 5, xstyle = 5,/irr, xrange=[0,24], yrange = [0,max(ae)]

oplot, [6,6],[0,2000]
oplot, [12,12],[0,2000]
oplot, [18,18],[0,2000]

ctpos = [pos(2)+0.01,pos(1),pos(2)+0.05,pos(3)]
plotct, ncolors, ctpos, mm(levels), '% Diff Hall', /right

closedevice

end
