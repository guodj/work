
;Values = fltarr(nFiles,4,3,nVars,nAlts)

filelist = findfile('./summary_??.save')
nFiles = n_elements(filelist)

for iFile = 0,nFiles-1 do begin

    print, 'Restoring file : ',filelist(iFile)
    restore, filelist(iFile)

    if (iFile eq 0) then begin

        nAlts = n_elements(values(0,0,0,0,*))
        nVars = n_elements(values(0,0,0,*,0))
        nTimes = n_elements(values(*,0,0,0,0))
        all = fltarr(nFiles, nTimes, 4, 3, nVars, nAlts)
        iMean_ = 0
        iStd_  = 1
        iMin_  = 2
        iMax_  = 3
        iGlobal_ = 0
        iSouth_  = 1
        iNorth_  = 2 

    endif

    all(iFile, *,*,*,*,*) = values

endfor

imsmm = fix(ask('mean, std, min, max (0,1,2,3)','0'))
igsn = fix(ask('global, south, north (0,1,2)','0'))
iNorm = fix(ask('base results (0) or normalized results (1)','0'))

smsmm = ['mean','std','min','max']
sgsn = ['global','south','north']
sNorm = ['base','norm']

vars(33) = 'eden'
vars(18) = 'vertv'
display,vars
if (n_elements(iVar) eq 0) then iVar = 9
iVar = fix(ask('variable to plot',tostr(iVar)))

alts = reform(values(0,0,0,2,*))/1000.0
;display, string(alts,format='(f5.1)')
;if (n_elements(iAlt) eq 0) then iAlt = 30
;iAlt = fix(ask('altitude top plot',tostr(iAlt)))
iAlt = [14,24,33,40,47]

time = findgen(nTimes)*24.0/(nTimes-1)

psfile = vars(iVar)+'_'+sNorm(iNorm)+'_'+smsmm(imsmm)+'_'+sgsn(igsn)+'.ps'
setdevice, psfile, 'p', 5
plotdumb
makect,'bry'

if (iNorm eq 1) then begin
    for iFile = 0, nFiles-1 do begin
        all(iFile, *, imsmm,igsn,iVar, *) = $
          (all(iFile, *,imsmm,igsn,iVar, *) - $
           all(nFiles-1, *,imsmm,igsn,iVar,*)) / $
          all(nFiles-1, *,imsmm,igsn,iVar,*) * 100.0
    endfor
endif

ppp = n_elements(iAlt)
space = 0.02
pos_space, ppp, space, sizes, ny = ppp

for i = 0, n_elements(iAlt)-1 do begin

    get_position, ppp, space, sizes, i, pos, /rect

    pos(0) = pos(0) + 0.1

    mini = min(all(*, *,imsmm,igsn,iVar,iAlt(i)))
    maxi = max(all(*, *,imsmm,igsn,iVar,iAlt(i)))
    r = maxi-mini
    yrange = [mini-r*0.1,maxi+r*0.1]

    if i eq n_elements(iAlt)-1 then xtn = strarr(10) $
    else xtn = strarr(10)+' '

    plot, time, all(nFiles-1, *,imsmm,igsn,iVar,iAlt(i)), /noerase, $
      yrange = yrange, ystyle = 1, xstyle = 1, $
      thick = 3, pos = pos, ytitle = vars(iVar), $
      xtickname = xtn, xticks = 9, xtickv = findgen(9)*24/8, xminor = 3

    for iFile = 0,nFiles-2 do begin
        oplot, time,all(iFile, *, imsmm,igsn,iVar,iAlt(i)), $
          color = 255-(iFile+1)*250.0/nFiles
    endfor
    iFile = nFiles-1
    oplot, time,all(iFile, *, imsmm,igsn,iVar,iAlt(i)), $
      color = 0, thick = 3

    xyouts, time(nTimes-2), maxi, $
      string(alts(iAlt(i)), format='(f5.1)')+' km', $
      alignment = 1

    if (igsn eq iSouth_) then begin
        oplot, [6 , 6], yrange, linestyle = 0
        xyouts, 6.1, mini-0.08*r, 'Noon'
        oplot, [18,18], yrange, linestyle = 1
        xyouts, 18.1, mini-0.08*r, 'Midnight'
        if (i eq 0) then begin
            xyouts, 0, maxi+r*0.12, 'Southern Hemisphere'
        endif
    endif

    if (igsn eq iNorth_) then begin
        oplot, [18,18], yrange, linestyle = 0
        xyouts, 18.1, mini-0.08*r, 'Noon'
        oplot, [ 6, 6], yrange, linestyle = 1
        xyouts, 6.1, mini-0.08*r, 'Midnight'
        if (i eq 0) then begin
            xyouts, 0, maxi+r*0.12, 'Northern Hemisphere'
        endif
    endif

    if (i eq 0) then begin
        if (iNorm eq 1) then xyouts, 24, maxi+r*0.12, 'Normalized',align=1
        xyouts, 12, maxi+r*0.12, smsmm(imsmm), align=0.5
    endif
endfor

closedevice

end
