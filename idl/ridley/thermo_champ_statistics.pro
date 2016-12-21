
filelist = findfile('*/data/gitm_champ*.save')
;filelist = findfile('./gitm_champ*.save')
nFiles = n_elements(filelist)

for iFile = 0L, nFiles-1 do begin

   file = filelist(iFile)
   print, 'Restoring ',file
   restore, file

   nTimes = n_elements(GitmTime)
   dl = ((GitmLocalTime(1:nTimes-1) - GitmLocalTime(0:nTimes-2)) * $
         (GitmLats(1:nTimes-1) - GitmLats(0:nTimes-2)))

   for i=2L,nTimes-4 do begin
      if (dl(i) lt -0.2 and $
          abs(dl(i-1)) lt 0.05 and $
          dl(i+1) gt 0.2 and $
          abs(dl(i+2)) lt 0.05) then begin
         print, 'bad : ',gitmlocaltime(i:i+2)
         gitmlocaltime(i+1) = -9999.0
         gitmlats(i+1) = -9999.0
      endif
      if (dl(i) gt 0.2 and $
          abs(dl(i-1)) lt 0.05 and $
          dl(i+1) lt -0.2 and $
          abs(dl(i+2)) lt 0.05) then begin
         print, 'bad: ',gitmlocaltime(i:i+2)
         gitmlocaltime(i+1) = -9999.0
         gitmlats(i+1) = -9999.0
      endif
   endfor

   l = where(gitmlats gt -100, c)
   if (c gt 0) then begin

      gitmlats(l) = gitmlats(l) / !dtor

      if (iFile eq 0) then begin
         GitmTimeAll      = GitmTime
         GitmLatsAll      = GitmLats
         GitmLonsAll      = GitmLons
         GitmLocalTimeAll = Gitmlocaltime
         GitmIntRhoAll    = GitmIntRho
         ChampIntRhoAll   = ChampIntRho
      endif else begin
         GitmTimeAll      = [GitmTimeAll,GitmTime]
         GitmLatsAll      = [GitmLatsAll,GitmLats]
         GitmLonsAll      = [GitmLonsAll,GitmLons]
         GitmLocalTimeAll = [GitmLocalTimeAll,Gitmlocaltime]
         GitmIntRhoAll    = [GitmIntRhoAll,GitmIntRho]
         ChampIntRhoAll   = [ChampIntRhoAll,ChampIntRho]
      endelse

      maxi = max([gitmintrho,champintrho])
      fac = 1.0
      l2 = where(gitmintrho gt maxi*fac, c2)
      while (c2 lt 0.005*nTimes) do begin
         fac = fac - 0.01
         l2 = where(gitmintrho gt maxi*fac, c2)
      endwhile
      maxi = maxi*fac

      t = mean(gitmtime)
      c_r_to_a, itime, t
      c_a_to_ymd, itime, ymd
      ym = strmid(ymd,0,6)

      psfile = 'champ_gitm_'+ym+'.ps'

      setdevice, psfile, 'p', 5
      plotdumb

      ppp = 6
      space = 0.04
      pos_space, ppp, space, sizes

      ;----------------------------------------------------------------
      ; Get f107
      ;----------------------------------------------------------------

      f107_read, f107time, f107

      l = where(f107time ge min(gitmtime) and $
                f107time lt max(gitmtime), c)

      f107time = f107time(l)
      f107 = f107(l)

      ;----------------------------------------------------------------
      ; All Data
      ;----------------------------------------------------------------

      l = where(gitmlats gt -100, c)

      plotnum = 0
      get_position, ppp, space, sizes, plotnum, pos
      
      plot, gitmintrho(l), champintrho(l), psym = 3, $
            xrange = [0,maxi], yrange=[0,maxi], $
            xstyle = 1, ystyle = 1, $
            ytitle = 'CHAMP Rho (1e-12 kg/m3)', pos = pos
      oplot, [0,maxi], [0,maxi], thick = 3.0

      r = c_correlate(gitmintrho(l), champintrho(l), 0.0)
      rms = sqrt(mean((gitmintrho(l)-champintrho(l))^2))
      rmsn = sqrt(mean(((gitmintrho(l)-champintrho(l))/champintrho(l))^2))

      xyouts, maxi*0.05, maxi*0.90, '(A) All Data', $
              charsize = 0.7
      xyouts, maxi*0.05, maxi*0.82, 'Corr: '+string(r,format='(f4.2)'), $
              charsize = 0.7
      xyouts, maxi*0.3, maxi*0.12, 'RMS: '+string(rms,format='(f6.3)')+$
              'x10!U-12!N kg/m!U3!N', charsize = 0.7
      xyouts, maxi*0.3, maxi*0.05, 'RMSN: '+string(r,format='(f4.2)'), $
              charsize = 0.7

      ;----------------------------------------------------------------
      ; Drivers
      ;----------------------------------------------------------------

      plotnum = 1
      get_position, ppp, space, sizes, plotnum, pos
      dx = pos(2)-pos(0)
      dy = pos(3)-pos(1)
      pos(0) = pos(0) + dx*0.12
      pos(1) = pos(1) + dy*0.12

      stime = min(gitmtime)
      etime = max(gitmtime)
      time_axis, stime, etime, btr, etr, xtickname, $
                 xtitle, xtickv, xminor, xtickn

      minf = min(f107)
      maxf = max(f107)
      r = maxf-minf
      minf = minf-r*0.1
      maxf = maxf+r*0.1

      plot, f107time-stime, f107, yrange = [minf,maxf], $
            pos = pos, /noerase, thick = 4, $
            xtickname = xtickname, $
            xtitle = xtitle, xminor = xminor, $
            xstyle = 1, ystyle = 1, xrange = [btr,etr], $
            xtickv = xtickv, xticks = xtickn, $
            ytitle = 'F!D10.7!N (SFU)'

      xyouts, btr+(etr-btr)*0.05, maxf-r*0.1, '(B) F10.7', charsize=0.7

      ;----------------------------------------------------------------
      ; North
      ;----------------------------------------------------------------

      l = where(gitmlats gt 60, c)
      
      plotnum = 2
      get_position, ppp, space, sizes, plotnum, pos
      
      plot, gitmintrho(l), champintrho(l), psym = 3, $
            xrange = [0,maxi], yrange=[0,maxi], $
            xstyle = 1, ystyle = 1, $
            ;ytitle = 'CHAMP Rho (1e-12 kg/m3)', $
            pos = pos, /noerase
      oplot, [0,maxi], [0,maxi], thick = 3.0

      r = c_correlate(gitmintrho(l), champintrho(l), 0.0)
      rms = sqrt(mean((gitmintrho(l)-champintrho(l))^2))
      rmsn = sqrt(mean(((gitmintrho(l)-champintrho(l))/champintrho(l))^2))

      xyouts, maxi*0.05, maxi*0.90, '(C) North Polar', $
              charsize = 0.7
      xyouts, maxi*0.05, maxi*0.82, 'Corr: '+string(r,format='(f4.2)'), $
              charsize = 0.7
      xyouts, maxi*0.3, maxi*0.12, 'RMS: '+string(rms,format='(f6.3)')+$
              'x10!U-12!N kg/m!U3!N', charsize = 0.7
      xyouts, maxi*0.3, maxi*0.05, 'RMSN: '+string(r,format='(f4.2)'), $
              charsize = 0.7

      ;----------------------------------------------------------------
      ; South
      ;----------------------------------------------------------------

      l = where(gitmlats gt -100 and gitmlats lt -60, c)
      
      plotnum = 3
      get_position, ppp, space, sizes, plotnum, pos
      
      plot, gitmintrho(l), champintrho(l), psym = 3, $
            xrange = [0,maxi], yrange=[0,maxi], $
            xstyle = 1, ystyle = 1, $
            ;ytitle = 'CHAMP Rho (1e-12 kg/m3)', $
            pos = pos, /noerase
      oplot, [0,maxi], [0,maxi], thick = 3.0

      r = c_correlate(gitmintrho(l), champintrho(l), 0.0)
      rms = sqrt(mean((gitmintrho(l)-champintrho(l))^2))
      rmsn = sqrt(mean(((gitmintrho(l)-champintrho(l))/champintrho(l))^2))

      xyouts, maxi*0.05, maxi*0.90, '(D) South Polar', $
              charsize = 0.7
      xyouts, maxi*0.05, maxi*0.82, 'Corr: '+string(r,format='(f4.2)'), $
              charsize = 0.7
      xyouts, maxi*0.3, maxi*0.12, 'RMS: '+string(rms,format='(f6.3)')+$
              'x10!U-12!N kg/m!U3!N', charsize = 0.7
      xyouts, maxi*0.3, maxi*0.05, 'RMSN: '+string(r,format='(f4.2)'), $
              charsize = 0.7

      ;----------------------------------------------------------------
      ; Low Lats, Ascending Phase
      ;----------------------------------------------------------------

      dl = gitmlats(1:nTimes-1)-gitmlats(0:nTimes-2)

      l = where(gitmlats gt -45 and gitmlats lt 45 and dl gt 0.0, c)
      
      plotnum = 4
      get_position, ppp, space, sizes, plotnum, pos
      
      plot, gitmintrho(l), champintrho(l), psym = 3, $
            xrange = [0,maxi], yrange=[0,maxi], $
            xstyle = 1, ystyle = 1, $
            ytitle = 'CHAMP Rho (1e-12 kg/m3)', $
            xtitle = 'Gitm Rho (1e-12 kg/m3)', $
            pos = pos, /noerase
      oplot, [0,maxi], [0,maxi], thick = 3.0

      r = c_correlate(gitmintrho(l), champintrho(l), 0.0)
      rms = sqrt(mean((gitmintrho(l)-champintrho(l))^2))
      rmsn = sqrt(mean(((gitmintrho(l)-champintrho(l))/champintrho(l))^2))

      xyouts, maxi*0.05, maxi*0.90, '(E) Low Lat. LT: '+$
              string(mean(gitmlocaltime(l)),format='(f5.2)'), $
              charsize = 0.7
      xyouts, maxi*0.05, maxi*0.82, 'Corr: '+string(r,format='(f4.2)'), $
              charsize = 0.7
      xyouts, maxi*0.3, maxi*0.12, 'RMS: '+string(rms,format='(f6.3)')+$
              'x10!U-12!N kg/m!U3!N', charsize = 0.7
      xyouts, maxi*0.3, maxi*0.05, 'RMSN: '+string(r,format='(f4.2)'), $
              charsize = 0.7

      ;----------------------------------------------------------------
      ; Low Lats, Decending Phase
      ;----------------------------------------------------------------

      dl = gitmlats(1:nTimes-1)-gitmlats(0:nTimes-2)

      l = where(gitmlats gt -45 and gitmlats lt 45 and dl lt 0.0, c)
      
      plotnum = 5
      get_position, ppp, space, sizes, plotnum, pos
      
      plot, gitmintrho(l), champintrho(l), psym = 3, $
            xrange = [0,maxi], yrange=[0,maxi], $
            xstyle = 1, ystyle = 1, $
            ;ytitle = 'CHAMP Rho (1e-12 kg/m3)', $
            xtitle = 'Gitm Rho (1e-12 kg/m3)', $
            pos = pos, /noerase
      oplot, [0,maxi], [0,maxi], thick = 3.0

      r = c_correlate(gitmintrho(l), champintrho(l), 0.0)
      rms = sqrt(mean((gitmintrho(l)-champintrho(l))^2))
      rmsn = sqrt(mean(((gitmintrho(l)-champintrho(l))/champintrho(l))^2))

      xyouts, maxi*0.05, maxi*0.90, '(F) Low Lat. LT: '+$
              string(mean(gitmlocaltime(l)),format='(f5.2)'), $
              charsize = 0.7
      xyouts, maxi*0.05, maxi*0.82, 'Corr: '+string(r,format='(f4.2)'), $
              charsize = 0.7
      xyouts, maxi*0.3, maxi*0.12, 'RMS: '+string(rms,format='(f6.3)')+$
              'x10!U-12!N kg/m!U3!N', charsize = 0.7
      xyouts, maxi*0.3, maxi*0.05, 'RMSN: '+string(r,format='(f4.2)'), $
              charsize = 0.7

      closedevice

      psfile = 'champ_gitm_vs_'+ym+'.ps'

      setdevice, psfile, 'p', 5
      plotdumb

      ppp = 4
      space = 0.06
      pos_space, ppp, space, sizes

      ;----------------------------------------------------------------
      ; Low Lats, Ascending Phase vs alt
      ;----------------------------------------------------------------

      dl = gitmlats(1:nTimes-1)-gitmlats(0:nTimes-2)

      maxlat = 60.0

      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl gt 0.0, c)

      minalt = float(fix(min(champintalt(l))))
      maxalt = float(fix(max(champintalt(l)))+1)

      nA = maxalt-minalt + 1
      denchamp = fltarr(nA) - 9999.0
      stdchamp = fltarr(nA) - 9999.0
      dengitm = fltarr(nA) - 9999.0
      stdgitm = fltarr(nA) - 9999.0
      alts = findgen(nA)+minalt

      for i=0,nA-1 do begin
         a = minalt + i
         l = where(gitmlats gt -maxlat and gitmlats lt maxlat and $
                   dl gt 0.0 and $
                   champintalt ge a-1.0 and $
                   champintalt lt a+1.0, c)
         if (c gt 1) then begin
            denchamp(i) = mean(champintrho(l))
            stdchamp(i) = stdev(champintrho(l))
            dengitm(i) = mean(gitmintrho(l))
            stdgitm(i) = stdev(gitmintrho(l))
         endif
      endfor

      plotnum = 0
      get_position, ppp, space, sizes, plotnum, pos
      
      maxi = max([denchamp,dengitm])+max([stdchamp,stdgitm])
      mini = min([denchamp,dengitm])-max([stdchamp,stdgitm])
      dx = maxi-mini
      maxi = maxi + dx*0.05
      mini = mini - dx*0.05
      dy = maxalt - minalt
      maxalt = maxalt + 0.2*dy
      minalt = minalt - 0.05*dy

      plot, denchamp, alts, $
            pos = pos, /noerase, $
            thick = 4, $
            ystyle = 1, xstyle = 1, $
            yrange = [minalt, maxalt], $
            xrange = [mini, maxi], psym = -5, $
            ytitle = 'Altitude (km)', $
            xtitle = 'Rho (1e-12 kg/m3)'

      oplot, dengitm, alts, thick = 4, psym = -4, linestyle = 2

      for i=0,nA-1 do begin
         oplot, denchamp(i)+[-stdchamp(i),stdchamp(i)], $
                [alts(i),alts(i)]-0.1, thick = 2
         oplot, denchamp(i)+[stdchamp(i),stdchamp(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2
         oplot, denchamp(i)-[stdchamp(i),stdchamp(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2

         oplot, dengitm(i)+[-stdgitm(i),stdgitm(i)], $
                [alts(i),alts(i)]+0.1, thick = 2, linestyle=2
         oplot, dengitm(i)+[stdgitm(i),stdgitm(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2, linestyle=2
         oplot, dengitm(i)-[stdgitm(i),stdgitm(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2, linestyle=2
      endfor

      oplot, mini+[dx*0.1,dx*0.3], maxalt-[dy*0.08,dy*0.08], thick = 4
      xyouts, mini+dx*0.31, maxalt-dy*0.08, 'CHAMP', charsize=0.8
      oplot, mini+[dx*0.1,dx*0.3], maxalt-[dy*0.15,dy*0.15], $
             thick = 4, linestyle = 2
      xyouts, mini+dx*0.31, maxalt-dy*0.15, 'GITM', charsize=0.8
      xyouts, mini+dx*0.09, maxalt-dy*0.12, '(A)', charsize=0.8, $
              align = 1.0
 
      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl lt 0.0, c)
      xyouts, maxi-dx*0.09, maxalt-dy*0.12, 'LT: '+$
              string(mean(gitmlocaltime(l)),format='(f5.2)'), $
              charsize = 0.8, alignment = 1.0
 
      ;----------------------------------------------------------------
      ; Low Lats, Decending Phase vs alt
      ;----------------------------------------------------------------

      dl = gitmlats(1:nTimes-1)-gitmlats(0:nTimes-2)

      maxlat = 60.0

      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl lt 0.0, c)

      minalt = float(fix(min(champintalt(l))))
      maxalt = float(fix(max(champintalt(l)))+1)

      nA = maxalt-minalt + 1
      denchamp = fltarr(nA) - 9999.0
      stdchamp = fltarr(nA) - 9999.0
      dengitm = fltarr(nA) - 9999.0
      stdgitm = fltarr(nA) - 9999.0
      alts = findgen(nA)+minalt

      for i=0,nA-1 do begin
         a = minalt + i
         l = where(gitmlats gt -maxlat and gitmlats lt maxlat and $
                   dl lt 0.0 and $
                   champintalt ge a-1.0 and $
                   champintalt lt a+1.0, c)
         if (c gt 1) then begin
            denchamp(i) = mean(champintrho(l))
            stdchamp(i) = stdev(champintrho(l))
            dengitm(i) = mean(gitmintrho(l))
            stdgitm(i) = stdev(gitmintrho(l))
         endif
      endfor

      plotnum = 1
      get_position, ppp, space, sizes, plotnum, pos
      
      l2 = where(denchamp gt 0.0)
      l3 = where(dengitm gt 0.0)
      maxi = max([denchamp(l2),dengitm(l3)])+max([stdchamp(l2),stdgitm(l3)])
      mini = min([denchamp(l2),dengitm(l3)])-max([stdchamp(l2),stdgitm(l3)])
      dx = maxi-mini
      maxi = maxi + dx*0.05
      mini = mini - dx*0.05
      dy = maxalt - minalt
      maxalt = maxalt + 0.2*dy
      minalt = minalt - 0.05*dy

      plot, denchamp(l2), alts(l2), $
            pos = pos, /noerase, $
            thick = 4, $
            ystyle = 1, xstyle = 1, $
            yrange = [minalt, maxalt], $
            xrange = [mini, maxi], psym = -5, $
            xtitle = 'Rho (1e-12 kg/m3)'
      oplot, dengitm(l3), alts(l3), thick = 4, psym = -4, linestyle = 2

      for i=0,nA-1 do begin
         oplot, denchamp(i)+[-stdchamp(i),stdchamp(i)], $
                [alts(i),alts(i)]-0.1, thick = 2
         oplot, denchamp(i)+[stdchamp(i),stdchamp(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2
         oplot, denchamp(i)-[stdchamp(i),stdchamp(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2

         oplot, dengitm(i)+[-stdgitm(i),stdgitm(i)], $
                [alts(i),alts(i)]+0.1, thick = 2, linestyle=2
         oplot, dengitm(i)+[stdgitm(i),stdgitm(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2, linestyle=2
         oplot, dengitm(i)-[stdgitm(i),stdgitm(i)], $
                [alts(i)-0.25,alts(i)+0.25], thick = 2, linestyle=2
      endfor

      oplot, mini+[dx*0.1,dx*0.3], maxalt-[dy*0.08,dy*0.08], thick = 4
      xyouts, mini+dx*0.31, maxalt-dy*0.08, 'CHAMP', charsize=0.8
      oplot, mini+[dx*0.1,dx*0.3], maxalt-[dy*0.15,dy*0.15], $
             thick = 4, linestyle = 2
      xyouts, mini+dx*0.31, maxalt-dy*0.15, 'GITM', charsize=0.8
      xyouts, mini+dx*0.09, maxalt-dy*0.12, '(B)', charsize=0.8, $
              align = 1.0
 
      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl lt 0.0, c)
      xyouts, maxi-dx*0.09, maxalt-dy*0.12, 'LT: '+$
              string(mean(gitmlocaltime(l)),format='(f5.2)'), $
              charsize = 0.8, alignment = 1.0
 
      ;----------------------------------------------------------------
      ; Low Lats, Ascending Phase vs alt
      ;----------------------------------------------------------------

      f107int = interpolate_mine(gitmtime, f107, f107time)

      dl = gitmlats(1:nTimes-1)-gitmlats(0:nTimes-2)

      maxlat = 60.0

      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl gt 0.0, c)

      minf107 = float(fix(min(f107)))
      maxf107 = float(fix(max(f107))+1)

      nA = maxf107-minf107 + 1
      denchamp = fltarr(nA) - 9999.0
      stdchamp = fltarr(nA) - 9999.0
      dengitm = fltarr(nA) - 9999.0
      stdgitm = fltarr(nA) - 9999.0
      f107s = findgen(nA)+minf107

      for i=0,nA-1 do begin
         f = minf107 + i
         l = where(gitmlats gt -maxlat and gitmlats lt maxlat and $
                   dl gt 0.0 and $
                   f107int ge f-1.0 and $
                   f107int lt f+1.0, c)
         if (c gt 1) then begin
            denchamp(i) = mean(champintrho(l))
            stdchamp(i) = stdev(champintrho(l))
            dengitm(i) = mean(gitmintrho(l))
            stdgitm(i) = stdev(gitmintrho(l))
         endif
     endfor

      plotnum = 2
      get_position, ppp, space, sizes, plotnum, pos
      
      l2 = where(denchamp gt 0.0)
      l3 = where(dengitm gt 0.0)
      maxi = max([denchamp(l2),dengitm(l3)])+max([stdchamp(l2),stdgitm(l3)])
      mini = min([denchamp(l2),dengitm(l3)])-max([stdchamp(l2),stdgitm(l3)])
      dy = maxi-mini
      maxi = maxi + dy*0.2
      mini = mini - dy*0.05
      dx = maxf107 - minf107
      maxf107 = maxf107 + 0.05*dx
      minf107 = minf107 - 0.05*dx

      plot,f107s(l2), denchamp(l2), $
            pos = pos, /noerase, $
            thick = 4, $
            ystyle = 1, xstyle = 1, $
            xrange = [minf107, maxf107], $
            yrange = [mini, maxi], psym = -5, $
            xtitle = 'f107 (km)', $
            ytitle = 'Rho (1e-12 kg/m3)'
      oplot, f107s, dengitm, thick = 4, psym = -4, linestyle = 2

      for i=0,nA-1 do begin
         oplot, [f107s(i),f107s(i)]-0.1, $
                denchamp(i)+[-stdchamp(i),stdchamp(i)], thick = 2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                denchamp(i)+[stdchamp(i),stdchamp(i)], thick = 2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                denchamp(i)-[stdchamp(i),stdchamp(i)], thick = 2

         oplot, [f107s(i),f107s(i)]+0.1, $
                dengitm(i)+[-stdgitm(i),stdgitm(i)], $
                thick = 2, linestyle=2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                dengitm(i)+[stdgitm(i),stdgitm(i)], $
                thick = 2, linestyle=2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                dengitm(i)-[stdgitm(i),stdgitm(i)], $
                thick = 2, linestyle=2
      endfor

      oplot, minf107+[dx*0.1,dx*0.3], maxi-[dy*0.08,dy*0.08], thick = 4
      xyouts, minf107+dx*0.31, maxi-dy*0.08, 'CHAMP', charsize=0.8
      oplot, minf107+[dx*0.1,dx*0.3], maxi-[dy*0.15,dy*0.15], $
             thick = 4, linestyle = 2
      xyouts, minf107+dx*0.31, maxi-dy*0.15, 'GITM', charsize=0.8
      xyouts, minf107+dx*0.09, maxi-dy*0.12, '(C)', charsize=0.8, $
              align = 1.0
 
      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl gt 0.0, c)
      xyouts, maxf107-dx*0.09, maxi-dy*0.12, 'LT: '+$
              string(mean(gitmlocaltime(l)),format='(f5.2)'), $
              charsize = 0.8, alignment = 1.0
 
      ;----------------------------------------------------------------
      ; Low Lats, Ascending Phase vs alt
      ;----------------------------------------------------------------

      f107int = interpolate_mine(gitmtime, f107, f107time)

      dl = gitmlats(1:nTimes-1)-gitmlats(0:nTimes-2)

      maxlat = 60.0

      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl lt 0.0, c)

      minf107 = float(fix(min(f107)))
      maxf107 = float(fix(max(f107))+1)

      nA = maxf107-minf107 + 1
      denchamp = fltarr(nA) - 9999.0
      stdchamp = fltarr(nA) - 9999.0
      dengitm = fltarr(nA) - 9999.0
      stdgitm = fltarr(nA) - 9999.0
      f107s = findgen(nA)+minf107

      for i=0,nA-1 do begin
         f = minf107 + i
         l = where(gitmlats gt -maxlat and gitmlats lt maxlat and $
                   dl lt 0.0 and $
                   f107int ge f-1.0 and $
                   f107int lt f+1.0, c)
         if (c gt 1) then begin
            denchamp(i) = mean(champintrho(l))
            stdchamp(i) = stdev(champintrho(l))
            dengitm(i) = mean(gitmintrho(l))
            stdgitm(i) = stdev(gitmintrho(l))
         endif
      endfor

      plotnum = 3
      get_position, ppp, space, sizes, plotnum, pos
      
      l2 = where(denchamp gt 0.0)
      l3 = where(dengitm gt 0.0)
      maxi = max([denchamp(l2),dengitm(l3)])+max([stdchamp(l2),stdgitm(l3)])
      mini = min([denchamp(l2),dengitm(l3)])-max([stdchamp(l2),stdgitm(l3)])
      dy = maxi-mini
      maxi = maxi + dy*0.2
      mini = mini - dy*0.05
      dx = maxf107 - minf107
      maxf107 = maxf107 + 0.05*dx
      minf107 = minf107 - 0.05*dx

      plot,f107s(l2), denchamp(l3), $
            pos = pos, /noerase, $
            thick = 4, $
            ystyle = 1, xstyle = 1, $
            xrange = [minf107, maxf107], $
            yrange = [mini, maxi], psym = -5, $
            xtitle = 'f107 (km)', $
            ytitle = 'Rho (1e-12 kg/m3)'
      oplot, f107s(l3), dengitm(l3), thick = 4, psym = -4, linestyle = 2

      for i=0,nA-1 do begin
         oplot, [f107s(i),f107s(i)]-0.1, $
                denchamp(i)+[-stdchamp(i),stdchamp(i)], thick = 2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                denchamp(i)+[stdchamp(i),stdchamp(i)], thick = 2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                denchamp(i)-[stdchamp(i),stdchamp(i)], thick = 2

         oplot, [f107s(i),f107s(i)]+0.1, $
                dengitm(i)+[-stdgitm(i),stdgitm(i)], $
                thick = 2, linestyle=2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                dengitm(i)+[stdgitm(i),stdgitm(i)], $
                thick = 2, linestyle=2
         oplot, [f107s(i)-0.25,f107s(i)+0.25], $
                dengitm(i)-[stdgitm(i),stdgitm(i)], $
                thick = 2, linestyle=2
      endfor

      oplot, minf107+[dx*0.1,dx*0.3], maxi-[dy*0.08,dy*0.08], thick = 4
      xyouts, minf107+dx*0.31, maxi-dy*0.08, 'CHAMP', charsize=0.8
      oplot, minf107+[dx*0.1,dx*0.3], maxi-[dy*0.15,dy*0.15], $
             thick = 4, linestyle = 2
      xyouts, minf107+dx*0.31, maxi-dy*0.15, 'GITM', charsize=0.8
      xyouts, minf107+dx*0.09, maxi-dy*0.12, '(D)', charsize=0.8, $
              align = 1.0

      l = where(gitmlats gt -maxlat and gitmlats lt maxlat and dl lt 0.0, c)
      xyouts, maxf107-dx*0.09, maxi-dy*0.12, 'LT: '+$
              string(mean(gitmlocaltime(l)),format='(f5.2)'), $
              charsize = 0.8, alignment = 1.0
 
      closedevice

   endif

endfor

avediff = fltarr(24)
avediffNorm = fltarr(24)

for iTime=0,23 do begin

   l = where(gitmlocaltimeall gt itime-1 and gitmlocaltimeall le itime+1 and $
             gitmlatsall ge -45.0 and gitmlatsall le 45.0, c)
   if (c gt 0) then begin
      avediff(itime) = mean(gitmintrhoall(l) - champintrhoall(l))
      avediffNorm(itime) = mean(gitmintrhoall(l) - champintrhoall(l))/ $
                           mean(champintrhoall(l))
   endif

endfor

c_r_to_a, itime, min(gitmtimeall)

itime(1) = 1
itime(2) = 1
itime(3:5) = 0
c_a_to_r, itime, stime

DoY = (gitmtimeall - stime)/(24.0*3600.0)

avediffN = fltarr(52)
avediffS = fltarr(52)
avediffNNorm = fltarr(52)
avediffSNorm = fltarr(52)

for iTime=0,51 do begin

   l = where(DoY gt (itime-1)*7 and DoY le (itime+1)*7 and $
             gitmlatsall gt 60.0, c)
   if (c gt 0) then begin
      avediffN(itime) = mean(gitmintrhoall(l) - champintrhoall(l))
      avediffNNorm(itime) = mean(gitmintrhoall(l) - champintrhoall(l))/ $
                            mean(champintrhoall(l))
   endif

   l = where(DoY gt (itime-1)*7 and DoY le (itime+1)*7 and $
             gitmlatsall lt -60.0 and gitmlatsall gt -100.0, c)
   if (c gt 0) then begin
      avediffS(itime) = mean(gitmintrhoall(l) - champintrhoall(l))
      avediffSNorm(itime) = mean(gitmintrhoall(l) - champintrhoall(l))/ $
                            mean(champintrhoall(l))
   endif

endfor


end
