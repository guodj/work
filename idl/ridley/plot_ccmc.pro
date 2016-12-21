
;filelist = findfile('*.ccmc')
filelist = findfile('*.ccmc')
nFiles = n_elements(filelist)

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, file
   shift = 0
   ndata = 1
   ccmc_results, file, time, lat, lon, data, shift, ndata

   c_r_to_a, itime, (time(0)+24.0*3600.0)
   yyyy = tostr(itime(0))
   ddd  = tostr(jday(itime(0), itime(1), itime(2)),3)

   radar = ''
   if (strpos(file,'electron') gt -1) then begin
      ccmcfile = '../CCMC.fism/observations/Eden/CHAMP.Eden.'+yyyy+'.'+ddd+'.dat'
      shift = 2
      ndata = 1
      ccmc_results, ccmcfile, ccmctime, ccmclat, ccmclon, ccmcdata,shift, ndata
      ytitle = 'Electron Density (#/m!U3!N)'
      radar = 'CHAMP'
   endif

   if (strpos(file,'density') gt -1) then begin
      ccmcfile = '../CCMC.fism/observations/Nden/CHAMP.Nden.'+yyyy+'.'+ddd+'.dat'
      shift = 3
      ndata = 3
      ccmc_results, ccmcfile, ccmctime, ccmclat, ccmclon, ccmcdata,shift, ndata
      ccmcdata = reform(ccmcdata(*,1))
      ytitle = 'Mass Density (kg/m!U3!N)'
      radar = 'CHAMP'
   endif

   if (strpos(file,'eis') gt -1) then radar = 'EISCAT_Svalbard'
   if (strpos(file,'msh') gt -1) then radar = 'Millstone_Hill'
   if (strpos(file,'jic') gt -1) then radar = 'Jicamarca_JRO'
   if (strpos(file,'pkr') gt -1) then radar = 'Poker_Flat'
   if (strpos(file,'stf') gt -1) then radar = 'Sondrestrom'

   if (strpos(file,'vi') gt -1) then begin
      ccmcfile = '../CCMC.fism/observations/Vdrift/'+radar+'.Vdrift.'+yyyy+'.'+ddd+'.dat'
      shift = 5
      ndata = 1
      ccmc_results, ccmcfile, ccmctime, ccmclat, ccmclon, ccmcdata,shift, ndata
      ccmcdata = reform(ccmcdata(*,0))
      ytitle = 'Vertical Drift (m/s)'
   endif

   if (strpos(file,'hmf2') gt -1) then begin
      ccmcfile = '../CCMC.fism/observations/HmF2/'+radar+'_HmF2.'+yyyy+'.'+ddd+'.dat'
      shift = 3
      ndata = 1
      ccmc_results, ccmcfile, ccmctime, ccmclat, ccmclon, ccmcdata,shift, ndata
      ccmcdata = reform(ccmcdata(*,0))
      ytitle = 'HmF2 (km)'
   endif

   if (strpos(file,'nmf2') gt -1) then begin
      ccmcfile = '../CCMC.fism/observations/NmF2/'+radar+'_NmF2.'+yyyy+'.'+ddd+'.dat'
      shift = 3
      ndata = 1
      ccmc_results, ccmcfile, ccmctime, ccmclat, ccmclon, ccmcdata,shift, ndata
      ccmcdata = reform(ccmcdata(*,0))
      ytitle = 'NmF2 (km)'
   endif

   if (n_elements(ccmctime) gt 1) then begin

      tm = time - ccmctime(0)
      l = where(tm ge 0,c)

      if (c gt 0) then begin

         time = time(l)
         lat = lat(l)
         lon = lon(l)
         data = data(l)

         outx = time
         inx = ccmctime
         iny = ccmcdata
         measure = interpolate_mine(outx, iny, inx)

         nPts = n_elements(time)
         
         for i=0,nPts-1 do begin
            d = abs(time(i) - ccmctime)
            if min(d) gt 600.0 then measure(i) = 1.0e32
         endfor

         diff = measure-data
         l = where(measure lt 1.0e31,c)
         if (c gt 0) then begin
            rms  = sqrt(mean(diff(l)^2))
            rmsd = sqrt(mean(measure(l)^2))
            pe = 1.0 - rms/rmsd
         endif

         psfile = file+'.ps'
         setdevice, psfile, 'p', 5

         stime = min(time)
         etime = max(time)
         time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn
         
         ppp = 2
         space = 0.07
         pos_space, ppp, space, sizes, ny = ppp
         
         plotdumb
         
         plotnum = 0
         get_position, ppp, space, sizes, plotnum, pos, /rect
         pos(0) = 0.1

         amax = max([abs(data),abs(measure(l))])*1.2
         if (min(data) lt 0.0) then amin = -amax else amin = 0.0

         dtime = time-stime
         plot, dtime, data, $
               xtickname = xtickname,			$
               xtitle = xtitle,			$
               xtickv = xtickv,			$
               xminor = xminor,			$
               xticks = xtickn,   $
               pos = pos, /noerase, thick = 4, $
               yrange = [amin,amax],ystyle=1, ytitle = ytitle, linestyle = 1
         
         oplot, dtime, measure, thick = 4, linestyle = 0, max_val = 1.0e31
         
         xyouts, etr*0.95, amax*1.02, radar, alignment = 1.0, charsize = 1.2
         
         oplot, [0.05,0.1]*etr, [0.93,0.93]*amax
         xyouts, 0.11*etr, 0.93*amax, radar
         oplot, [0.05,0.1]*etr, [0.86,0.86]*amax, linestyle = 1
         xyouts, 0.11*etr, 0.86*amax, 'GITM'
         
         plotnum = 1
         get_position, ppp, space, sizes, plotnum, pos, /rect
         pos(0) = 0.1
         
         amax = max(abs(diff(l)))*1.1
         
         plot, dtime, diff, $
               xtickname = xtickname,			$
               xtitle = xtitle,			$
               xtickv = xtickv,			$
               xminor = xminor,			$
               xticks = xtickn,   $
               pos = pos, /noerase, thick = 4, $
               yrange = [-amax,amax], ystyle = 1, ytitle = 'Diff in '+ytitle, $
               max_val = 1.0e31
         
         oplot, dtime, measure*0.0, thick = 4, linestyle = 1
         
         xyouts, etr*0.95, amax*0.85, 'PE : '+string(pe,format='(f6.3)'), $
                 alignment = 1.0
         xyouts, etr*0.95, amax*0.75, 'RMS : '+string(rms,format='(e13.4)'), $
                 alignment = 1.0
         
         closedevice

      endif else print, "No GITM points found to match data!!!!"

   endif

endfor

end
