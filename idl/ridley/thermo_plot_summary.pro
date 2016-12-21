
filelist = findfile('??/data/summary.save')
nfiles = n_elements(filelist)

for iFile=0,nFiles-1 do begin

   restore, filelist(iFile)

   if (iFile eq 0) then begin
      nTimes = n_elements(time_all)
      time   = fltarr(nFiles, nTimes)
      rhon   = fltarr(nFiles, 5, nTimes)
      rhos   = fltarr(nFiles, 5, nTimes)
      tempn  = fltarr(nFiles, 5, nTimes)
      temps  = fltarr(nFiles, 5, nTimes)
      alts   = global_alt(0,0,*)/1000.0
      iAlts  = intarr(5)
      l = where(alts gt 120.0)
      iAlts(0) = l(0)
      l = where(alts gt 200.0)
      iAlts(1) = l(0)
      l = where(alts gt 300.0)
      iAlts(2) = l(0)
      l = where(alts gt 400.0)
      iAlts(3) = l(0)
      l = where(alts gt 500.0)
      iAlts(4) = l(0)
      
   endif

   time(iFile,*) = time_all
   for iAlt = 0,4 do begin
      rhos(iFile,iAlt,*)  = hem_rho(0,0,*,iAlts(iAlt))
      rhon(iFile,iAlt,*)  = hem_rho(0,1,*,iAlts(iAlt))
      temps(iFile,iAlt,*) = hem_temp(0,0,*,iAlts(iAlt))
      tempn(iFile,iAlt,*) = hem_temp(0,1,*,iAlts(iAlt))
   endfor

endfor


end
