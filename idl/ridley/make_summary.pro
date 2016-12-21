
if (n_elements(hour) eq 0) then hour = '01'
hour = ask('hour (01, 02, etc)',hour)

print, 'Reading data from : Run'+hour+'/data/'
filelist = findfile('Run'+hour+'/data/3DALL*.bin')

nFiles = n_elements(filelist)

for iFile = 0,nFiles-1 do begin

  File = filelist(iFile)
  print, 'Reading file :', File
  gitm_read_bin, File, data, time, nVars, Vars, version

  if (iFile eq 0) then begin
      dLon = data(0,1,0,0) - data(0,0,0,0)
      dLat = data(1,0,1,0) - data(1,0,0,0)
      r = cRe_ + reform(data(2,*,*,*))
      area = dLon * dLat * r^2 * cos(reform(data(1,*,*,*)))
      nAlts = n_elements(data(0,0,0,*))
      iMean_ = 0
      iStd_  = 1
      iMin_  = 2
      iMax_  = 3
      iGlobal_ = 0
      iSouth_  = 1
      iNorth_  = 2 
      Values = fltarr(nFiles,4,3,nVars,nAlts)
      lglobal = where(data(1,*,*,0) gt -!pi/2 and data(1,*,*,0) lt !pi/2)
      lNorth = where(data(1,*,*,0) gt !pi/4 and data(1,*,*,0) lt !pi/2)
      lSouth = where(data(1,*,*,0) gt -!pi/2 and data(1,*,*,0) lt -!pi/4)
  endif

  for iAlt = 0,nAlts-1 do begin
      for iVar = 0,nVars-1 do begin
          Var = reform(data(iVar,*,*,iAlt))
          are = reform(area(*,*,iAlt))
          Values(iFile,iMean_,iGlobal_,iVar,iAlt) = $
            total(Var(lGlobal) * are(lGlobal)) / total(are(lGlobal))
          Values(iFile,iMean_,iSouth_,iVar,iAlt) = $
            total(Var(lSouth) * are(lSouth)) / total(are(lSouth))
          Values(iFile,iMean_,iNorth_,iVar,iAlt) = $
            total(Var(lNorth) * are(lNorth)) / total(are(lNorth))

          Values(iFile,iMin_,iGlobal_,iVar,iAlt) = $
            min(Var(lGlobal))
          Values(iFile,iMin_,iSouth_,iVar,iAlt) = $
            min(Var(lSouth))
          Values(iFile,iMin_,iNorth_,iVar,iAlt) = $
            min(Var(lNorth))

          Values(iFile,iMax_,iGlobal_,iVar,iAlt) = $
            max(Var(lGlobal))
          Values(iFile,iMax_,iSouth_,iVar,iAlt) = $
            max(Var(lSouth))
          Values(iFile,iMax_,iNorth_,iVar,iAlt) = $
            max(Var(lNorth))

          Values(iFile,iStd_,iGlobal_,iVar,iAlt) = $
            stddev(Var(lGlobal))
          Values(iFile,iStd_,iSouth_,iVar,iAlt) = $
            stddev(Var(lSouth))
          Values(iFile,iStd_,iNorth_,iVar,iAlt) = $
            stddev(Var(lNorth))

      endfor
  endfor

endfor

save, file = './summary_'+hour+'.save', values, vars

end
