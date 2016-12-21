
if (n_elements(filelist) eq 0) then filelist = findfile("-t *.bin")

filelist = ask('filename to plot',filelist(0))
filelist = findfile(filelist)

nfiles = n_elements(filelist)

if (nFiles eq 1 and strlen(filelist(0)) eq 0) then begin
   print, 'can not find file!'
   stop
endif

if (n_elements(psfile) eq 0) then psfile = filelist(0)+'.ps'
psfile = ask('ps file name',psfile)

diffdir = ask('directory to use to make differences (return for no diffs)','')

file = filelist(0)
VarsToGet=[0,1,2]
gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                        VarsToGet = VarsToGet

alts = reform(data(2,0,0,*)) / 1000.0
lats = reform(data(1,*,*,0))
lons = reform(data(0,*,*,0))
nAlts = n_elements(alts)
nLons = n_elements(lons(*,0))
nLats = n_elements(lats(0,*))

CalcTEC = 0
display, vars
print,tostr(nVars),'. TEC'
if (n_elements(iVar) eq 0) then iVar = '9' else iVar = tostr(iVar)
iVar = fix(ask('which var to plot',iVar))
if (iVar eq nVars) then begin
   CalcTEC = 1
   iVar = 33
endif

for i=0,nalts-1 do print, tostr(i)+'. '+string(alts(i))
if (n_elements(iAlt) eq 0) then iAlt='0' else iAlt=tostr(iAlt)
iAlt = fix(ask('which altitude to plot',iAlt))

AllValues = fltarr(nFiles, nLons, nLats)
AllTimes = dblarr(nFiles)

VarsToGet=[iVar]

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)

   print, 'Reading file ',filename

   file = filelist(iFile)
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet
   AllTimes(iFile) = time

   if (CalcTEC eq 1) then begin
      tec = fltarr(nLons,nLats)

      gitm_read_bin_1var, file, He, timeHe, nVars, Vars, version, $
                          VarsToGet = [31,32]

      for iAlt = 2,nAlts-3 do begin
         tec = tec + (data(0,*,*,iAlt)-He(0,*,*,iAlt)-He(1,*,*,iAlt)) * $
               (alts(iAlt)-alts(iAlt-1))*1000.0
      endfor
      AllValues(iFile,*,*) = tec/1e16
   endif else AllValues(iFile,*,*) = data(0,*,*,iAlt)

   if (strlen(diffdir) gt 0) then begin
      l1 = strpos(file,'_t')
      l2 = strpos(file,'.bin')
      filename = diffdir+'/'+strmid(file,l1-5,l2-2)+'??.bin'
      gitm_read_bin_1var, filename, data_diff, time, nVars, Vars, version, $
                          VarsToGet = VarsToGet
      if (CalcTEC eq 1) then begin
         tec = fltarr(nLons,nLats)
         for iAlt = 2,nAlts-3 do begin
            tec = tec + data_diff(0,*,*,iAlt) * $
                  (alts(iAlt)-alts(iAlt-1))*1000.0
         endfor
         AllValues(iFile,*,*) = AllValues(iFile,*,*) - tec/1e16
      endif else AllValues(iFile,*,*) = AllValues(iFile,*,*) - data_diff(0,*,*,iAlt)
   endif

endfor

determine_min_max, AllValues, mini, maxi

if (CalcTEC) then colortitle = 'TEC (TECU)' else colortitle = Vars(iVar)
maxrange = 40.0

for iFile = 0,nFiles-1 do begin

   if (nFiles gt 1) then begin
      p = strpos(psfile,'.ps')
      if (p gt -1) then psfile = strmid(psfile,0,p)
      psfile_final = psfile+'_'+tostr(iFile,4)+'.ps'
   endif else begin
      psfile_final = psfile
   endelse

   value = reform(AllValues(iFile,*,*))
   c_r_to_a, itime, AllTimes(iFile)
   c_a_to_s, itime, stime
   title = colortitle + ' at '+stime

   print, 'Writing file ',psfile_final
   thermo_threeplot, psfile_final, $
                     value, AllTimes(iFile), lons, lats, mini, maxi, $
                     title, colortitle, maxrange

endfor

end

