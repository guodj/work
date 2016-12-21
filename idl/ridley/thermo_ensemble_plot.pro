

directories = findfile('-d Ccmc_Run0[24]/data Ccmc_Run1[159]/data')
nDirs = n_elements(directories)

nFilesMax = 1440

Files = strarr(nDirs, nFilesMax)

nFm = 0

for iDir = 0,nDirs-1 do begin

   filelist = findfile(directories(iDir)+'/'+'3DALL*.bin')
   nFiles = n_elements(filelist)
   if (nFiles gt nFm) then nFm = nFiles
   Files(iDir, 0:nFiles-1) = filelist

endfor

nFilesMax = nFm

file = Files(0,0)
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

r = 6372000.0 + alts(iAlt)*1000.0
dlat = lats(0,1) - lats(0,0)
dlon = lons(1,0) - lons(0,0)
area = dlat * dlon * r * cos(lats)
totalarea = total(area)

UseBaseDir = 1
if (n_elements(basedir) eq 0) then basedir = 'none'
basedir = ask('base directory to compare to (none for no comparison)',basedir)
if (strpos(basedir,'none') ne 0) then begin
   BaseFiles = findfile(basedir+'/3DALL*.bin')
   nFilesBase = n_elements(BaseFiles)
   BaseData = fltarr(nFilesBase, nLons, nLats)
   BaseGlobalMean = fltarr(nFilesBase)
   BaseTime = dblarr(nFilesBase)
endif

VarsToGet=[iVar]

AllData = fltarr(nDirs, nFilesMax, nLons, nLats)
MeanData = fltarr(nFilesMax, nLons, nLats)
StdData  = fltarr(nFilesMax, nLons, nLats)

AllTimes = dblarr(nFilesMax)

IsGood = intarr(nFilesMax)

for iT = 0, nFilesMax-1 do begin

   filelist = Files(*,iT)
   IsGoodTime = 1
   for iDir = 0,nDirs-1 do if (strlen(filelist(iDir)) eq 0) then IsGoodTime = 0

   if (IsGoodTime) then begin

      IsGood(iT) = 1

      for iDir = 0,nDirs-1 do begin

         file = Files(iDir,iT)
         if (iDir eq 0) then print, 'Reading File : ',file
         gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                             VarsToGet = VarsToGet
         if (iDir eq 0) then AllTimes(iT) = time
         AllData(iDir, iT, *, *) = reform(data(0,*,*,iAlt))

      endfor

      for iLon = 0,nLons-1 do for iLat=0,nLats-1 do begin
         MeanData(iT,iLon,iLat) = mean(AllData(*,iT,iLon,iLat))
         StdData(iT,iLon,iLat) = stdev(AllData(*,iT,iLon,iLat))
      endfor

   endif

endfor

for iT = 0, nFilesBase-1 do begin
   file = BaseFiles(iT)
   print, 'reading file ',file
   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet
   BaseGlobalMean(iT) = total(data(0,*,*,iAlt)*area)/totalarea
   BaseTime(it) = time
endfor


l = where(IsGood,nGood)
mini = min(StdData(l,*,*))
maxi = max(StdData(l,*,*))

globalmean = fltarr(nGood)
globalstd  = fltarr(nGood)

for iGood = 0, nGood-1 do begin
   globalmean(iGood) = total(MeanData(l(iGood),*,*)*area)/totalarea
   globalstd(iGood)  = total(StdData(l(iGood),*,*)*area)/totalarea
endfor

stime = min(AllTimes(l))
etime = max(AllTimes(l))

t = AllTimes(l)-stime

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice,'globalstd.ps','p',4
plot, t, globalmean, pos=[0.1,0.3,0.95,0.8], thick = 4, $
      xtickname = xtickname,			$
      xtitle = xtitle,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn, $
      yrange = [min(globalmean)-max(globalstd), $
                max(globalmean)+max(globalstd)], $
      ytitle = vars(iVar)

for iGood = 0, nGood-1 do begin
   y = globalmean(iGood) + globalstd(iGood)*[-1.0,1.0]
   oplot, t(iGood)*[1,1],y, thick = 2
endfor

if (UseBaseDir) then oplot, BaseTime-stime, BaseGlobalMean, thick = 4, linestyle = 2

closedevice

stop

; std plots

IsPercent = 1

;l = where(IsGood)
;MinMean = min(MeanData(l,*,*))
;if MinMean lt 0 then begin
;   IsPercent = 1
;endif else begin
;   StdData(l,*,*) = 100.0 * StdData(l,*,*) / meanData(l,*,*)
;   vars(iVar) = vars(iVar) + ' (% diff)'
;endelse

l = where(IsGood)
mini = min(StdData(l,*,*))
maxi = max(StdData(l,*,*))

psfile = 'test.ps'
maxrange = 40.0

for iFile = 0,nFilesMax-1 do begin

   if (IsGood(iFile)) then begin

      if (nFiles gt 1) then begin
         p = strpos(psfile,'.ps')
         if (p gt -1) then psfile = strmid(psfile,0,p)
         psfile_final = psfile+'_'+tostr(iFile,4)+'.ps'
      endif else begin
         psfile_final = psfile
      endelse

      value = reform(stdData(iFile,*,*))
      c_r_to_a, itime, AllTimes(iFile)
      c_a_to_s, itime, stime
      title = vars(iVar) + ' at '+stime

      print, 'Writing file ',psfile_final
      thermo_threeplot, psfile_final, $
                        value, AllTimes(iFile), lons, lats, mini, maxi, $
                        title, vars(iVar), maxrange

   endif

endfor

end

