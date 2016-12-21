
if (n_elements(filelist) eq 0) then begin
   filelist = findfile("-t /*.save")
   if (strlen(filelist(0)) eq 0) then filelist = findfile("-t *.bin")
endif

filelist = ask('filename to plot',filelist(0))
filelist = findfile(filelist)

nfiles = n_elements(filelist)

if (nFiles eq 1 and strlen(filelist(0)) eq 0) then begin
   print, 'can not find file!'
   stop
endif

for iFile = 0, nFiles-1 do begin

   filename = filelist(iFile)

   print, 'Reading file ',filename

   read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
                           vars, data, rb, cb, bl_cnt, iTime, Version

   l = strpos(filename,'.')

   ncfile = strmid(filename,0,l+1)+'nc'
   print, 'nc file -->'+ncfile+'<--'

   id = ncdf_create('test.nc',/clobber)

   ncdf_control, id, /fill

   lonid = ncdf_dimdef(id, 'Longitude', nlons)
   latid = ncdf_dimdef(id, 'Latitude',  nlats)
   altid = ncdf_dimdef(id, 'Altitude', nalts)
   timeid = ncdf_dimdef(id, 'time', 1)

   nVars = n_elements(vars)

   TID = ncdf_vardef(id, 'Year', [timeid], /long)
   TID = [TID,ncdf_vardef(id, 'Month', [timeid], /long)]
   TID = [TID,ncdf_vardef(id, 'Day', [timeid], /long)]
   TID = [TID,ncdf_vardef(id, 'Hour', [timeid], /long)]
   TID = [TID,ncdf_vardef(id, 'Minute', [timeid], /long)]
   TID = [TID,ncdf_vardef(id, 'Second', [timeid], /long)]

   for iVar = 0, nVars-1 do begin

      x = byte(vars(iVar))
      n = n_elements(x)
      ni = 0
      i = 0
      while i lt n do begin
         if ((x(i) ge byte('a') and x(i) le byte('z')) or $
             (x(i) ge byte('A') and x(i) le byte('Z')) or $
             (x(i) ge byte('1') and x(i) le byte('9'))) then begin
            if (ni eq 0) then v = x(i) else v = [v,x(i)]
            ni++
         endif
         if (x(i) eq byte('+')) then $
            v = [v,byte('_'),byte('p'),byte('l'),byte('u'),byte('s')]
         if (x(i) eq byte('!')) then i++
         i++
      endwhile
      var = string(v)
      ;print, vars(iVar), '--->'+var+'<--'
      vid = ncdf_vardef(id, var, [lonid,latid,altid], /float)
      if (iVar eq 0) then VarID = [vid] else VarID = [VarID,vid]
   endfor

   ncdf_control, id, /endef

   for i=0,5 do ncdf_varput, id, TID(i), long(reform(itime(i)))

   for iVar = 0, nVars-1 do $
      ncdf_varput, id, VarID(iVar), float(reform(data(iVar,*,*,*)))

   ncdf_close,id

   spawn, 'mv test.nc '+ncfile

endfor


end

