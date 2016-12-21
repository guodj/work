
;filelist1 = findfile('3DALL_t12110[789]_*bin')  
;filelist2 = findfile('3DALL_t12111[012]_*bin')
;filelist = [filelist1,filelist2]
filelist = findfile('3DALL_t13031[012345]_??000?.bin')

month = 'nov'

nfiles = n_elements(filelist)

VarsToGet = [0,1,2,4,6,16,17]

for ifile = 0, nFiles-1 do begin

   file = filelist(iFile)

   gitm_read_bin_1var, file, data, time, nVars, Vars, version, $
                       VarsToGet = VarsToGet

   p = strpos(file, 'bin')
   if (p eq 0) then p = strlen(file)
   on2file = strmid(file,0,p)+'on2'

   print, on2file

;   setdevice, on2file+'.ps','p',5

   Alt = reform(data(2,0,0,*))

   iAlt300 = where(Alt gt 300.0*1000)
   iAlt300 = iAlt300(0)

   nLons = n_elements(data(0,*,0,0))
   nLats = n_elements(data(1,0,*,0))
   nAlts = n_elements(data(2,0,0,*))

   io_ = 3
   in2_ = 4
   iue_ = 5
   iun_ = 6

   o      = fltarr(nLons, nLats)
   n2     = fltarr(nLons, nLats)
   AltInt = fltarr(nLons, nLats)

   MaxValN2 = 1.0e21

   for iLon = 0, nLons-1 do begin
      for iLat = 0, nLats-1 do begin

         iAlt = nAlts-1
         Done = 0
         if (max(data(in2_,iLon,iLat,*)) eq 0.0) then Done = 1
         while (Done eq 0) do begin
            dAlt = (Alt(iAlt)-Alt(iAlt-1))
            n2Mid = (data(in2_,iLon,iLat,iAlt) + $
                     data(in2_,iLon,iLat,iAlt-1)) /2.0
            oMid  = (data(io_,iLon,iLat,iAlt) + $
                     data(io_,iLon,iLat,iAlt-1)) /2.0

            if (n2(iLon,iLat) + n2Mid*dAlt lt MaxValN2) then begin
               n2(iLon,iLat) = n2(iLon,iLat) + n2Mid*dAlt
               o(iLon,iLat)  =  o(iLon,iLat) +  oMid*dAlt
               iAlt = iAlt - 1
            endif else begin
               dAlt = (MaxValN2 - n2(iLon,iLat)) / n2Mid
               n2(iLon,iLat) = n2(iLon,iLat) + n2Mid*dAlt
               o(iLon,iLat)  =  o(iLon,iLat) +  oMid*dAlt
               AltInt(iLon,iLat) = Alt(iAlt) - dAlt
               Done = 1
            endelse
         endwhile

      endfor
   endfor

   ratio = n2*0.0
   loc = where(n2 gt 0.0,count)
   if (count gt 0) then ratio(loc) = o(loc)/n2(loc)
   loc = where(n2 eq 0.0,count)

   openw,2,on2file

   printf,2,'lon lat on2ratio Ue Un'

   printf,2, nLons-4, nLats-4

   for iLon = 2, nLons-3 do for iLat = 2, nLats-3 do begin

      printf,2,data(0,iLon,iLat,iAlt300)/!dtor, $
             data(1,iLon,iLat,iAlt300)/!dtor, $
             ratio(iLon,iLat), $
             data(iUe_,iLon,iLat,iAlt300), $
             data(iUn_,iLon,iLat,iAlt300)

   endfor

   close,2

endfor

end
