
pro fpi_read_makela_all, filelist, fpi, all

  nFiles = n_elements(filelist)
  for iFile = 0, nFiles-1 do begin
     file = filelist(iFile)
     fpi_read_makela, file, onefpi, oneall

     if (iFile eq 0) then begin
        fpi = onefpi 
        all = oneall
     endif else begin

        all2 = { $
               time: [all.time, oneall.time], $
               azimuth: [all.azimuth, oneall.azimuth], $
               zenith: [all.zenith, oneall.zenith], $
               temperature: [all.Temperature, oneall.Temperature], $
               velocity: [all.Velocity, oneall.Velocity], $
               FlagT: [all.flagT, oneall.flagT], $
               FlagV: [all.flagV, oneall.flagV]}

        fpi2 = {NorthTime  : [fpi.NorthTime, onefpi.NorthTime], $
                SouthTime  : [fpi.SouthTime, onefpi.SouthTime], $
                EastTime   : [fpi.EastTime, onefpi.EastTime], $
                WestTime   : [fpi.WestTime, onefpi.WestTime], $
                ZenithTime : [fpi.ZenithTime, onefpi.ZenithTime], $
                north      : [fpi.north, onefpi.north], $
                south      : [fpi.south, onefpi.south], $
                east       : [fpi.east, onefpi.east], $
                west       : [fpi.west, onefpi.west], $
                northTemp  : [fpi.northTemp, onefpi.northTemp], $
                southTemp  : [fpi.southTemp, onefpi.southTemp], $
                eastTemp   : [fpi.eastTemp, onefpi.eastTemp], $
                westTemp   : [fpi.westTemp, onefpi.westTemp], $
                northTflag : [fpi.northTflag, onefpi.northTflag], $
                southTflag : [fpi.southTflag, onefpi.southTflag], $
                eastTflag  : [fpi.eastTflag, onefpi.eastTflag], $
                westTflag  : [fpi.westTflag, onefpi.westTflag], $
                northVflag : [fpi.northVflag, onefpi.northVflag], $
                southVflag : [fpi.southVflag, onefpi.southVflag], $
                eastVflag  : [fpi.eastVflag, onefpi.eastVflag], $
                westVflag  : [fpi.westVflag, onefpi.westVflag], $
                fpitime    : [fpi.fpitime, onefpi.fpitime]}

        fpi = fpi2
        all = all2

     endelse
  endfor

end
