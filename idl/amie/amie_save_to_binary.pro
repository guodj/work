
filelist = findfile('b*.save')
display,filelist
if (n_elements(iFile) eq 0) then iFile = 0
iFile = fix(ask('file to transform to binary format',tostr(iFile)))

file = filelist(iFile)
restore,file

data(*,0,*,*) = data(*,0,*,*)*1000.0

c_r_to_a, itime, mean(time)
c_a_to_ymd, itime, ymd

FileOut = 'b'+ymd
amie_write_binary, fileout, fields, lats, mlts, time, data, $
                   imf = imf, ae = ae, dst = dst, hpi = hp, cpcp = cpcp, $
                   Version = Version

end
