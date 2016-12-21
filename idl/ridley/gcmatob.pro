
ntmax = 800

filein = ''
print, 'Enter ascii data file from TIEGCM run :'
read,filein

fileout = ''
print, 'Enter output binary file name : '
read,fileout

print, 'Reading ascii data'
rd_tgcm_ascii, filein, data, vars, utdes, axdes, xax, yax, ntmax

print, 'Writing binary file'
write_tgcm_bin, fileout, data, vars, utdes, axdes, xax, yax, ntmax

end