
if (n_elements(filein) eq 0) then begin
   filein = findfile('fismflux*.dat')
   filein = filein(0)
endif

filein = ask('file to plot',filein)

fism_read_input, filein, time, data
l = strpos(filein,'.dat')
file_front = strmid(filein,0,l)

if n_elements(nEnsembles) eq 0 then nEnsembles = 5
nEnsembles = fix(ask('number of ensembles',tostr(nEnsembles)))

if n_elements(uncertainty) eq 0 then uncertainty = 0.1
uncertainty = float(ask('uncertainty (0-1)',string(uncertainty)))

nTimes = n_elements(time)

mini = min(data(58,*)) * 0.8
maxi = max(data(58,*)) * 1.2

cUncertainty = string(uncertainty,format='(f4.2)')

setdevice, 'fism_'+cUncertainty+'.ps','l',5

stime = min(time)
plot, (time-stime)/3600.0, data(58,*), yrange = [mini,maxi], $
         xtitle = 'UT hours', ytitle = 'flux', thick = 4

nBins = 59
nTimes = n_elements(time)
spawn,'pwd',pwd

for i=0,nEnsembles-1 do begin

   d = data

   for iT = 0, nTimes-1 do begin
      v = 1.0 + uncertainty*randomn(s,nBins)
      l = where(v lt 0.1,c)
      if (c gt 0) then v(l) = 0.1
      l = where(v gt 2.0,c)
      if (c gt 0) then v(l) = 2.0
      d(*,iT) = v*data(*,iT)
   endfor

   oplot, (time-stime)/3600.0, d(58,*)

   file = file_front+'_'+cUncertainty+'_'+ $
          chopr('00'+tostr(i),3)+'.dat'

   print, 'writing file : ',file
   openw,1,file

   printf,1,''
   printf,1,'fism_emsemble.pro'
   printf,1,filein
   printf,1,pwd
   printf,1,''
   printf,1,'#START'

   for iT=0,nTimes-1 do begin
      c_r_to_a, itime, time(iT)
      printf,1,itime, d(*,iT),format='(6i5,'+tostr(nBins)+'e11.4)'
   endfor
   close,1

endfor

closedevice

end
