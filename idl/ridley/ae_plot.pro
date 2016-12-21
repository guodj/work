
if (n_elements(file) eq 0) then begin
  filelist = findfile('sme*.dat')
  file = filelist(0)
endif
file = ask('filename to plot',file)

ae_read, file, time, ae

setdevice, 'ae.ps','p',4

titles = ['AE (nT)','AL (nT)','AU (nT)']

nVars = 3
nPts = n_elements(time)
data = fltarr(nVars,nPts)
data(0,*) = ae(*,0)
data(1,*) = ae(*,1)
data(2,*) = ae(*,2)
overplots = [0,0,0]

l = where(abs(data) gt 5000.0,c)
if (c gt 0) then data(l) = 0.0

plot_data, data, time, nvars, titles, overplots

closedevice

end
