
dtmax = 24.0*3600.0

filein = ask('file to plot','today.csv')
penguin_read_data, filein, vars, data, time

iVars = [25,23,25,31]
over  = [ 0, 1, 1, 0]

nV = n_elements(iVars)
nT = n_elements(data(0,*))

datatoplot = data(iVars,*)
Titles = Vars(ivars(where(over eq 0)))

l = where(datatoplot(0,*) gt 0.0, c)
if (c gt 0) then datatoplot(0,l) = 0.0
l = where(datatoplot(1,*) gt 0.0, c)
if (c gt 0) then datatoplot(1,l) = 0.0
l = where(datatoplot(2,*) gt 0.0, c)
if (c gt 0) then datatoplot(2,l) = 0.0
l = where(datatoplot(3,*) gt 15.0, c)
if (c gt 0) then datatoplot(3,l) = 15.0

psfile = 'housekeeping.ps'

setdevice, psfile,'p', 5
plot_data, datatoplot, time, nV, Titles, over, /day, dt = dtmax
closedevice

end
