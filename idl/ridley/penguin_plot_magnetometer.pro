dtmax = 24.0*3600.0

filein = ask('file to plot','today.csv')
penguin_read_data, filein, vars, data, time

iVars = [7,8,9]
over  = [0,0,0]

nV = n_elements(iVars)
nT = n_elements(data(0,*))

datatoplot = data(iVars,*)
Titles = Vars(ivars(where(over eq 0)))

if (strpos(filein,'today') gt -1) then psfile = 'today.ps' $
else psfile = 'magnetometer.ps'

print, psfile

setdevice, psfile,'p', 5
plot_data, datatoplot, time, nV, Titles, over, /day, dt=dtmax
closedevice

end
