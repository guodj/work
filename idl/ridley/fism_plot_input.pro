
if (n_elements(filein) eq 0) then begin
   filein = findfile('fismflux*.dat')
   filein = filein(0)
endif

filein = ask('file to plot',filein)

fism_read_input, filein, time, data

stime = min(time)
etime = max(time)

stime = min(time)
etime = max(time)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice, 'fism.ps','p',5

for i=0,58 do begin

   xtn = strarr(10)+' '
   xt = ' '
   if (i mod 10 eq 0) then plotdumb
   if (i mod 10 eq 9) then begin
      xtn = xtickname
      xt = xtitle
   endif

   frac = (i mod 10)/10.0
   top = 1.0-frac
   bot = top-0.09
   pos = [0.1, bot, 0.95, top]

   plot, (time-stime), data(i,*), ystyle = 1, xstyle = 1, $
         xtickname=xtn, xtitle=xt, xtickv=xtickv, $
         xminor=xminor, xticks=xtickn, thick = 4, $
         pos = pos, /noerase

endfor

closedevice

end
