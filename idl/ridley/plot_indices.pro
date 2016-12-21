
run = 1

if run eq 1 then begin

  title = 'No Outflow - Inner Boundary = 10 cm!E-3!N'
  t1 = 18
  t3 = 25.0
  y = 40.0
  y2 = 50.0

endif

if run eq 2 then begin

  title = 'No Outflow - Inner Boundary = 100 cm!E-3!N'
  t1 = 18
  t3 = 35.0
  y = 60.0
  y2 = 70.0

endif

if run eq 3 then begin

  title = 'No Outflow - Inner Boundary = 1000 cm!E-3!N'
  t1 = 20
  t3 = 4980.0/60.0
  y = 65.0
  y2 = 72.0

endif

if run eq 4 then begin

  title = 'No Outflow - Inner Boundary = 25 - 250 cm!E-3!N'
  t1 = 4
  t3 = 43.0
  y = 68.0
  y2 = 75.0

endif

if run eq 5 then begin

  title = 'No Outflow - Inner Boundary = 100 - 1000 cm!E-3!N'
  t1 = 7
  t3 = 71.0
  y = 68.0
  y2 = 75.0

endif

filelist = findfile('*.save')
nfiles = n_elements(filelist)

if (nfiles gt 0) then begin

  file = filelist(nfiles-1)
  restore, file

  cpcp = indices(0,*,1,1) - indices(0,*,1,0)

  ntimes = n_elements(indices(0,*,1,1))
  time   = findgen(ntimes)/3.0+1.0/3.0

  setdevice, 'cpcp.ps', 'l',4,0.75

  plot, time, cpcp, xtitle = 'Time (minutes)', ytitle = 'Potential (kV)', $
          title = title, $
          xrange = [0,120], yrange = [0,100]

  oplot, [t3,t3], [0,1000]

  if (t1 lt 0) then begin
    loc = where(cpcp eq min(cpcp))
    t1 = loc(0)
  endif

  oplot, [time(t1),time(t1)], [0,1000]

  i = t1+10

  while (cpcp(i) gt cpcp(i-1)) do i = i + 1

  t2 = i-1

  oplot, [time(t2),time(t2)], [0,1000]

  dt = tostr(fix(time(t2) - time(t1)))+' min'
  oplot, [time(t1)+0.5, time(t2)-0.5], [y-1.0, y-1.0]
  oplot, [time(t1)+0.5, time(t1)+1.0], [y-1.0, y-2.0]
  oplot, [time(t1)+0.5, time(t1)+1.0], [y-1.0, y-0.0]
  oplot, [time(t2)-0.5, time(t2)-1.0], [y-1.0, y-2.0]
  oplot, [time(t2)-0.5, time(t2)-1.0], [y-1.0, y-0.0]
  xyouts, time(t1)+1, y+1, dt

  dt = tostr(fix(t3 - time(t1)))+' min'
  oplot, [time(t1)+0.5, t3-0.5], [y2-1.0, y2-1.0]
  oplot, [time(t1)+0.5, time(t1)+1.0], [y2-1.0, y2-2.0]
  oplot, [time(t1)+0.5, time(t1)+1.0], [y2-1.0, y2-0.0]
  oplot, [(t3)-0.5, (t3)-1.0], [y2-1.0, y2-2.0]
  oplot, [(t3)-0.5, (t3)-1.0], [y2-1.0, y2-0.0]
  xyouts, time(t1)+1, y2+1, dt


  print, time(t1),time(t2)

  closedevice

endif


end