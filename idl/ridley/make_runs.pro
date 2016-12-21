
vx =  [ 315.0, 315.0, 386.0, 480.0, 563.0, 480.0, 480.0, 480.0, 480.0, 690.0, 345.0, 400.0]
bz =  [  -0.1,  -0.1,  -0.1,  -0.1,  -0.1,  -6.0,  -3.0,   3.0,   6.0,  -0.1,  -0.1,   0.0]
n  =  [   3.0,   6.0,   8.0,  10.0,  18.0,  10.0,  10.0,  10.0,  10.0,   5.0,  20.0,   5.0]

vx = -1.0*vx
vy = vx*0.0
vz = vx*0.0

bx = bz*0.0
by = bz*0.0

t = vx*0.0 + 90000.0


ntimes = n_elements(vx)

spawn, 'mv run base'

openw,2,'run_all'
printf,2,'#!/bin/tcsh'
printf,2,''

for i=0,ntimes-1 do begin

  dirname = 'Run_'+chopr('0'+tostr(i),2)

  spawn, 'make rundir_rel'
  spawn, 'mv run '+dirname

  spawn, 'cp base/PARAM.in '+dirname

  openw,1,dirname+'/IMF'

  printf,1,''
  printf,1,'#SOLARWIND
  printf,1,strcompress(n(i))
  printf,1,strcompress(t(i))
  printf,1,strcompress(vx(i))
  printf,1,strcompress(vy(i))
  printf,1,strcompress(vz(i))
  printf,1,strcompress(bx(i))
  printf,1,strcompress(by(i))
  printf,1,strcompress(bz(i))
  printf,1,''

  close,1

  openw,1,dirname+'/input'
  printf,1,''
  close,1

  printf,2,'cd '+dirname
  printf,2,'Launch.pl < input > run_'+chopr('0'+tostr(i),2)+'.log'
  printf,2,'cd ..'
  printf,2,''

endfor

close,2
spawn,'chmod a+x run_all'

spawn,'mv base run'

end

