
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

if (n_elements(initial_guess) eq 0) then begin
   initial_guess = findfile('-t b*')
   iFile = 0
   help, initial_guess,iFile
   while ((strpos(initial_guess(iFile),"data") gt 0 or $
           strpos(initial_guess(iFile),"ps") gt 0 or $
           strpos(initial_guess(iFile),"sum") gt 0) and $
          iFile lt n_elements(initial_guess)-1) do begin
      iFile = iFile + 1
      print, iFile
   endwhile
   print, iFile
   initial_guess = initial_guess(iFile)
   if strlen(initial_guess) eq 0 then initial_guess='b970101'
endif

amie_file = ask('AMIE binary file name',initial_guess)
initial_guess = amie_file

if (n_elements(psfile) eq 0) then psfile = amie_file+'.ps'
psfile = ask('ps file',psfile)

read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp, version

nFields = n_elements(fields)
nLats = n_elements(lats)
nMlts = n_elements(mlts)
nTimes = n_elements(time)

;--------------------------------------------------------------
; Need to figure out what to plot, so list to fields to the user
;--------------------------------------------------------------

for i=0,nfields-1 do print, tostr(i+1)+'. '+fields(i)

;--------------------------------------------------------------
; Get field to be contoured
;--------------------------------------------------------------

type_1 = fix(ask('field to plot','1'))-1
if (type_1 lt 0) or (type_1 gt nfields-1) then type_1 = 0

;--------------------------------------------------------------
; Get start time and end time, with defaults as the file times
;--------------------------------------------------------------

c_r_to_a, itime_start, min(time)
c_r_to_a, itime_end, max(time)

c_a_to_s, itime_start, stime_start
c_a_to_s, itime_end, stime_end

start_time = ask('start time of plotting',strmid(stime_start,0,15))

if (strlen(start_time) lt 9) then $
  start_time = strmid(stime_start,0,9)+' '+start_time

;--------------------------------------------------------------
; I got sick of typing in the ending date, so if the date is
; the same, assume the user just wants to enter the time
;--------------------------------------------------------------

sdate = strmid(start_time,0,9)
if strpos(stime_end,sdate) gt -1 then 					$
  end_time_default = strmid(stime_end,10,5)				$
else end_time_default = strmid(stime_end,0,15)

end_time   = ask('end time of plotting',end_time_default)

;--------------------------------------------------------------
; If the user entered a short string, assume it is just a time
; and add the date on the front
;--------------------------------------------------------------

if (strlen(end_time) lt 9) then end_time = strmid(start_time,0,9)+' '+end_time

;--------------------------------------------------------------
; Now figure out where in the file these things are, with the
; default to give the user everything
;--------------------------------------------------------------

c_s_to_a, itime_start, start_time
c_a_to_r, itime_start, rtime_start

c_s_to_a, itime_end, end_time
c_a_to_r, itime_end, rtime_end

n_start = where(time ge rtime_start)
if n_start(0) ge 0 then n_start = n_start(0) else n_start = 0

n_end = where(time ge rtime_end)
if n_end(0) ge 0 then n_end = n_end(0) else n_end = n_elements(time)-1


nTimes = n_end-n_start+1
subdata = reform(data(n_start:n_end,type_1,0,*))

nlats = n_elements(lats)
nmlts = n_elements(mlts)

lat2d = fltarr(nTimes,nlats)
tim2d = dblarr(nTimes,nlats)
for i=0,nTimes-1 do lat2d(i,*) = lats
for i=0,nTimes-1 do tim2d(i,*) = time(i+n_start)-time(n_start)

subtime = time(n_start:n_end)
stime = time(n_start)
etime = time(n_end)

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

setdevice, psfile, 'p', 5

makect, 'bristow'
ncolors = 255

pos = [0.1, 0.1, 0.9, 0.5]

mini = 0.0
maxi = max(subdata)*1.1
levels = findgen(31)/30.0 * (maxi-mini) + mini

l = where(subdata lt levels(1),c)
if (c gt 0) then subdata(l) = levels(1)

yrange = [60,71]

contour, subdata, tim2d, lat2d, ystyle = 1, $
         xrange = [btr,etr], xstyle = 1, $
         ytitle = 'Latitude (deg)',		$
         xtickname = xtickname,			$
         xtitle = xtitle,			$
         xtickv = xtickv,			$
         xminor = xminor,			$
         xticks = xtickn,   $
         pos = pos, $
         yrange = yrange, $
         thick = 3, levels = levels, /fill

contour, subdata, tim2d, lat2d, ystyle = 5, $
         xrange = [btr,etr], xstyle = 5, $
         xtickname = xtickname,			$
         xtitle = xtitle,			$
         xtickv = xtickv,			$
         xminor = xminor,			$
         xticks = xtickn,   $
         pos = pos, $
         yrange = yrange, $
         thick = 3, nlevels = 5, /noerase

ctpos=pos
ctpos(0) = pos(2)+0.01
ctpos(2) = ctpos(0) + 0.05
plotct, ncolors, ctpos, mm(levels), fields(type_1), /right

pos = [0.1, 0.525, 0.9, 0.725]
plot, time-stime, hp(*,0), pos = pos, $
         xrange = [btr,etr], xstyle = 1, $
         xtickname = xtickname,			$
         xtickv = xtickv,			$
         xminor = xminor,			$
         xticks = xtickn,   $
         thick = 3, /noerase, ytitle = 'Hemispheric Power (GW)'

powerfile = 'power_2012.txt'
power_read, powerfile, powertime, power, North, nPts, notes

powertime = powertime-stime
l = where(powertime ge btr and powertime le etr and North eq 1)
oplot, powertime(l), power(l), psym = 4


pos = [0.1, 0.75, 0.9, 0.95]
plot, time-stime, cpcp/1000.0, pos = pos, $
      xrange = [btr,etr], xstyle = 1, $
      xtickname = xtickname,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      thick = 3, /noerase, $
      ytitle = 'CPCP (kV)'


closedevice

end
