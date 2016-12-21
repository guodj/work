
acedir = getenv('ACE_DATA')

if (n_elements(year) eq 0) then  year  = '2000'
if (n_elements(month) eq 0) then month = '01'
if (n_elements(day) eq 0) then   day   = '01'

year  = tostr(fix(ask('year', string(year))))
month = chopr('0'+tostr(fix(ask('month',string(month)))),2)
day   = chopr('0'+tostr(fix(ask('day',  string(day)))),2)

file = acedir+'/'+year+'/ace'+year+month+day+'.save'

restore, file

if (n_elements(average) eq 0) then average = '-1'
average = fix(ask('time to average over (-1 for no averaging)', $
          strcompress(string(average))))

fileout = 'imf'+year+month+day+'.dat'
fileout = ask('output file',fileout)

if (average le 0) then begin
  ntimes = n_elements(time)
endif else begin
  stime = min(time)
  etime = max(time)
  c_r_to_a, istime, stime
  c_r_to_a, ietime, etime
  istime(3:5) = 0
  ietime(3)   = 23
  ietime(4)   = 59
  ietime(5)   = 59
  c_a_to_r, istime, stime
  c_a_to_r, ietime, etime
  ntimes = fix((etime - stime)/average) + 1
endelse

openw,1,fileout

printf,1,''
printf,1,'Data taken from file : ', file
printf,1,''
printf,1,'#ZEROBX'
printf,1,'T'
printf,1,''
printf,1,'#TIMEDELAY'
printf,1,delay
printf,1,''
if (isgse) then begin
    printf,1,'#COOR'
    printf,1,'GSE'
    printf,1,''
endif
printf,1,'#START'

print, 'ntimes : ',ntimes

for it=0,ntimes-1 do begin

    if (average le 0) then begin
        c_r_to_a, itime, time(it)
        i = it
    endif else begin
        currenttime = stime + double(average) * double(it)
        c_r_to_a, itime, currenttime
    endelse

    year   = itime(0)
    month  = itime(1)
    day    = itime(2)
    hour   = itime(3)
    minute = itime(4)
    second = itime(5)
    milli  = 0

    if (average le 0) then begin

        loc = it

    endif else begin

        loc = where(time ge currenttime - average/2 and $
                    time lt currenttime + average/2, count)

        if (count eq 0) then begin
            
            d = abs(time - currenttime)
            loc = where(d eq min(d))

        endif

    endelse

    bx   = mean(b(0,loc))
    by   = mean(b(1,loc))
    bz   = mean(b(2,loc))
    vx   = mean(v(0,loc))
    vy   = mean(v(1,loc))
    vz   = mean(v(2,loc))
    ; vz and vy should not be within 50% of vx.
    if ((abs(vx)-abs(vy))/abs(vx) lt 0.5) then vy = 0.0
    if ((abs(vx)-abs(vz))/abs(vx) lt 0.5) then vz = 0.0
    den  = mean(n(loc))
    temp = mean(t(loc))

    printf,1, format = '(i5,5i3,i4,3f8.2,4f9.2,f10.1)', $
              year, month, day, hour, minute, second, milli, $
              bx, by, bz, vx, vy, vz, den, temp

endfor

close,1


end


