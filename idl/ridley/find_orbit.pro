
if (n_elements(stime) eq 0) then begin
    itime = [2001,05,04,12,00,00]
    c_a_to_s, itime, stime
endif

stime = ask('date',stime)

c_s_to_a, itime, stime
c_a_to_r, itime, rtime

y = chopr('0'+tostr(itime(0)),4)

ymd = chopr('0'+tostr(itime(0)),4) + $
  chopr('0'+tostr(itime(1)),2) + $
  chopr('0'+tostr(itime(2)),2)

file = findfile(y+'/*'+ymd+'*.cdf')

if (strlen(file(0)) gt 0) then begin
    read_cdf_pos, file(0), xyz, time

    d = abs(time - rtime)
    loc = where(d eq min(d))

    s = size(xyz)

    r = sqrt(xyz(0,loc(0),0)^2+xyz(1,loc(0),0)^2+xyz(2,loc(0),0)^2)
    if (r gt 300) then xyz = xyz/6372.0

    print, strmid(stime,0,15)
    print, format='(i4,a,i4,a,i4)', $
      round(xyz(0,loc(0),0)),',',round(xyz(1,loc(0),0)),',',$
            round(xyz(2,loc(0),0))

endif

end
