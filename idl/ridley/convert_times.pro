
PRO convert_times, timeint, starttime, endtime, ftime, 			$
		   bftime, eftime, time, dt

  nfile = n_elements(timeint(*,0))

  starttime = timeint(0:nfile-1,0)
  endtime = timeint(0:nfile-1,1)

  dum1=intarr(6)
  dum2=intarr(6)
  mintime = double(1.0e32)
  maxtime = double(0.0)

  for i=0,nfile-1 do begin

    strtime = starttime(i)
    c_s_to_a, dum1, strtime
    c_a_to_r, dum1, tdum1

    strtime = endtime(i)
    c_s_to_a, dum2, strtime
    c_a_to_r, dum2, tdum2

    if tdum1 lt mintime then mintime=tdum1
    if tdum2 gt maxtime then maxtime=tdum2

  endfor

  c_r_to_a, dum1, mintime
  if (dum1(0) lt 1900) then dum1(0) = dum1(0) + 1900
  if (dum1(0) lt 1965) then dum1(0) = dum1(0) + 100
  c_r_to_a, dum2, maxtime
  if (dum2(0) lt 1900) then dum2(0) = dum2(0) + 1900
  if (dum2(0) lt 1965) then dum2(0) = dum2(0) + 100

;  dum1=[dum1(0:2),0,0,0]
;  dum2=[dum2(0:2),23,59,59]

  ftime=[dum1,dum2]

  c_a_to_r, ftime(0:5), bftime
  c_a_to_r, ftime(6:11), eftime
  dt = eftime - bftime

  dum1=[dum1(0:2),0,0,0]

  dumtime=[dum1,dum2]

  c_a_to_r, dumtime(0:5), bdumtime
  c_a_to_r, dumtime(6:11), edumtime

  diff_time, bdumtime, edumtime, dum3
  time=[dum1,dum2(0:2),0,0,0,dum3(0:2),0,0,0]

  return

end

