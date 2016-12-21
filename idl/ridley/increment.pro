PRO Increment, times, addtimes

  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]

  times(0,0,*) = times(0,1,*)

  times(0,1,*) = times(0,1,*) + addtimes

  if times(0,1,5) gt 59 then begin
    times(0,1,5) = times(0,1,5) - 60
    times(0,1,4) = times(0,1,4) + 1
  endif

  if times(0,1,4) gt 59 then begin
    times(0,1,4) = times(0,1,4) - 60
    times(0,1,3) = times(0,1,3) + 1
  endif

  if times(0,1,3) gt 23 then begin
    times(0,1,3) = times(0,1,3) - 24
    times(0,1,2) = times(0,1,2) + 1
  endif

  eday = dayofmon
  if (times(0,1,1) eq 2) and (times(0,1,0) mod 4 eq 0) then eday(1) = eday(1)+1
  if times(0,1,2) gt eday(times(0,1,1)-1) then begin
    times(0,1,2) = times(0,1,2) - eday(times(0,1,1)-1)
    times(0,1,1) = times(0,1,1) + 1
  endif

  if times(0,1,1) gt 12 then begin
    times(0,1,1) = times(0,1,1) - 12
    times(0,1,0) = times(0,1,0) + 1
  endif

  if times(0,0,5) gt 59 then begin
    times(0,0,5) = times(0,0,5) - 60
    times(0,0,4) = times(0,0,4) + 1
  endif

  if times(0,0,4) gt 59 then begin
    times(0,0,4) = times(0,0,4) - 60
    times(0,0,3) = times(0,0,3) + 1
  endif

  if times(0,0,3) gt 23 then begin
    times(0,0,3) = times(0,0,3) - 24
    times(0,0,2) = times(0,0,2) + 1
  endif

  eday = dayofmon
  if (times(0,0,1) eq 2) and (times(0,0,0) mod 4 eq 0) then eday(1) = eday(1)+1
  if times(0,0,2) gt eday(times(0,0,1)-1) then begin
    times(0,0,2) = times(0,0,2) - eday(times(0,0,1)-1)
    times(0,0,1) = times(0,0,1) + 1
  endif

  if times(0,1,1) gt 12 then begin
    times(0,0,1) = times(0,0,1) - 12
    times(0,0,0) = times(0,0,0) + 1
  endif

  n = n_elements(times(*,0,0))
  if n ge 1 then for i=1,n-1 do times(i,*,*)=times(0,*,*)

  RETURN

END
	
