pro diff_time, time1, time2, timearray

  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]

  speryear = 31536000.0
  sperday  = 86400.0
  sperhour = 3600.0
  spermin  = 60.0
  timearray = intarr(5)

  diff = abs(time2-time1)

  timearray(0)=fix(diff/speryear)
  timearray(1)=fix((diff-timearray(0)*speryear)/sperday)
  timearray(2)=fix((diff-timearray(0)*speryear-timearray(1)*sperday)/sperhour)
  timearray(3)=fix((diff-timearray(0)*speryear-timearray(1)*sperday-	$
		    timearray(2)*sperhour)/spermin)
  timearray(4)=fix(diff-fix(diff/spermin)*spermin)

  RETURN

END
