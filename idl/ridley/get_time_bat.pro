pro get_time_bat, time

  stime = time(0:5)
  etime = time(6:11)

  c_a_to_s, stime, sst
  c_a_to_s, etime, set

  print, 'Current start time : ',sst
  print, 'Current end time : ',set

  nsst = ''
  nset = ''

  print, 'Enter new start time :'
  read, nsst
  if strlen(nsst) gt 0 then c_s_to_a, stime, nsst

  print, 'Enter new end time :'
  read, nset
  if strlen(nset) gt 0 then c_s_to_a, etime, nset

  c_a_to_r, stime, sr
  c_a_to_r, etime, er

  diff_time, sr, er, dt

  time = [stime, etime, dt]

  return

end
