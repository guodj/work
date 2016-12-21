
pro power_write, file, time, power, North, notes

  ; sample: '2011-01-01 00:38:02 NOAA-17 (S)  3   4.80   0.90'
  midline = ' NOAA-XX '

  openw,1,file

  nNotes = n_elements(notes)
  for i=0,nNotes-1 do printf,1,notes(i)
  printf,1,''

  nPts = n_elements(time)
  for i=0,nPts-1 do begin

     if (North(i)) then hem = '(N)' else hem = '(S)'

     c_r_to_a, itime, time(i)

     iHp = 2.09 * alog(power(i)+0.0001) * 1.0475
     if (iHp lt 1) then iHp = 1
     printf,1,itime(0),'-',itime(1),'-',itime(2),' ', $
              itime(3),':',itime(4),':',itime(5),midline,hem, $
              iHp, power(i),1.0, $
              format = '(i0,a,i02,a,i02,a,i02,a,i02,a,i02,a,a,i3,2f7.2)'

  endfor

  close,1

end
