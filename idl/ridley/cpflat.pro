;  This procedure is designed to take output files from imf_dbdt.for and
;  find the times in a magnetometer file, then output 1 hour before the
;  events and 2 hours after the events.

pro cpflat

  yymm = ''
  print,'Enter YYMM of files to convert : '
  read,yymm
  loc = ''
  print,'Enter location of the greenland magnetometer files : '
  read,loc
  fileout = ''
  print,'Enter file to output : '
  read,fileout
  prewin = 0.0
  postwin = 0.0
  print,'Enter window size (in minutes), before & after event (pre,post) : '
  read, prewin, postwin
  prewin = prewin*60.0
  postwin = postwin*60.0
  outimfdbdt=''
  print, 'Enter name of output file from imf_dbdt.for :'
  read, outimfdbdt

  magfile = loc+'gma'+yymm
  setup_var_menu, 1, magfile, timeint, name, unit, nvar, nrows
  magfile = pad_string(magfile, 100)
  fileout = pad_string(fileout, 100)

  names = strarr(nvar(0)+1)
  units = strarr(nvar(0)+1)
  names(0:nvar(0)) = name(0,0:nvar(0))
  names = [names,'Delay']
  units(0:nvar(0)) = unit(0,0:nvar(0))
  units = [units,'seconds']
  source = strarr(nvar(0)+2)+'IDL program for IMF changes'

  fdnfile = 0.0
  fdids = fltarr(6,6)
  data = fltarr(100)
  inunit = 0.0

; open all of the data files :
;  first one, is the new one, which we have to create, brand new:

  outunit = ffcreate(fileout,names,units,source,fdnfile,fdids)

;  The next one is the old magnetometer data file:

  z = call_external('fopen_exe','fopen',		$
	magfile, fdnfile, fdids, inunit)

;  the final one is the outyymm.dat file, which contains the times of the
;    events which we want to take out of the magnetometer file :

  close,10
  openr,10,outimfdbdt

  line = ''
  uttime = dblarr(100)
  utsingle = double(0.0)
  itime = intarr(6)
  rows = lonarr(100,2)
  stringtime = strarr(100,2)
  row = 1.0
  noe = 0

  while not eof(10) do begin

    readf,10,line

    bl = byte(line)
    bmin = byte('0')
    bmax = byte('9')
    both = byte(' ')

    loc = where(((bl gt bmax(0)) or 		$
		(bl lt bmin(0))) and 		$
		(bl ne both(0)), count)

    if count eq 0 then begin

      itime(0) = fix(strmid(line,3,2))
      itime(1) = fix(strmid(line,8,2))
      itime(2) = fix(strmid(line,13,2))
      itime(3) = fix(strmid(line,18,2))
      itime(4) = fix(strmid(line,23,2))
      itime(5) = fix(strmid(line,28,2))
      
      print, itime

      c_a_to_r, itime, utsingle
      uttime(noe) = utsingle

      utsingle = utsingle - double(prewin)
      row = 0.0
      z = call_external('fsearch_exe','fsearch',		$
		fdnfile, fdids, inunit, utsingle, row)
      rows(noe,0) = long(row)
      loc = where(rows gt rows(noe,0), count)
      if count gt 0 then rows(noe,0) = max(rows(loc))+1

      utsingle = utsingle - 5.0*60.0
      c_r_to_a, midtime, utsingle
      c_a_to_s, midtime, tempstring
      stringtime(noe,0) = tempstring

      c_a_to_r, itime, utsingle
      utsingle = utsingle + double(postwin)
      row = 0.0
      z = call_external('fsearch_exe','fsearch',		$
		fdnfile, fdids, inunit, utsingle, row)
      rows(noe,1) = long(row)

      utsingle = utsingle + 5.0*60.0
      c_r_to_a, endtime, utsingle
      c_a_to_s, endtime, tempstring
      stringtime(noe,1) = tempstring

      noe = noe + 1

    endif

  endwhile

  close,10

  roww = 1.0

  que = ''

  print, 'Do you wish to create a command file to plot all of the data ?'
  read, que
  if (strmid(que,0,1) eq 'y') or (strmid(que,0,1) eq 'Y') then batchf = 1  $
  else batchf = 0

  if batchf eq 1 then begin

    print, 'Enter name of batch file : '
    filebatch = ''
    read, filebatch

    openw, 10, filebatch

    spawn, 'show def', def

    printf, 10, '$ set default '+strcompress(def(0),/remove_all)

  endif

  for i=0,noe-1 do begin

    print, 'writing ',stringtime(i,0)

    for j = rows(i,0), rows(i,1) do begin

;  actually reading from the data file, and writing to the new file

      rowr = float(j)

      z = call_external('fget_exe','fget',			$
	  fdnfile, fdids, inunit, data, utsingle, rowr)

      data(nvar(0)+3) = (utsingle - uttime(i))/60.0

      z = call_external('fput_exe','fput',			$
	  fdnfile,fdids,outunit,data,utsingle,roww)

      roww = roww + 1.0

    endfor

    if batchf eq 1 then begin

      printf, 10, '$ idl'
      printf, 10, 'batplot'
      printf, 10, '9          ;  magnetometer quick views'
      printf, 10, strcompress(fileout, /remove_all)
      printf, 10, ''
      printf, 10,							$
	'zasu$dkb500:[software.plotting.idl.data_files]quick_all_stat.dat'
      printf, 10, 'var'
      printf, 10, '1,22'
      printf, 10, '1,23'
      printf, 10, '1,24'
      printf, 10, '1,25'
      printf, 10, ''
      printf, 10, 'sm'
      printf, 10, '0.5'

      c_r_to_a, itime, uttime(i)
      sd = '0'+tostr(itime(2))
      sd = chopr(sd,2)
      sh = '0'+tostr(itime(3))
      sh = chopr(sh,2)
      sm = '0'+tostr(itime(4))
      sm = chopr(sm,2)
      psfile = strcompress(fileout+sd+sh+sm+'.ps',/remove_all)
      printf, 10, 'psname'
      printf, 10, psfile
      printf, 10, 'time'
      printf, 10, stringtime(i,0)
      printf, 10, stringtime(i,1)
      printf, 10, 'psl'
      printf, 10, 'exit'
      printf, 10, 'exit'

    endif

  endfor

  if batchf eq 1 then close, 10

  z = call_external('fclose_exe','fclose',fdnfile,fdids,inunit)
  z = call_external('fclose_exe','fclose',fdnfile,fdids,outunit)

end

