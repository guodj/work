
;
;	This procedure creates a widget which displays files and directories,
;	allowing the user to step through the directories, selecting
;	files which are on display.
;
;	datafile  - output - contains list of files selected
;		    input  - must be set up as a strarr(maxnumfil)
;	numoffil  - output - contains number of files selected
;	maxnumfil - input  - maximum number of files allowed to be selected
;			     6 absolute max
;	filtitle  - input  - 

pro get_files, datafile, numoffil, maxnumfil, filtitle, initpath, exten

  common fileblock, starttime, endtime, varcount, ftime 
  common tempfile, dirlid, fileid, dataid, dirid, path, filename, mesid,      $
		   flatlist2, dirlist2, maxfil, ext, ver, pre, pref

  if n_elements(maxnumfil) eq 0 then maxfil = 0 else maxfil = maxnumfil - 1
  if n_elements(initpath) eq 0 then initpath = ''
  if n_elements(filtitle) eq 0 then filtitle = ''
  if !version.os eq 'vms' then begin
    if n_elements(exten) eq 0 then ext = '*' else ext = exten
    pre = '*.'
    pref = ''
  endif else begin
    pref = '-d '
    if n_elements(exten) eq 0 then begin
      ext = '' 
      pre = ''
    endif else begin
      ext = exten
      pre = '*.'
    endelse
  endelse
  if !version.os eq 'vms' then ver = ';0' else ver=''

  if maxfil gt 5 then maxfil = 5
             
  if strlen(initpath) eq 0 then begin

    pwd = ['pwd']
    if !version.os eq 'vms' then spawn, '$ show def',paths		$
    else spawn, pwd, paths, /noshell

    if !version.os eq 'vms' then 					$
      paths(0) = strmid(paths(0),2,strlen(paths(0))-2)			$
    else paths(0) = paths(0) + '/'

  endif else begin

    paths = strarr(1)
    paths(0) = initpath

  endelse

  if !version.os eq 'vms' then dirlist = findfile(paths(0)+'*.dir')	$
  else begin
    filelist = findfile('-Fa')
    n = 0
    for i=0,n_elements(filelist)-1 do begin
      if strpos(filelist(i),'/') gt 0 then begin
        if n eq 0 then 							$
	  dirlist = [chopl(filelist(i),strlen(filelist(i))-1)]		$
	else								$
	  dirlist = [dirlist,chopl(filelist(i),strlen(filelist(i))-1)]
	n = n + 1
      endif
    endfor
  endelse
  loc = where((byte(ext) eq 44) or (byte(ext) eq 32), count)
  if count eq 0 then begin
;    if !version.os eq 'vms' then 					$
      flatlist = findfile(pref+paths(0)+pre+ext+ver)		;	$
;    else spawn, 'ls -d '+paths(0)+pre+ext+ver,flatlist
  endif else begin
    tcount = 0
    for i=0,count do begin
      if i eq 0 then first = 0						$
      else first = loc(i-1)+1
      if i eq count then len = strlen(ext) - first			$
      else len = loc(i) - first
      fsearch = paths(0)+pre+strmid(ext,first,len)+ver
;      if !version.os eq 'vms' then 					$
	fl = findfile(pref+fsearch, count = nfiles)		;	$
;      else begin
;	spawn, 'ls -d '+fsearch,fl
;	if n_elements(fl) gt 1 then count = n_elements(fl)		$
;	else if strlen(fl(0)) gt 0 then count = 1 else count = 0
;      endelse
      if (nfiles gt 0) or ((tcount eq 0) and (i eq count)) then begin
	if tcount eq 0 then flatlist = fl				$
	else flatlist = [flatlist,fl]
	tcount = 1
      endif
    endfor
  endelse

  if strlen(dirlist(0)) gt 0 then begin

    j=1


    if !version.os eq 'vms' then begin

      while (j lt strlen(dirlist(0))) and 	$
            (strmid(dirlist(0),j,1) ne ']') do j=j+1

      path_len = j+1
      path=strmid(dirlist(0),0,path_len)

      dirlist2 = strarr(n_elements(dirlist)+1)
      dirlist2(0) = '[-]'

      for i=0, n_elements(dirlist)-1 do begin

        j=path_len

        while (j lt strlen(dirlist(i))) and 	$
              (strmid(dirlist(i),j,4) ne '.DIR') do j=j+1

        dirlist2(i+1) = strmid(dirlist(i),path_len,j-path_len)

      endfor

    endif else dirlist2 = dirlist

  endif else begin

    if !version.os eq 'vms' then dirlist2 = ['[-]'] 			$
    else dirlist2 = ['.','..']
    path = paths(0)

  endelse

  if strlen(flatlist(0)) gt 0 then begin

    if !version.os eq 'vms' then begin

      j=1

      while (j lt strlen(flatlist(0))) and 	$
            (strmid(flatlist(0),j,1) ne ']') do j=j+1

      path_len = j+1

    endif else begin

      bf = byte(flatlist(0))
      loc = where(bf eq 47, count)
      if count gt 0 then path_len = loc(count-1)+1 else path_len = 0

    endelse

    flatlist2 = strarr(n_elements(flatlist))

    n = 0
    for i=0, n_elements(flatlist)-1 do begin

      j=path_len

      if !version.os eq 'vms' then begin

	if strlen(ext) gt 1 then begin

          while (j lt strlen(flatlist(i))) 			$
       	    and (strmid(flatlist(i),j,1) ne '.') do j=j+1

          flatlist2(i) = strmid(flatlist(i),path_len,j-path_len)

	endif else 						$
	  flatlist2(i) = strmid(flatlist(i),path_len,		$
			 strlen(flatlist(i))-path_len)

      endif else begin

	if strlen(ext) gt 1 then begin

	  found = 0

	  loc = where((byte(ext) eq 44) or (byte(ext) eq 32), count)
	  j = 0
	  while not found do begin

	    if j eq 0 then first = 0					$
	    else first = loc(j-1)+1
	    if j eq count then len = strlen(ext) - first		$
	    else len = loc(j) - first
	    exts = strmid(ext,first,len)
	    el = strlen(ext)+1
	    k = path_len
	    while (k lt strlen(flatlist(i)))	$
	      and (strmid(flatlist(i),k,el) ne '.'+exts) do k=k+1

	    if (k lt strlen(flatlist(i))) or	$
	       (j eq count) then found = 1
	  endwhile

          flattemp = strmid(flatlist(i),path_len,k-path_len)

	endif else							$
	  flattemp = strmid(flatlist(i),path_len,			$
		     strlen(flatlist(i))-path_len)

	found = 0
	j = 0
	while (not found) and (j lt n_elements(dirlist2)) do begin

	  if flattemp eq dirlist2(j) then found = 1
	  j = j + 1

	endwhile

	if not found then begin
	  flatlist2(n) = flattemp
	  n=n+1
	endif

      endelse

    endfor

  endif else flatlist2 = [' ']

  secondb = widget_base(title=filtitle)

  dataid = lonarr(6)

  font = '-adobe-times-bold-r-*-*-14-140*'

  for i=0, maxfil do begin

    dum = widget_label(secondb, value='File'+strcompress(string(i+1))+' : ', $
		  	yoff = i*50 + 12, font=font)
    dataid(i) = widget_text(secondb, /editable, xsize=50, ysize=1, 	$
			yoff = i*50 + 5, xoff = 65, /frame, font=font)

  endfor

  ybo = (maxfil+1)*50

  t8 = widget_base(secondb, /column, yoff=ybo, xsize=535)

  font = '-adobe-helvetica-bold-r-*-*-14-140*'

  t7 = widget_button(t8, value="Done", uvalue = "DONE", font=font)

  path = paths(0)

  font = '-adobe-times-bold-r-*-*-14-140*'

  dirid = widget_text(secondb, /editable, value=path, xsize=57, 	$
			ysize=1, yoff=ybo+50, font=font)

  t8 = widget_base(secondb, /column, yoff=ybo+90, xsize=265)
  font = '-adobe-times-bold-r-*-*-14-140*'
  b1 = widget_label(t8, value='Directories:',font=font)

  t38 = widget_base(secondb, /COLUMN, yoff=ybo+125,/frame, xsize=265)

  font = '-adobe-courier-medium-r-*-*-14-140*'
  dirlid  = WIDGET_LIST(t38, uvalue=dirlist2, value=dirlist2, 	$
			ysize=20, font=font)

  t8 = widget_base(secondb, /column, yoff=ybo+90, xsize=265, xoff=270)
  font = '-adobe-times-bold-r-*-*-14-140*'
  b1 = widget_label(t8, value='Files:',font=font)

  t38 = widget_base(secondb, /COLUMN, yoff=ybo+125, xoff=270,/frame, xsize=265)
  
  font = '-adobe-courier-medium-r-*-*-14-140*'
  fileid = WIDGET_LIST(t38, uvalue=flatlist2, value=flatlist2, 	$
			ysize=20,font=font)

  dum = widget_label(secondb, value='Messages : ', $
		      yoff = ybo+530)
  mesid = widget_text(secondb, xsize=40, ysize=1, 	$
		      yoff = ybo+520, xoff = 100, /frame)

  widget_control, secondb, /realize

  filename=strarr(6)

  XMANAGER, 'FILE', secondb,event_handler='FILE_EVENT', /MODAL   

  k=0
  datafile=strarr(6)

  for i=0,maxfil do begin

     if filename(i) ne '' then begin

	datafile(k)=filename(i)
	k=k+1

     endif

  endfor

  numoffil=k

  RETURN

END

