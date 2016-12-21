
;-----------------------------------------------------------------------------
;
;  Get File Name Widget
;
;-----------------------------------------------------------------------------

pro file_event, ev

  common tempfile, dirlid, fileid, dataid, dirid, path, filename, mesid,      $
		   flatlist2, dirlist2, maxfil, ext, ver, pre, pref

  widget_control, ev.id, get_uvalue=value

  if (n_elements(value) eq 0) then value = ''
  name = strmid(tag_names(ev, /structure_name), 7, 1000)
  newdir=0
  mess = ''

  widget_control, /hourglass

  case (name) of
  "BUTTON": begin
	      if (value eq "DONE") then begin

	        WIDGET_CONTROL, /destroy, ev.top
	        doneflag = 1
	        return

	      endif
            end
  "LIST": begin

            if (ev.id eq fileid) then begin

	      k=0
	      done = 0
	      iden = 0
	      while done eq 0 do begin

		if strlen(filename(k)) eq 0 then done = 1
		if filename(k) eq path(0)+flatlist2(ev.index) then begin
		  done = 1
		  iden = 1
		  mess = 'Identical file selected!!'
		endif
		if strmid(filename(k),0,1) eq ' ' then done = 1
		if k eq maxfil then begin
		  if done ne 1 then mess = 'Too Many Files!!'
		  done = 1
		endif
		if done eq 0 then k=k+1 

	      endwhile

	      if (iden ne 1) and (filename(k) eq '') then 	$
		filename(k) = path(0)+flatlist2(ev.index)

	    endif	

            if (ev.id eq dirlid) then begin

	      if !version.os eq 'vms' then begin

		if dirlist2(ev.index) eq '[-]' then diradd='-' else 	$
		diradd=dirlist2(ev.index)

		tpath=strmid(path(0),0,strlen(path(0))-1)+'.'+diradd+']'

	      endif else begin

		if dirlist2(ev.index) eq '.' then tpath = path(0)	$
		else if dirlist2(ev.index) eq '..' then begin
		  lslash = where(byte(path(0)) eq 47, nslash)
		  if nslash le 1 then tpath = '/'			$
		  else tpath = strmid(path(0),0,lslash(nslash-2)+1)
		endif else tpath = path(0) + dirlist2(ev.index) + '/'

	      endelse

	      newdir=1

	    endif

	  end

  else: begin

          if ev.id eq dirid then begin

	    newdir=1
	    tpath=''
	    widget_control, ev.id, get_value=tpath

	    print, tpath
            if (strmid(tpath(0),strlen(tpath(0))-1,1) ne '/') and	$
	       (!version.os ne 'vms') then tpath(0) = tpath(0)+'/'

          endif

	  idn = -1

	  for i=0, maxfil do if dataid(i) eq ev.id then idn = i

	  if idn ne -1 then begin

	    dummy = ''
	    widget_control, dataid(idn), get_value=dummy
	    filename(idn)=dummy

	  endif

        endelse

  endcase

  if newdir eq 1 then begin

    if !version.os eq 'vms' then dirlist = findfile(tpath(0)+'*.dir')	$
    else begin
      filelist = findfile('-Fa '+tpath(0))
      n = 0
      for i=0,n_elements(filelist)-1 do begin
        if strpos(filelist(i),'/') gt 0 then begin
         if n eq 0 then 						$
	    dirlist = [chopl(filelist(i),strlen(filelist(i))-1)]	$
	  else								$
	    dirlist = [dirlist,chopl(filelist(i),strlen(filelist(i))-1)]
	  n = n + 1
        endif
      endfor
    endelse
    loc = where((byte(ext) eq 44) or (byte(ext) eq 32), count)
    if count eq 0 then begin
;      if !version.os eq 'vms' then 					$
	flatlist = findfile(pref+tpath(0)+pre+ext+ver)	;		$
;      else spawn,'ls -d '+tpath(0)+pre+ext+ver,flatlist
    endif else begin
      tcount = 0
      for i=0,count do begin
	if i eq 0 then first = 0					$
	else first = loc(i-1)+1
	if i eq count then len = strlen(ext)-first			$
	else len = loc(i) - first
	fsearch = tpath(0)+pre+strmid(ext,first,len)+ver
;	if !version.os eq 'vms' then 					$
	  fl = findfile(pref+fsearch,count = nfiles)		;	$
;        else begin
;	  spawn,'ls -d '+fsearch,fl
;	  if n_elements(fl) gt 1 then count = n_elements(fl) 		$
;	  else if strlen(fl(0)) gt 0 then count = 1 else count = 0
;	endelse
        if (nfiles gt 0) or ((tcount eq 0) and (i eq count)) then begin
	  if tcount eq 0 then flatlist = fl 				$
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

      endif

      k = 0

      if !version.os eq 'vms' then begin
        dirlist2 = strarr(n_elements(dirlist)+1)
        dirlist2(0) = '[-]'

        for i=0, n_elements(dirlist)-1 do begin

          j=path_len

          while (j lt strlen(dirlist(i))) and 	$
                (strmid(dirlist(i),j,4) ne '.DIR') do j=j+1

          dirlist2(i+1) = strmid(dirlist(i),path_len,j-path_len)

        endfor

      endif else begin

	dirlist2 = dirlist
	path = tpath(0)

      endelse

    endif else begin

; this section of code is wrong!!!

      if !version.os eq 'vms' then dirlist = findfile(tpath(0)+'*.dir')	$
      else begin
        filelist = findfile('-Fa '+tpath(0))
        n = 0
        for i=0,n_elements(filelist)-1 do begin
          if strpos(filelist(i),'/') gt 0 then begin
           if n eq 0 then 						$
	      dirlist = [chopl(filelist(i),strlen(filelist(i))-1)]	$
	    else							$
	      dirlist = [dirlist,chopl(filelist(i),strlen(filelist(i))-1)]
	    n = n + 1
          endif
        endfor
      endelse

      if strlen(dirlist(0)) eq 0 then begin
        if !version.os eq 'vms' then dirlist2 = ['[-]'] 		$
	else dirlist2 = ['.','..']
        path = tpath(0)
      endif else mess = 'Invalid Path!!!'

    endelse

    if strlen(flatlist(0)) gt 0 then begin

      k = 0

      if !version.os eq 'vms' then begin

	j=1

	while (j lt strlen(flatlist(0))) and 	$
	      (strmid(flatlist(0),j,1) ne ']') do j=j+1

	path_len = j+1

      endif else begin

        loc = where(byte(flatlist(0)) eq 47, count)
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
	      el = strlen(exts)+1
	      k = path_len
	      while (k lt strlen(flatlist(i)))	$
	        and (strmid(flatlist(i),k,el) ne '.'+exts) do begin
		k=k+1
	      endwhile

	      if (k lt strlen(flatlist(i))) or	$
	         (j eq count) then found = 1

	      j = j + 1

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

    widget_control, dirid, set_value=path
    widget_control, dirlid, set_value=dirlist2
    widget_control, fileid, set_value=flatlist2

  endif

  widget_control, mesid, set_value=mess
  for i=0, maxfil do widget_control, dataid(i), set_value=filename(i)

end

