
PRO Get_Range_Int, maxminque, rangepl, numofvar, numofrang, 		$
	titlenam, rangescale, screen

  common crange, ranges, baseid

  rangepl = intarr(numofvar+1)

  if (maxminque eq 2) or (maxminque eq 4) then begin

      if (maxminque eq 1) or (maxminque eq 2) then 	$
	rs = 'Range' else rs = 'Scale'

      rangemenu = ['Choose Variables for Same '+rs+' ('+rs+' - 01) :',	$
		   '                  '+titlenam,       $
                   '                  Next '+rs+'                  ']

      numselected = 0
      currange    = 0
      i           = 0

      while numselected lt numofvar+1 do begin

        if n_elements(screen) eq 0 then 				$
	  ranum = wmenu(rangemenu, title=0, init=1)-1			$
	else begin
	  ranum = 0
	  while (ranum le 0) or (ranum gt n_elements(rangemenu)-1) do begin
	    print, rangemenu(0)
	    for j=1,n_elements(rangemenu)-1 do 				$
	      print, tostr(j),'. ',rangemenu(j)
	    print, 'Enter choice:'
	    read, ranum
	  endwhile
	  ranum = ranum - 1
	endelse

	if ranum eq numofvar+1 then begin

          if i ne 0 then begin

            i=0
	    currange = currange + 1
	    rangemenu(0) = 'Choose Variables for Same '+rs+' ('+rs+	$
		' - '+strcompress(string(currange+1),/remove_all)+') :'

	  endif

	endif else begin

	  rangemenu(ranum+1) = rs+' '+					$
		strcompress(string(currange+1),/remove_all)+		$
		'          '+titlenam(ranum)
	  numselected = numselected + 1
          i=1
	  rangepl(ranum) = currange

	endelse

      endwhile

      numofrang = currange

  endif else begin

    if (maxminque eq 1) or (maxminque eq 3) then			$
      numofrang = 0							$
    else numofrang = numofvar

  endelse

  if (maxminque le 2) or (maxminque eq 5) then rs = 'Range' 		$
  else rs = 'Scale'
  ranges = strarr(numofrang+1,2)

  if n_elements(screen) eq 0 then begin

    baseid = setup_ranges(rs,numofrang+1)
    widget_control, baseid.base, /realize
    xmanager, 'ranges', baseid.base, event_handler = 'mm_event', /modal

  endif else begin

    print, 'Have not completed at this time....'

  endelse

  for i=0,numofrang do rangescale(i,*) = ranges(i,*)

  RETURN

END
