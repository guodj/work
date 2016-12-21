;-----------------------------------------------------------------------------
;
;  Main Widget
;    Controls all aspects of the plot, such as :
;	- Plot times
;	- Title
;	- Range and Scales
;	- Vector Lengths
;	- Velocities for Magnetometer plots
;	- Max Latitude for Dwell plots
;
;-----------------------------------------------------------------------------

pro widplote, ev, screen = screen

  common timeblock, sid, time, bftime, eftime
  common varblock, nfile, lid, titles, units, lfn
  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common rep_data, change
  common type, plottype, types
  common data, single, t, dt, basetime, numofvar, numofpts
  common text_id, intext, inps, instring, inmaxx, inmaxy, inminx, inminy
  common control, getplottype, getfiles, stoppro
  common flag_structure, flags
  common comdis, displays

  if n_elements(screen) eq 0 then begin

    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 1000)

  endif else begin

    name = 'BUTTON'
    value = ev.value

  endelse

  case (name) of

  "BUTTON": begin

	if value eq 'PLOT' then begin
	  pvalue = 'PLOT'
	  value = 'PLOT'
	endif
	if value eq 'ANIMATE' then begin
	  pvalue = 'ANIMATE'
	  value = 'PLOT'
	endif
	if value eq 'Make PS File (landscape)' then begin
	  pvalue = 'PSL'
	  value = 'PLOT'
	endif
	if value eq 'Make EPS File (landscape)' then begin
	  pvalue = 'EPSL'
	  value = 'PLOT'
	endif
	if value eq 'Make PS File (portrait)' then begin
	  pvalue = 'PSP'
	  value = 'PLOT'
	endif
	if value eq 'Make EPS File (portrait)' then begin
	  pvalue = 'EPSP'
	  value = 'PLOT'
	endif
	if value eq 'Print plot to oosik (landscape)' then begin
	  pvalue = 'OL'
	  value = 'PLOT'
	endif
	if value eq 'Print plot to oosik (Portrait)' then begin
	  pvalue = 'OP'
	  value = 'PLOT'
	endif
        if value eq 'SFLD' then begin
	  pvalue = value
	  value = 'PLOT'
	  if n_elements(psfile) gt 0 then psfile_old = psfile 	$
	  else psfile_old = ''
	  if n_elements(screen) eq 0 then begin
	    get_files, dumfile,n,1,'Enter output file name'
	    psfile = dumfile(0)
	  endif else begin
	    print, 'Enter output file name : '
	    psfile = ''
	    read, psfile
	  endelse
	endif

        case (value) of

	  "DONE" : begin

	     stoppro = 1

	     WIDGET_CONTROL, /destroy, ev.top
	     return

	   end

	  "New File" : begin

	     stoppro = 0
	     getfiles = 1

	     WIDGET_CONTROL, /destroy, ev.top
	     return

	   end

	  "New Plot Type" : begin

	     stoppro = 0
	     getplottype = 1

	     WIDGET_CONTROL, /destroy, ev.top
	     return

	   end

;	  "SFLD"  : compute, grafvar, titlenam, unitnam, times, value
	  "SELF"  : compute, grafvar, titlenam, unitnam, times, value
	  "FFT"   : compute, grafvar, titlenam, unitnam, times, value

	  "EXTRAS" : extra

	  "PLOT" : begin

            c_a_to_r, time(0:5), btime
            c_a_to_r, time(6:11), etime

	    done = 0
	    ani = 0
            if (pvalue ne 'PLOT') and (pvalue ne 'ANIMATE') and		$
	      (strmid(pvalue,0,1) ne 'O') then begin
	      pslen = strlen(psfile)
	      if pslen eq 0 then done = 1
	    endif 
	    if btime eq etime then begin
	      value='Times are identical... can not continue'
	      done = 1
	    endif
	    if numofvar lt 0 then begin
	      value = 'No Variables selected... can not continue'
	      done = 1
	    endif

             if ((btime lt bftime) and (etime lt bftime)) or	      $
                ((btime gt eftime) and (etime gt eftime)) then begin
	       value = 'Time period Selected is out of Range'	
	       done = 1
	     endif

	     if done eq 0 then begin

               get_var_and_pos, grafvar, titlenam, unitnam, times, screen 
	       if pvalue eq 'PLOT' then begin
	         set_plot, 'X'
		 xsize = 1.33*float(change.y_win_size)
		 if !d.window eq -1 then					$
	           window, 1, xsize = xsize, ysize = change.y_win_size
	       endif else begin
		 case (pvalue) of
		   'PSL' : setdevice, psfile, 'landscape', psfont, 	$
				percent
		   'PSP' : setdevice, psfile, 'portrait', psfont, 	$
				percent
		   'EPSL' : setdevice, psfile, 'landscape', psfont, 	$
				percent, /eps
		   'EPSP' : setdevice, psfile, 'portrait', psfont, 	$
				percent, /eps
		   'OL' : begin
		     get_temp_ps, tempps
		     setdevice, tempps, 'landscape', psfont, percent
		   end
		   'OP' : begin
		     get_temp_ps, tempps
		     setdevice, tempps, 'portrait', psfont, percent
		   end
		   'ANIMATE' : ani = 1
		   else :
		 endcase
	       endelse

	       !P.thick = 1.0

	       if n_elements(screen) eq 0 then widget_control, /hourglass

	       if pvalue eq 'SFLD' then begin
		 plottype_old = plottype
		 plottype = -1
	       endif
	       plotvar, grafvar, titlenam, unitnam, times, screen, ani
	       if pvalue eq 'SFLD' then plottype = plottype_old
               if draw eq 1 then 					$
		 draw_final, final_x, final_y, ext.opos
	       draw_arr_final, arrow_x, arrow_y, ext.opos
	       draw_arr_final, ext.line_x, ext.line_y, ext.opos, 1
               str_put, str_x, str_y, str_s, str_o, str_l, ext.opos

	       if (pvalue ne 'PLOT') and 		$
		  (ani eq 0) and 			$
		  (pvalue ne 'SFLD') then begin
		 device, /close
		 set_plot, 'X'
		 plotwhere = 2
		 if (pvalue eq 'OP') or (pvalue eq 'OL') then		$
		   spawn, 'print /que=poosik '+tempps, list
	       endif

	     endif

	     if value ne 'PLOT' then print, value

             if pvalue eq 'SFLD' then psfile=psfile_old

	  end

	  "CHANGE RANGE" : change_range

	  "SELECT"                               : begin

	    if checkselect eq 1 then begin

	      if n_elements(screen) eq 0 then			$
		widget_control, ev.id, set_value='Select all'

	      for i=0,nfile-1 do begin

		k = 0

	        for j=0,n_elements(titles(0,*))-1 do begin

		  if strlen(titles(i,j)) gt 0 then begin

	            titles(i,j) = strmid(titles(i,j),0,14)
		    k=k+1

                  endif

		endfor

		if n_elements(screen) eq 0 then begin
		  n = lid(i)
		  widget_control, n, set_value = titles(i,0:k-1)
		endif

	      endfor

	      change.data = 1
	      change.place = 1
	      numofvar = -1
              checkselect = 0

	    endif else begin

	      if n_elements(screen) eq 0 then			$
		widget_control, ev.id, set_value='Deselect all'

	      numofvar = -1

	      for i=0,nfile-1 do begin

		k = 0

	        for j=0,n_elements(titles(0,*))-1 do begin

		  if strlen(titles(i,j)) gt 0 then begin

	            titles(i,j) = strmid(titles(i,j),0,14)+     $
		      '(selected)'
	            numofvar = numofvar + 1
		    k=k+1

                  endif

		endfor

	        if n_elements(screen) eq 0 then	begin
		  n = lid(i)
	          widget_control, n, set_value = titles(i,0:k-1)
		endif

	      endfor

	      change.data = 1
	      change.place = 1
              checkselect = 1

            endelse

	  end

	  "INVAR"                               : begin

	    if checkinvar eq 2 then begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Add Invariant Latitudes'
	      checkinvar = 1
	    endif else begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Take off Invariant Latitudes'
	      checkinvar = 2
	      change.invar = 1
            endelse

            change.data = 1

	  end

	  "LOGX"                               : begin

	    if flags.logx eq 1 then begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Log X Axis'
	      flags.logx = 0
	    endif else begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Linear X Axis'
	      flags.logx = 1
            endelse

	  end

	  "LOGY"                               : begin

	    if flags.logy eq 1 then begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Log Y Axis'
	      flags.logy = 0
	    endif else begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Linear Y Axis'
	      flags.logy = 1
            endelse

	  end

	  "STATION"				: begin

	    if checkstat eq 2 then begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Add Stations'
	      checkstat = 1
	    endif else begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Take off Stations'
	      checkstat = 2
            endelse

	  end

	  "MEAN"				: begin

	    if remmean eq 1 then begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Remove Mean'
	      remmean = 0
	    endif else begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Add Mean'
	      remmean = 1
            endelse

	  end

	  "PUBLISH"				: begin

	    if publish eq 1 then begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Take Date Stamp Off'
	      publish = 0
	    endif else begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Put Date Stamp On Plot'
	      publish = 1
            endelse

	  end

	  "CLR"					: begin

	    if clr eq 1 then begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Color Plots'
	      clr = 0
	    endif else begin
	      if n_elements(screen) eq 0 then			$
                widget_control, ev.id, set_value='Black and White Plots'
	      clr = 1
            endelse

	  end

	  "PERCENT"				: begin

	    percent = get_percent([0,100], percent)

	  end

	  "PERFONT"				: begin

	    fsize = change.fontsize
	    fsize = get_percent([0,300],fsize)
	    change.fontsize = fsize

	  end

	  "PERTITLE"				: begin

	    tsize = change.titlesize
	    tsize = get_percent([0,300],tsize)
	    change.titlesize = tsize

	  end

	  "WINSIZE"				: begin

	    tsize = float(change.y_win_size)/100.0
	    tsize = get_percent([0,1000],tsize)
	    change.y_win_size = fix(tsize*100)

	  end

	  "New Symbol" 	  			: new_symbol
	  "Change Position"			: change.place = 1

          else : begin
            
	    for i=0,nfile-1 do if value eq lfn(i) then n=i

	    if n_elements(n) gt 0 then begin
	      setup_file, fbase, lfn(n), titles, n
	      get_var, fbase, titles, n, numofvar
	      widget_control, lid(n), set_value = titles(n,*)
	      change.data = 1
	      change.place = 1
	    endif else print, 'Not an Option!'

	  end

	  endcase
      end

  "SLIDER": begin
              value = ev.value

	      for i=0,n_elements(sid)-1 do if ev.id eq sid(i) then n=i

              case (n) of 

	        17 : begin
		  if value gt 0 then unconpts = value 		$
		  else unconpts = 1000000.0
	        end

	        18 : columns = value

	        19 : offdiv = 1.0/float(value)

	        20 : offmul = float(value)

	        21 : begin

		       if value le 200 then 				      $
			 maxran = float(fix(17.0*(value-60)/140.0 + 65.5))
		       if value gt 200 then 				      $
			 maxran = float(fix(4.0*(value-200)/6.0 + 400.5))

                       widget_control, sid(26), 			      $
			 set_value=strcompress(string(maxran),/remove_all)

                     end

                22 : nb = value

                23 : begin

		  velocity = float(value/10.0)
		  vstring = strcompress(string(velocity),/remove_all)
		  vstring = strmid(vstring,0,5)
                  widget_control, sid(28),				$
		    set_value='Velocity of Event : '+vstring

		end		  

                24 : beam = value

                25 : step = value

                26 : print, 'Impossible'

                27 : inarow = value

		28 : print, 'Impossible'

		29 : change.minvel = value

	        else : begin

                  dumtime=time

	          if n ge 12 then begin

		    time(n) = value

	            dumtime(6) = time(0)+time(12)
		    dumtime(7) = time(1)
	            dumtime(8:11) = dumtime(2:5)+time(13:16)

		    if dumtime(11) gt 59 then begin
		      dumtime(11) = dumtime(11) - 60
		      dumtime(10) = dumtime(10) + 1
		    endif

		    if dumtime(10) gt 59 then begin
		      dumtime(10) = dumtime(10) - 60
		      dumtime(9) = dumtime(9) + 1
		    endif

		    if dumtime(9) gt 23 then begin
		      dumtime(9) = dumtime(9) - 24
		      dumtime(13)=dumtime(13)+1
		    endif

                    daysofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
		    daysofyear = 365
                    if ((dumtime(0) mod 4) eq 0) then begin
                      daysofmon(1) = daysofmon(1) + 1
		      daysofyear=daysofyear+1
		    endif

	 	    i = dumtime(1)
		    if n ne 13 then dumtime(8) = dumtime(2)+dumtime(13)

		    while (dumtime(8) gt daysofmon(i-1)) do begin

		      i=i+1
		      dumtime(8) = dumtime(8) - daysofmon(i-2)

		      if i eq 13 then begin

                        if ((dumtime(6) mod 4) eq 0) then begin
                          daysofmon(1) = daysofmon(1) - 1
		          daysofyear=daysofyear-1
		        endif

		        dumtime(6) = dumtime(6)+1

                        if ((dumtime(6) mod 4) eq 0) then begin
                          daysofmon(1) = daysofmon(1) + 1
		          daysofyear=daysofyear+1
		        endif

		        i=1

		      endif

		    endwhile

		    dumtime(7) = i

	          endif

                  if n lt 12 then begin

                    dumtime(n)=value

                    if (n le 2) then begin

                      daysofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
                      if ((dumtime(0) mod 4) eq 0) then                $
                        daysofmon(1) = daysofmon(1) + 1

                      if dumtime(2) gt daysofmon(dumtime(1)-1) then          $
                        dumtime(2) = daysofmon(dumtime(1)-1)

                    endif else begin

                      if (n ge 6) and (n le 8) then begin

                        daysofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
                        if ((dumtime(6) mod 4) eq 0) then              $
                          daysofmon(1) = daysofmon(1) + 1

                        if dumtime(8) gt daysofmon(dumtime(7)-1) then        $
                          dumtime(8) = daysofmon(dumtime(7)-1)

                      endif

                    endelse

	          endif

                  c_a_to_r, dumtime(0:5), btime
                  c_a_to_r, dumtime(6:11), etime
		  c_r_to_a, ftime, eftime
                  c_a_to_r, [ftime(0),ftime(1),ftime(2),23,59,59], endftime

                    if (btime gt etime) then begin

                      if (n ge 6) then dumtime(n) = time(n)         $
                      else begin
                        dumtime(6:11) = dumtime(0:5)
                        time(6:11) = time(0:5)
                        change.data = 1
                      endelse

                    endif else change.data = 1

                    if (etime gt endftime) then begin

                        dumtime(6:11) = [ftime(0)+1900,ftime(1),	$
					 ftime(2),23,59,59]
                        change.data = 1

                    endif

                  time(0:11) = dumtime(0:11)

	          c_a_to_r, time(0:5), btime
	          c_a_to_r, time(6:11), etime
	          dt = etime-btime

	          time(12) = fix(dt/(365.0*24.0*60.0*60.0))
	          time(13) = fix(dt/(24.0*60.0*60.0)) - 365*time(12)
	          time(14) = fix(dt/(60.0*60.0)) - 24*time(13) - 365*24*time(12)
	          time(15) = fix(dt/60.0) - 60*time(14) - 24*60*time(13) - $
			     365*24*60*time(12)
	          time(16) = fix(dt) - fix(dt/60.0)*60

                  loc = where(sid(0:16) ne 0, count)

                  for i=0,count-1 do          $
                    widget_control, sid(loc(i)), set_value=time(loc(i))

                endelse

              endcase

            end

  "TEXT_CH" : begin

	    widget_control, ev.id, get_value=value
	    if ev.id eq intext then titleofplot = value(0)
	    if ev.id eq inps then psfile = value(0)
	    if ev.id eq inmaxx then begin
	      if strlen(value(0)) gt 0 then flags.maxx = float(value(0))      $
	      else flags.maxx = -1.0e-30
	    endif
	    if ev.id eq inmaxy then begin
	      if strlen(value(0)) gt 0 then flags.maxy = float(value(0))      $
	      else flags.maxy = -1.0e-30
	    endif
	    if ev.id eq inminx then begin
	      if strlen(value(0)) gt 0 then flags.minx = float(value(0))      $
	      else flags.minx = -1.0e-30
	    endif
	    if ev.id eq inminy then begin
	      if strlen(value(0)) gt 0 then flags.miny = float(value(0))      $
	      else flags.miny = -1.0e-30
	    endif

	  end

  "LIST" : begin
	    for i=0,nfile-1 do if ev.id eq lid(i) then n=i 
	    if strmid(titles(n,ev.index),14,10) eq '(selected)' then begin
	      titles(n,ev.index) = strmid(titles(n,ev.index),0,14)
	      numofvar = numofvar - 1
	    endif else begin
	      titles(n,ev.index) = strmid(titles(n,ev.index),0,14)+	$
				   '(selected)'
	      numofvar = numofvar + 1
	    endelse
	    widget_control, ev.id, set_value = titles(n,*)
	    value=value(ev.index)
	    change.data = 1
	    change.place = 1
	  end
  else:
  endcase

  if n_elements(screen) eq 0 then begin
    WIDGET_CONTROL, ev.top, get_uvalue=out_text
    widget_control, out_text, set_value=string(value)
  endif

end

