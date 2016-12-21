
; -------------------------------------------------------------------
;  MAIN WIDGET
; -------------------------------------------------------------------


pro set_up_main, base, datafile

  common fileblock, starttime, endtime, varcount, ftime 
  common timeblock, sid, time, bftime, eftime
  common varblock, nfile, lid, titles, units, lfn
  COMMON type, plottype, types
  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common rep_data, change
  common data, single, t, dt, basetime, numofvar, numofpts
  common text_id, intext, inps, instring, inmaxx, inmaxy, inminx, inminy

    base = WIDGET_BASE(title=types(0), /row)
    b6= widget_base(base,/frame,/column)
    butbase = WIDGET_BASE(base,/frame,/column)
    b1 = WIDGET_BASE(base,/frame,/column)
    b3 = widget_base(base,/frame)
    b4 = widget_base(base,/frame)
    b5 = widget_base(base,/frame)
    b7 = widget_base(base,/frame)

;
; Left Most Widget:
;
;	Listing of the Variables to choose from
; type_8
; type_9
; type_13
;

    if plottype eq 8 or							      $
       plottype eq 13 or						      $
       plottype eq 9 then						      $
      ysize = fix(32.0/float(nfile) - nfile/3 + nfile/6) - 3		      $
    else								      $
      ysize = fix(34.0/float(nfile) - nfile/3 + nfile/6) - 3

    if nfile eq 1 then ysize = min([ysize, n_elements(titles(0,*))])

    lid = lonarr(nfile)
    lfn = strarr(nfile)

    for i=0,nfile-1 do begin

      j = 0
      if !version.os eq 'vms' then begin
        while (strmid(datafile(i),j,1) ne ']') and 			      $
	      (j lt strlen(datafile(i))) do j=j+1
      endif else begin
	lslash = where(byte(datafile(i)) eq 47, nslash)
	if nslash gt 0 then j = lslash(nslash-1)			$
	else j = strlen(datafile(i))
      endelse
      if j lt strlen(datafile(i)) then j=j+1 else j=0
      fn = strmid(filename(i),j,strlen(datafile(i))-j)
      lfn(i) = fn
      font = '-adobe-times-bold-r-*-*-18-180*'
      f1 = widget_button(b6, value=fn, uvalue = fn, font=font)
      font = '-adobe-times-medium-r-*-*-12-120*'
      f2 = widget_label(b6, value=starttime(i),font=font)
      f3 = widget_label(b6, value=endtime(i),font=font)

      t38 = widget_base(b6, /COLUMN)
      font = '-adobe-times-medium-r-*-*-14-140*'
      lid(i) = WIDGET_LIST(t38, ysize=ysize, uvalue=titles(i,0:varcount(i)),  $
			  value=titles(i,0:varcount(i)), xsize=15,font=font)

    endfor

;
;   type_8
;

    if plottype eq 8 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      f1 =widget_label(b6, value='Choose 5 Variables:',font=font)
      f2 =widget_label(b6, value='Range, Azimuth, Velocity,',font=font)
      f2 =widget_label(b6, value='Elevation, Event',font=font)

    endif

;
;   type_13
;

    if plottype eq 13 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      f1 =widget_label(b6, value='Choose 6 Variables:', font=font)
      f2 =widget_label(b6, value='Range, Azimuth, Radial Vel,',font=font)
      f2 =widget_label(b6, value='Angular Vel, Elevation, Event',font=font)

    endif

;
;   type_9
;

    if plottype eq 9 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      f1 =widget_label(b6, value='Choose 4 Variables:',font=font)
      f2 =widget_label(b6, value='Range, Azimuth, Velocity',font=font)
      f2 =widget_label(b6, value='Elevation',font=font)

    endif

;   type_8
;   type_9
;   type_13

    if (plottype ne 8) and 						      $
       (plottype ne 9) and 						      $
       (plottype ne 13) then begin

      font = '-adobe-times-medium-r-*-*-14-140*'
      t11 = widget_button(b6, value="Select All", 		$
		uvalue="SELECT",font=font, xsize=15)

    endif
;
; Second to Left Most Widget:
;
;	Main Controller  -  Contains PLOT and DONE
;

    t1 = WIDGET_TEXT(b1, ysize=1, value=[' '])
    WIDGET_CONTROL, base, set_uvalue=t1

    if plottype ne 0 then t11 = widget_base(b1, /row) else t11 = b1
    font = '-adobe-helvetica-medium-r-*-*-18-180*'
    t7 = widget_button(butbase, value="Done", uvalue = "DONE", font=font)
    font = '-adobe-helvetica-medium-r-*-*-14-140*'
    t7 = widget_button(butbase, value="Get New File", 		$
	uvalue = "New File", font=font)
    if plottype ne 0 then begin

      t7=widget_button(butbase,value="Get New Plot type",	$
	uvalue="New Plot Type", font=font)

      t3 = widget_label(butbase, value = ' - ')

      t3 = widget_label(b1, value = 'Title of Plot : ',font=font)
      font = '-adobe-times-medium-r-*-*-14-140*'
      intext = widget_text(b1, /editable, xsize=30, ysize=1, 	$
	/frame,font=font)

    endif else begin

      intext = 0

    endelse

;  type_1   Simple Plots
;  type_2   Multiple variables on single plot

    if plottype eq 1 or							      $
       plottype eq 2 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t40 = widget_button(b1, value = 'Change Range Selections',	$
		uvalue = 'CHANGE RANGE',font=font)

    endif

;  type_4   
;  type_22   

    if plottype eq 4 or							      $
       plottype eq 22 then begin

      t45 = widget_base(b1, /row)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t3 = widget_label(t45, value = 'Min X : ',font=font)
      font = '-adobe-times-medium-r-*-*-14-140*'
      inminx = widget_text(t45, /editable, xsize=7, ysize=1, 	$
	/frame, font=font)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t3 = widget_label(t45, value = 'Min Y : ',font=font)
      font = '-adobe-times-medium-r-*-*-14-140*'
      inminy = widget_text(t45, /editable, xsize=7, ysize=1, 	$
	/frame,font=font)

      t45 = widget_base(b1, /row)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t3 = widget_label(t45, value = 'Max X : ',font=font)
      font = '-adobe-times-medium-r-*-*-14-140*'
      inmaxx = widget_text(t45, /editable, xsize=7, ysize=1, 	$
	/frame,font=font)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t3 = widget_label(t45, value = 'Max Y : ',font=font)
      font = '-adobe-times-medium-r-*-*-14-140*'
      inmaxy = widget_text(t45, /editable, xsize=7, ysize=1, 	$
	/frame,font=font)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t11 = widget_button(butbase, value="Log X Axis", 	      $
	    uvalue="LOGX",font=font)

      t11 = widget_button(butbase, value="Log Y Axis", 	      $
	    uvalue="LOGY",font=font)
      
    endif else begin

      inminx = 0
      inminy = 0
      inmaxx = 0
      inmaxy = 0

    endelse


;  type_1   Simple Plots
;  type_3   Overlapping stack plots
;  type_4   
;  type_20  Magnetometer Quick Looks

    if plottype eq 1 or							      $
       plottype eq 3 or							      $
       plottype eq 4 or							      $
       plottype eq 20 or						      $
       plottype eq 22 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t11 = widget_button(b1, value='New Symbol', uvalue='New Symbol',	$
		font=font)

      a = findgen(31)*!pi/15.0
      usersym, 0.5*sin(a), 0.5*cos(a), /fill

    endif

;  type_8   Azimuth Scans
;  type_13

    if plottype eq 8 or							      $
       plottype eq 13 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t48 = widget_base(b1, /row)
      t11 = widget_button(butbase, value="Add Invariant Latitudes", 	$
		uvalue="INVAR", font=font)
      t13 = widget_button(butbase, value="Add Stations", 		$
		font=font, uvalue="STATION")

    endif

;  type_8   Azimuth Scans
;  type_9
;  type_13

    if (plottype eq 8) or						      $
       (plottype eq 9) or						      $
       (plottype eq 13) then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t46 = widget_base(b1, /row)
      t11 = widget_button(butbase, value="Black and White Plots",	$
	    uvalue="CLR", font=font)
      t47 = widget_base(b1, /row)

    endif

;  type_8   Azimuth Scans
;  type_9   Dwell Scans Scans
;  type_13

    if plottype eq 8 or							      $
       plottype eq 9 or							      $
       plottype eq 13 then begin

      maxran = 800.0

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t12 = widget_label(b1, value='Maximum Range', font=font)
      rl  = widget_label(b1, value=strcompress(string(maxran),	$
		/remove_all), font=font)
      mr  = widget_slider(b1, min=60, max=800, value=800, 	$
		/drag, /suppress, font=font)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t12 = widget_label(b1, value='Minimum Velocity', font=font)
      mv  = widget_slider(b1, min=0, max=1000, value=change.minvel, 	$
		font=font)

    endif else begin

      rl = long(0)
      mr = long(0)
      mv = long(0)

    endelse

;  type_8   Azimuth Scans
;  type_13 
;  type_21  Magnetometer Clock Dial

    if plottype eq 8 or							      $
       plottype eq 13 or						      $
       plottype eq 21 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t12 = widget_label(b1, value='Number of Plots Per Page',	font=font)
      pp  = widget_slider(b1, min=1, max=20, value=4, font=font)

    endif else begin

      pp = long(0)

    endelse

;  type_1   Simple Plots
;  type_2   Multiple variables on single plot
;  type_3   Overlapping stack plots
;  type_4
;  type_20  Magnetometer Quick Looks

    if plottype eq 1 or							      $
       plottype eq 2 or							      $
       plottype eq 3 or							      $
       plottype eq 4 or							      $
       plottype eq 20 then begin	

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t12 = widget_label(b1, value='Minimum Seconds for Data Gap', font=font)
      dg  = widget_slider(b1, min=0, max=10000, value=0, font=font)

    endif else begin

      dg = long(0)

    endelse

;  type_3   Overlapping stack plots
;  type_20  Magnetometer Quick Looks

    if plottype eq 3 or							      $
       plottype eq 20 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t11 = widget_button(butbase, value="Remove Mean", 	$
		uvalue="MEAN", font=font)

    endif

;  type_3   Overlapping stack plots
;  type_8   Azimuth Scans
;  type_9   Dwell Plots
;  type_10  Magnetometer Velocity
;  type_13
;  type_20  Magnetometer Quick Looks
;  type_21  Magnetometer Clock Dial

    if plottype eq 3 or							      $
       plottype eq 8 or							      $
       plottype eq 9 or							      $
       plottype eq 10 or						      $
       plottype eq 13 or						      $
       plottype eq 20 or						      $
       plottype eq 21 then begin	

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t12 = widget_label(b1, value='Scale Divider', font=font)
      sd  = widget_slider(b1, min=1, max=20, value=1, font=font)

      t12 = widget_label(b1, value='Scale Multiplier', font=font)
      sm  = widget_slider(b1, min=1, max=20, value=1, font=font)

    endif else begin

      sd = long(0)
      sm = long(0)

    endelse

;  type_21  Magnetometer Clock Dial

    if plottype eq 21 then begin	

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t12 = widget_label(b1, value='Number of Points to Skip',font=font)
      po  = widget_slider(b1, min=0, max=120, value=0,font=font)

    endif else begin

      po = long(0)

    endelse

;  type_10  Magnetometer Velocity

    if plottype eq 10 then begin	

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      vl = widget_label(b1, value='Velocity of Event : 0.0 Km/s', font=font)
      ve  = widget_slider(b1, min=-150, max=150, value=0,	$
		/drag, /suppress, font=font)

      t12 = widget_label(b1, value='Beam angle for Dwell', font=font)
      dw  = widget_slider(b1, min=0, max=361, value=361,font=font)

    endif else begin

      ve = long(0)
      vl = long(0)
      dw = long(0)

    endelse

    if plottype eq 1 or							      $
       plottype eq 2 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t12 = widget_label(b1, value='Number of Columns',font=font)
      nc  = widget_slider(b1, min=1, max=5, value=1, font=font)

    endif else begin

      nc = long(0)

    endelse

;  idlflat
;  type_10  Magnetometer Velocity
;  type_20  Magnetometer Quick Looks
;  type_21  Magnetometer Clock Dial

    if plottype ne 0 and						      $
       plottype ne 10 and						      $
       plottype ne 20 and						      $
       plottype ne 21 then begin

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t8=widget_button(butbase, value="Change Position", 	$
		uvalue="Change Position",font=font)

    endif

;
; PS file stuff - for All plots!!
;

    if plottype ne 0 then begin

      t48 = widget_base(b1, /row)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t11 = widget_button(butbase, value="Extra! Extra!", 	$
		uvalue="EXTRAS",font=font)

      t3 = widget_label(butbase, value = ' - ',font=font)

      t3 = widget_label(b1, value = 'PS File name : ',font=font)
      font = '-adobe-times-medium-r-*-*-14-140*'
      inps = widget_text(b1, /editable, xsize=30, ysize=1, 	$
		/frame,font=font)

      font = '-adobe-helvetica-medium-r-*-*-18-180*'
      t8 = widget_button(butbase, value="Plot", uvalue = "PLOT",font=font)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'
      t8 = widget_button(butbase, value="Animate", uvalue = "ANIMATE",	$
		font=font)

      font = '-adobe-helvetica-bold-r-*-*-14-140*'
      psbase = widget_base(b1, /row)
      t8 = widget_button(psbase, value = "Make PS (land)",		$
   		uvalue = "Make PS File (landscape)",font=font,		$
		xsize=130)

      t8 = widget_button(psbase, value = "Make PS (port)",		$
		uvalue = "Make PS File (portrait)",font=font,		$
		xsize=130)

      font = '-adobe-helvetica-medium-r-*-*-12-120*'

      epsbase = widget_base(b1, /row)
      t8 = widget_button(epsbase, value = "Make EPS (land)",		$
   		uvalue = "Make EPS File (landscape)",font=font,		$
		xsize = 130)                                  

      t8 = widget_button(epsbase, value = "Make EPS (port)",		$
		uvalue = "Make EPS File (portrait)",font=font,		$
		xsize = 130)

      font = '-adobe-helvetica-medium-r-*-*-14-140*'

      if plottype eq 1 then 						$
        t8 = widget_button(butbase, value="Write Ascii to File",	$
	   uvalue = "SFLD")

      if !version.os eq 'vms' then begin

        t8 = widget_button(butbase, value = "print b/w (land)",	      $
   		uvalue = "Print plot to oosik (landscape)",font=font)

        t8 = widget_button(butbase, value = "print b/w (port)",		      $
		uvalue = "Print plot to oosik (Portrait)",font=font)
      endif

      t9 = widget_label(butbase, value = ' - ')

      t11 = widget_button(butbase, value='Take Date Stamp Off', 	$
	uvalue="PUBLISH",font=font)

      t11 = widget_button(butbase, value='Resize Postscript Image', 	$
	uvalue="PERCENT",font=font)

      t11 = widget_button(butbase, value='Resize Postscript Fonts', 	$
	uvalue="PERFONT",font=font)

      t11 = widget_button(butbase, value='Resize Postscript Title', 	$
	uvalue="PERTITLE",font=font)

      t11 = widget_button(butbase, value='Resize Window Screen', 	$
	uvalue="WINSIZE",font=font)

    endif else begin

      t8 = widget_button(butbase, value="Write Ascii to Screen or File", 	      $
	   uvalue = "SFLD")
      t8 = widget_button(butbase, value="Write new Flatfile from Variables",       $
	   uvalue = "SELF")
      t8 = widget_button(butbase, value="FFT Selected Variables",       $
	   uvalue = "FFT")

      t3 = widget_label(b1, value = 'New File name : ')
      inps = widget_text(b1, /editable, xsize=30, ysize=1, /frame)

    endelse

;
; Left Most Two Time Widgets:
;
;	Contain Starting Time and Ending Time 
;

    font = '-b&h-lucida-medium-r-normal-sans-14-100*'

    t10 = widget_label(b3, value='START', yoff=05, font=font)
    t11 = widget_label(b3, value='TIME', yoff=25, font=font)

    t24 = widget_label(b4, value='END', yoff=05, font=font)
    t25 = widget_label(b4, value='TIME', yoff=25, font=font)

    if ftime(0) eq ftime(6) then begin

      t12 = widget_label(b3, value='year :', yoff=50, font=font)
      t12 = widget_label(b3, value=strcompress(string(ftime(0))),    $
			  yoff=70, font=font)

      s1  = 0

      t12 = widget_label(b4, value='year :', yoff=50, font=font)
      t12 = widget_label(b4, value=strcompress(string(ftime(0))),    $
			  yoff=70, font=font)

      e1 = 0

      if ftime(1) eq ftime(7) then begin

        t14 = widget_label(b3, value='mon :', yoff=135, font=font)
	t14 = widget_label(b3, value=strcompress(string(ftime(1))), $
			    yoff=155, font=font)

        s2 = 0

        t14 = widget_label(b4, value='mon :', yoff=135, font=font)
	t14 = widget_label(b4, value=strcompress(string(ftime(1))), $
			    yoff=155, font=font)

        e2 = 0


        if ftime(2) eq ftime(8) then begin

          t14 = widget_label(b3, value='day :', yoff = 255, font=font)
	  t14 = widget_label(b3, value=strcompress(string(ftime(2))), $
			      yoff=275, font=font)

          s3 = 0

          t14 = widget_label(b4, value='day :', yoff = 255, font=font)
	  t14 = widget_label(b4, value=strcompress(string(ftime(2))), $
			      yoff=275, font=font)

          e3 = 0

        endif else begin

          t16 = widget_label(b3, value='day', yoff=255, font=font)
          s3 = widget_slider(b3, min=ftime(2), max=ftime(8), value=ftime(2), $
			/vertical, ysize=100, yoff=280, xoff=15, 	$
			font=font)

          t28 = widget_label(b4, value='day', yoff=255, font=font)
          e3 = widget_slider(b4, min=ftime(2), max=ftime(8), value=ftime(8), $
			/vertical, ysize=100, yoff=280, xoff=15, font=font)

        endelse

      endif else begin

        t14 = widget_label(b3, value='mon', yoff=135, font=font)
        s2 = widget_slider(b3, min=ftime(1), max=ftime(7), value=ftime(1), $
		/vertical, ysize=50, yoff=160, xoff=15, font=font)

        t14 = widget_label(b4, value='mon', yoff=135, font=font)
        e2 = widget_slider(b4, min=ftime(1), max=ftime(7), value=ftime(7), $
		/vertical, ysize=50, yoff= 160, xoff=15, font=font)

        t16 = widget_label(b3, value='day', yoff=255, font=font)
        s3 = widget_slider(b3, min=1, max=31, value=ftime(2),		$
		/vertical, ysize=100, yoff=280, xoff=15, font=font)

        t28 = widget_label(b4, value='day', yoff=255, font=font)
        e3 = widget_slider(b4, min=1, max=31, value=ftime(8), 	$
		/vertical, ysize=100, yoff=280, xoff=15, font=font)

      endelse

    endif else begin

      t12 = widget_label(b3, value='year', yoff=50, font=font)
      s1  = widget_slider(b3, min=ftime(0), max=ftime(6), value=ftime(0), $
		/vertical, ysize=50, yoff=75, font=font)

      t26 = widget_label(b4, value='year', yoff=50, font=font)
      e1  = widget_slider(b4, min=ftime(0), max=ftime(6), value=ftime(6), $
		/vertical, ysize=50, yoff=75, font=font)

      t14 = widget_label(b3, value='mon', yoff=135, font=font)
      s2 = widget_slider(b3, min=1, max=12, value=ftime(1), /vertical, 	$
		ysize=50, yoff=160, xoff=15, font=font)

      t28 = widget_label(b4, value='mon', yoff=135, font=font)
      e2 = widget_slider(b4, min=1, max=12, value=ftime(7), /vertical, 	$
		ysize=50, yoff=160, xoff=15, font=font)

      t16 = widget_label(b3, value='day', yoff=255, font=font)
      s3 = widget_slider(b3, min=1, max=31, value=ftime(2), /vertical, 	$
		ysize=100, yoff=280, xoff=15, font=font)

      t28 = widget_label(b4, value='day', yoff=255, font=font)
      e3 = widget_slider(b4, min=1, max=31, value=ftime(8), /vertical, 	$
		ysize=100, yoff=280, xoff=15, font=font)

    endelse

    t18 = widget_label(b3, value='hour', yoff=390, font=font)
    s4 = widget_slider(b3, min=0, max=23, value=0, /vertical, ysize=65,	      $
		yoff=415, xoff=15, font=font)
    t20 = widget_label(b3, value='min', yoff=490, font=font)
    s5 = widget_slider(b3, min=0, max=59, value=0, /vertical, ysize=100,      $
		yoff=515, xoff=15, font=font)
    t22 = widget_label(b3, value='sec', yoff=625, font=font)
    s6 = widget_slider(b3, min=0, max=59, value=0, /vertical, ysize=100,      $
		yoff=650, xoff=15, font=font)

    t32 = widget_label(b4, value='hour', yoff=390, font=font)
    e4 = widget_slider(b4, min=0, max=23, value=0, /vertical, ysize=65,      $
		yoff=415, xoff=15, font=font)
    t34 = widget_label(b4, value='min', yoff = 490, font=font)
    e5 = widget_slider(b4, min=0, max=59, value=0, /vertical, ysize=100,     $
		yoff=515, xoff=15, font=font)
    t36 = widget_label(b4, value='sec', yoff=625, font=font)
    e6 = widget_slider(b4, min=0, max=59, value=0, /vertical, ysize=100,     $
		yoff=650, xoff=15, font=font)

;
; Right Most Time Widget:
;
;	Contains Plot Length
;

    t10 = widget_label(b5, value='Plot', yoff=05, font=font)
    t11 = widget_label(b5, value='Length', yoff=25, font=font)

    nyear = fix(dt/(365.0*24.0*60.0*60.0))

    ndays = fix(dt/(24.0*60.0*60.0)) + 1
    nadays = ndays - nyear*365
    if ndays gt 365 then ndays = 365

    if nyear ge 1 then begin

      t12 = widget_label(b5, value='years', yoff=50, font=font)
      l1  = widget_slider(b5, min=0, max=nyear+1, value=nyear, /vertical,  $
		ysize=50, yoff=75, xoff=21, font=font)

    endif else l1 = 0

    if ndays ge 1 then begin

      t12 = widget_label(b5, value='days', yoff=135, font=font)
      pix = 220
      if ndays ge 100 then xoff=8 else xoff=15
      l2  = widget_slider(b5, min=0, max=ndays, value=nadays, 	$
		/vertical, ysize=pix, yoff=160, xoff=xoff, font=font)

    endif else l2 = 0

    t12 = widget_label(b5, value='hours', yoff=390, font=font)
    l3  = widget_slider(b5, min=0, max=24, value=0, /vertical, 	$
		ysize=65, yoff=415, xoff=15, font=font)
    t12 = widget_label(b5, value='mins', yoff=490, font=font)
    l4  = widget_slider(b5, min=0, max=60, value=0, /vertical, 	$
		ysize=100, yoff=515, xoff=15, font=font)
    t12 = widget_label(b5, value='secs', yoff=625, font=font)
    l5  = widget_slider(b5, min=0, max=60, value=0, /vertical, 	$
		ysize=100, yoff=650, xoff=15, font=font)

;
; Right Most Time Widget:
;
;	Contains Plot Length
;

    t10 = widget_label(b7, value='Number', yoff=05, font=font)
    t11 = widget_label(b7, value='of Plots', yoff=25, font=font)
    t11 = widget_label(b7, value='in a row', yoff=45, font=font)
    t11 = widget_label(b7, value='using Plot', yoff=65, font=font)
    t11 = widget_label(b7, value='Length', yoff=85, font=font)

    np  = widget_slider(b7, min=1, max=50, value=1, 		$
		ysize=100, yoff=125, xoff=15, /vertical, font=font)

    sid=[s1,s2,s3,s4,s5,s6,e1,e2,e3,e4,e5,e6,l1,l2,l3,l4,l5,		      $
	 dg,nc,sd,sm,mr,pp,ve,dw,po,rl,np,vl,mv]

  RETURN

END


