
; -------------------------------------------------------------------
;  MAIN WIDGET
; -------------------------------------------------------------------


pro d_menu_bat

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

  print, ' '
  print, 'VAR      -  Choose variable to plot'
  

;   type_8
;   type_9
;   type_13

    if (plottype ne 8) and 						      $
       (plottype ne 9) and 						      $
       (plottype ne 13) then begin

      print, 'SELECT   -  Select/Deselect All Variables for plotting'

    endif

   print, 'DISPLAY  -  Redisplay variables'
   print, 'HELP, ?  -  Redisplay this menu'
   print, 'QUIT     -  quits program (EXIT or Q also work)'
   print, ' '

;  type_1   Simple Plots
;  type_2   Multiple variables on single plot

    if plottype eq 1 or							      $
       plottype eq 2 then begin

      print, 'RANGE    -  Change ranges/scales'

    endif

;  type_4   
;  type_22   

    if plottype eq 4 or							      $
       plottype eq 22 then begin

      print, 'XYRAN    -  Change max and min on X and Y range'
      print, 'LOGX     -  Switch to change X axis to Log/Regular'
      print, 'LOGY     -  Switch to change Y axis to Log/Regular'

    endif

;  type_1   Simple Plots
;  type_3   Overlapping stack plots
;  type_4   
;  type_20  Magnetometer Quick Looks

    if plottype eq 1 or							      $
       plottype eq 3 or							      $
       plottype eq 4 or							      $
       plottype eq 20 or						      $
       plottype eq 22 then begin

      print, 'SYMBOL   -  Change symbols for plot

      a = findgen(31)*!pi/15.0
      usersym, 0.5*sin(a), 0.5*cos(a), /fill

    endif

;  type_8   Azimuth Scans
;  type_13

    if plottype eq 8 or							      $
       plottype eq 13 then begin

      print, 'INVAR    -  Add Invariant Latitudes to plots'
      print, 'STATION  -  Add Stations'

    endif

;  type_8   Azimuth Scans
;  type_9
;  type_13

    if (plottype eq 8) or						      $
       (plottype eq 9) or						      $
       (plottype eq 13) then begin

      print, 'CLR      -  Switch between black and white plots / color plots'

    endif

;  type_8   Azimuth Scans
;  type_9   Dwell Scans Scans
;  type_13

    if plottype eq 8 or							      $
       plottype eq 9 or							      $
       plottype eq 13 then begin

      maxran = 800.0

      print, 'MAXR     -  Set maximum range to plot' 

    endif

    rl = long(0)
    mr = long(0)

;  type_8   Azimuth Scans
;  type_13 
;  type_21  Magnetometer Clock Dial

    if plottype eq 8 or							      $
       plottype eq 13 or						      $
       plottype eq 21 then begin

      print, 'PPP      -  Set the number of plots per page'

    endif

      pp = long(0)

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

      print, 'DG       -  Set Minimum Seconds for Data Gap'

    endif

      dg = long(0)

;  type_3   Overlapping stack plots
;  type_20  Magnetometer Quick Looks

    if plottype eq 3 or							      $
       plottype eq 20 then begin

      print, 'MEAN     -  Remove Mean'

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

      print, 'SM       -  Set Scale multiplication factor'

    endif

      sd = long(0)
      sm = long(0)

;  type_21  Magnetometer Clock Dial

    if plottype eq 21 then begin	

      print, 'PTS      -  Set time spacing between plots'

    endif

      po = long(0)

;  type_10  Magnetometer Velocity

    if plottype eq 10 then begin	

      print, 'VEL      -  Set the phase velocity of the TCV'
      print, 'BEAM     -  Set the beam angle for imitation dwells'

    endif

      ve = long(0)
      dw = long(0)

    if plottype eq 1 or							      $
       plottype eq 2 then begin

      print, 'COL      -  Set number of Columns'

    endif

      nc = long(0)

;
; PS file stuff - for All plots!!
;

    if plottype ne 0 then begin

      print, 'TITLE    -  Enter/change title of plot'
      print, 'ROW      -  Select number of time periods in a row to plot'
      print, 'PERCENT  -  Resize PS image'
      print, 'PERFONT  -  Resize Font size'
      print, 'PERTITLE -  Resize Font for title'
      print, 'PUBLISH  -  Switch to take date stamp off'
      print, 'EXTRA    -  Labels, lines, fonts, curves, ect.'
      print, ' '
      print, 'PLOT     -  Plots specified variables in specified time'
      print, 'PSNAME   -  Sets PS File name 
      print, 'PSL      -  Makes a PS file (Landscape)'
      print, 'PSP      -  Makes a PS file (Portrait)'
      print, 'EPSL     -  Makes a EPS file (Landscape)'
      print, 'EPSP     -  Makes a EPS file (Portrait)'
      print, 'OL       -  Print straight to oosik (Landscape)'
      print, 'OP       -  Print straight to oosik (Portrait)'
      print, ' '
      print, 'TIME     -  Select start time and end time'
      print, ' '
      if plottype eq 1 then begin
        print, 'SFLD     -  Write selected columns to ASCII file'
        print, ' '
      endif
      print, 'DONE     -  Finished with plotting routine'
      print, 'FILE     -  Get New Files for plotting'
      print, 'PTYPE    -  Get New Plot type'

    endif

  RETURN

END


