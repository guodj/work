
pro get_plot_type, screen = screen

  COMMON type, plottype, types

  if n_elements(screen) eq 0 then screen = 0 else screen = 1

  typeofplot  = ['One or More Standard X vs. T Graph(s)',  	$         ; 1  
                 'Multiple Variables on single Plots',  	$         ; 2
                 'Over Lapping Stack Plot (Offset Plots)',  	$         ; 3
                 'One or More Standard X vs. Y Graph(s)',  	$         ; 4
;                 'Polar Contour',  				$         ; 5
;                 'Cartesian Contour',  			$         ; 6
;                 'Image',  					$         ; 7
                 'Polar Vector Plots',  			$         ; 8
                 'Dwell Vector Plots',  			$         ; 9
                 'Magnetometer Vector Plots',  			$         ; 10
;                 'Polar Vector Plots',  			$         ; 11
;                 'Eight Polar Contours on single page',  	$         ; 12
                 'Vector Field Plots',  			$         ; 13
;                 'Eight Polar Line of Sight Velocity Plots',  	$         ; 14
;                 'Special Polar Contours for Sondrestrom',  	$         ; 15
;                 'DMSP Horizontal Velocities',  		$         ; 16
;                 'Density Plots',  				$         ; 17
;                 'Latitudinal Variations of Magnetometers',  	$         ; 18
;                 'DMSP Horizontal Velocities (Portrait)',  	$         ; 19
                 'Magnetometer Quick Veiws',            	$	  ; 20
                 'Magnetometer Clock Dial Plot',		$   	  ; 21
                 'X vs Y with all same X variable']            		  ; 22

  plot_arr = [	1, 						$
		2, 						$
		3,                                              $
		4,						$
;		5,						$
;		6,						$
;		7,						$
		8,						$
		9,						$
		10,						$
;		11,						$
;		12,						$
		13,						$
;		14,						$
;		15,						$
;		16,						$
;		17,						$
;		18,						$
;		19,						$
		20,						$
		21,						$
		22]	

  if screen eq 0 then begin

    firstb = widget_base(title='Choose a Plot Type')

    t38 = widget_base(firstb, /COLUMN)

    font = '-adobe-helvetica-medium-r-*-*-14-140*'

    b1 = WIDGET_LIST(t38, uvalue=typeofplot, value=typeofplot,		      $
	  ysize=n_elements(typeofplot), font=font)

    widget_control, firstb, /realize, /hourglass 
    XMANAGER, 'get_plottype', firstb, event_handler='gettypee'

  endif else begin

   print, ' '
   plottype = 0
   while (plottype le 0) or (plottype gt n_elements(typeofplot)) do begin
     for i=0,n_elements(typeofplot)-1 do			$
       print, tostr(i+1),'. ',typeofplot(i)
     print, ' '
     print, 'Choose plot type :'
     read, plottype
   endwhile

  endelse
                              
  plottype = plot_arr(plottype-1)

  types = typeofplot(where(plot_arr eq plottype))

  RETURN

END

