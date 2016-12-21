pro change_circle, num, total_plots, dshape, dswitch, dstring, 	$
		   dchars, ddate

;  dshape = intarr(total_plots)
;  dswitch = intarr(total_plots)
;  dstring = strarr(total_plots,3)
;  ddate = 0

  psi = !pi*2.0*findgen(31)/30.0

;  x and y positions of a circle :

  xc = cos(psi)
  yc = sin(psi)

;  x and y positions of a square :

  xs = [-1.0,-1.0,1.0,1.0,-1.0]
  ys = [-1.0,1.0,1.0,-1.0,-1.0]

;  determine how many pages we have :
;    (although the number of pages option doesn't work too well)

  pages = total_plots/num
  if total_plots mod num gt 0 then pages = pages + 1

  for i = 0, pages - 1 do begin

;  determine what page we are on and what we should display at the bottom of
;  the page:

    if i eq pages-1 then begin

      plots = total_plots mod num 
      if plots eq 0 then plots = num
      botstring = 'Click in this Box to End'

    endif else begin

      plots=num
      botstring = 'Click in this Box to move to page Number '+		$
	strcompress(string(i+1),/remove_all)

    endelse 

;  Set up the plot area for the plots :

    set_up_circle, num, plots, pos, botstring, 				$
	dshape, dswitch, dstring, dchars, ddate

    out = 0

    while not out do begin

      get_circlen, pos, plotn, eventn

      if (plotn ge plots) then begin
	case (plotn) of 
	  plots : out = 1
	  plots + 1 : begin
	    if ddate then begin
	      xyouts, 0.9, 0.01, 'date off', alignment=0.5, /norm,	$
		color = 0
	      ddate = 0
	      xyouts, 0.9, 0.01, 'date on', alignment=0.5, /norm
	    endif else begin
	      xyouts, 0.9, 0.01, 'date on', alignment=0.5, /norm,	$
		color = 0
	      ddate = 1
	      xyouts, 0.9, 0.01, 'date off', alignment=0.5, /norm
	    endelse
	  end
	endcase 
      endif else begin

        if eventn eq 0 then begin

	  plot, [-1,1], [-1,1], /nodata, /noerase,			$
	    xstyle=5, ystyle=5, pos = pos(plotn,*)

          if dshape(i*num+plotn) eq 0 then begin
	    plots, xc, yc, color = 0
	    plots, xs, ys, linestyle = 2
	    dshape(i*num+plotn) = 1
	  endif else begin

            if dshape(i*num+plotn) eq 1 then begin
	      plots, xs, ys, color = 0
	      plots, xc, yc, linestyle = 2
	      dshape(i*num+plotn) = 0
	    endif

	  endelse

        endif

        if (eventn eq 1) or (eventn eq 3) or				$
	   ((eventn eq 2) and (dswitch(i*num+plotn) eq 1)) then begin

	  print, ''
	  print, 'Please Enter Text (Press Return to Finish)'
	  print, 'Back space does not work, so if you make a mistake,'
	  print, 'Press the back arrow key and it will be the back space.'

	  if eventn eq 2 then begin
	    xyouts,xyoutx,xyouty,'00:00 UT', color = 0, /norm, alignment=0.5
	    dswitch(i*num+plotn) = 2
          endif

	  dum$ = ''

	  by = byte('')
	  a = ''

	  repeat begin

	    if by(0) eq 27 then begin
	      a = fix(byte(get_kbrd(1)))
	      b = fix(byte(get_kbrd(1)))
              if (a(0) eq 91) and (b(0) eq 68) and 			$
	        (strlen(dum$) gt 0) then				$
	        dum$ = strmid(dum$,0,strlen(dum$)-1)
	    endif

	    if (by(0) ne 27) then dum$ = dum$ + a

	    xyouty = pos(plotn,3)+0.02

            if eventn eq 1 then begin

	      xyoutx = pos(plotn,0)
	      alignment = 0.0

	    endif

            if eventn eq 2 then begin

	      xyoutx = (pos(plotn,2)+pos(plotn,0))/2.0
	      alignment = 0.5
	      dswitch(i*num+plotn) = 2

	    endif

            if eventn eq 3 then begin

	      xyoutx = pos(plotn,2)
	      alignment = 1.0

	    endif

	    xyouts, xyoutx, xyouty,dstring((i*num+plotn),eventn-1), 	$
	      alignment=alignment, /norm, color = 0

	    dstring((i*num+plotn),eventn-1) = dum$

	    xyouts, xyoutx, xyouty,dstring((i*num+plotn),eventn-1), 	$
	      alignment=alignment, /norm

	    a = get_kbrd(1)

	    by = fix(byte(a))

	  endrep until by(0) eq 13

	  print, 'Done entering the string...'

	  dstring((i*num+plotn),eventn-1) = dum$

        endif else begin

          if eventn eq 2 then begin

	    xyoutx = (pos(plotn,2)+pos(plotn,0))/2.0
	    xyouty = pos(plotn,3)+0.02
  
	    case dswitch(i*num+plotn) of

              0 : begin
	        xyouts,xyoutx,xyouty,'00:00:00 UT',color = 0, 		$
		  /norm, alignment=0.5
	        xyouts,xyoutx,xyouty,'00:00 UT', /norm, alignment=0.5
	        dswitch(i*num+plotn) = 1
	      end

              2 : begin
	        xyouts,xyoutx,xyouty,dstring((i*num+plotn),1), /norm, 	$
		  color = 0, alignment=0.5
	        xyouts,xyoutx,xyouty,'00:00:00 UT', /norm, alignment=0.5
	        dswitch(i*num+plotn) = 0
	      end

	    endcase

          endif

        endelse

      endelse

    endwhile

    wdelete

  endfor

  return

end

