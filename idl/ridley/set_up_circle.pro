pro set_up_circle, num_per_page, num_total, pos, botstring,		$
	dshape, dswitch, dstring, dchars, ddate

  space = 0.08

  window, /free, xsize = 1100, ysize = 800.0

  pos_space, num_per_page, space, sizes

  plot, [0,1], /nodata,							$
	xstyle = 5,							$
	ystyle = 5,							$
	pos = [0,0,1,1]

  if num_total gt num_per_page then num_total=num_per_page

  pos = fltarr(num_total+2,4)

  psi = !pi*2.0*findgen(31)/30.0

  for i=0,num_total-1 do begin

    get_position, nb, space, sizes, i, dum_pos

    pos(i,*) = dum_pos

    if i gt 0 then begin

      plot, [0,1], /nodata,						$
	    xstyle = 5,							$
	    ystyle = 5,							$
	    pos = [0,0,1,1], /noerase

    endif

    dy = 0.03
    dx = 0.2*(dum_pos(2)-dum_pos(0))
    xs = dum_pos(0)+[0,0,1.0,1.0,0]*dx
    ys = dum_pos(3)+0.01+[0,1.0,1.0,0,0]*dy
    plots, xs, ys, linestyle = 2

    xstart = dum_pos(2) - dx
    xs = xstart+[0,0,1.0,1.0,0]*dx
    plots, xs, ys, linestyle = 2

    xstart = dum_pos(0) + dx
    dx = 0.6*dx/0.2
    xs = xstart+[0,0,1.0,1.0,0]*dx
    plots, xs, ys, linestyle = 2

    case dswitch(i) of

      0 : begin

	xyouts, xstart+dx/2.0, dum_pos(3)+0.02, '00:00:00 UT', /norm,	$
	  alignment = 0.5

      end

      1 : begin

	xyouts, xstart+dx/2.0, dum_pos(3)+0.02, '00:00 UT', /norm,	$
	  alignment = 0.5

      end

      2 : begin

	xyouts, xstart+dx/2.0, dum_pos(3)+0.02, dstring(i,1), /norm,	$
	  alignment = 0.5

      end

    endcase

    xyouts, dum_pos(0), dum_pos(3)+0.02, $
	dstring(i,0), /norm, alignment = 0.0

    xyouts, dum_pos(2), dum_pos(3)+0.02, $
	dstring(i,2), /norm, alignment = 1.0

    plot, [-1,1],[-1,1], /nodata,					$
	  xstyle = 5,							$
	  ystyle = 5,							$
	  pos = dum_pos, /noerase

    case dshape(i) of

      0 : begin

	xs = cos(psi)
	ys = sin(psi)

	plots, xs, ys, linestyle = 2

      end

      1 : begin

	xs = [-1.0,-1.0,1.0,1.0,-1.0]
	ys = [-1.0,1.0,1.0,-1.0,-1.0]

	plots, xs, ys, linestyle = 2

      end

      2 : begin

	xs = cos(psi)
	ys = sin(psi)

	plots, xs, ys, linestyle = 2
	plots, [-1.0,1.0],[0.0,0.0], linestyle=2
	plots, [0.0,0.0],[-1.0,1.0], linestyle=2

      end

    endcase

  endfor

  min = min([pos(0:num_total-1,1),0.05])

  plot, [0,1], /nodata,		     					$
	xstyle = 5,							$
	ystyle = 5,							$
	pos = [0,0,1,1], /noerase

  dy = min-0.01
  dx = 0.75
  xs = [0,0,1.0,1.0,0]*dx
  ys = [0,1.0,1.0,0,0]*dy
  plots, xs, ys, linestyle = 2

  pos(num_total,*) = [0.0,0.0,dx,min-0.01]

  xyouts, 0.5*dx, 0.01, botstring, alignment=0.5, /norm

  dx = 0.2
  xs = [0,0,1.0,1.0,0]*dx + (1.0-dx)
  plots, xs, ys, linestyle = 2
  if not ddate then 						$
    xyouts, 1.0-dx/2.0, 0.01, 'date on', alignment=0.5, /norm 	$
  else xyouts, 1.0-dx/2.0, 0.01, 'date off', alignment=0.5, /norm

  pos(num_total+1,*) = [0.8,0.0,0.8+dx,min-0.01]

end  

