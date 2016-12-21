function setup_symbol, plottype, screen

  if n_elements(screen) eq 0 then begin

    base = widget_base(title='Choose New Symbol', /column) 

    b = widget_button(base, value='Solid Lines', uvalue='Solid Lines')
    b = widget_button(base, value='Line With Crosses', 		$
	  uvalue='Line With Crosses')
    b = widget_button(base, value='Line With Stars', uvalue='Line With Stars')
    b = widget_button(base, value='Line With Diamonds', 		$
	  uvalue='Line With Diamonds')
    b = widget_button(base, value='Line With Triangles', 		$
	  uvalue='Line With Triangles')
    b = widget_button(base, value='Line With Squares', 		$
	  uvalue='Line With Squares')
    b = widget_button(base, value='Line With X', uvalue='Line With X')
    b = widget_button(base, value='Crosses', uvalue='Crosses')
    b = widget_button(base, value='Stars', uvalue='Stars')
    b = widget_button(base, value='Diamonds', uvalue='Diamonds')
    b = widget_button(base, value='Triangles', uvalue='Triangles')
    b = widget_button(base, value='Squares', uvalue='Squares')
    b = widget_button(base, value='Dots', uvalue='Dots')
    b = widget_button(base, value='Xs', uvalue='Xs')
    b = widget_button(base, value='Circles', uvalue='Circles')

    tp = 0
    info = {base : base, types : tp}

  endif else begin

    print, '1.  Solid Lines'
    print, '2.  Line With Crosses'
    print, '3.  Line With Stars'
    print, '4.  Line With Diamonds'
    print, '5.  Line With Triangles'
    print, '6.  Line With Squares'
    print, '7.  Line With X'
    print, '8.  Crosses'
    print, '9.  Stars'
    print, '10. Diamonds'
    print, '11. Triangles'
    print, '12. Squares'
    print, '13. Dots'
    print, '14. Xs'
    print, '15. Circles'

    tp = [0,-1,-2,-4,-5,-6,-7,1,2,4,5,6,3,7,8]
    info = {base : 0, types : tp}

  endelse

  return, info

end

