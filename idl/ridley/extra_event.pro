pro extra_event, ev

  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common comdis, displays
  common text_id, intext, inps, instring, inmaxx, inmaxy, inminx, inminy

  widget_control, ev.id, get_uvalue=value
  if (n_elements(value) eq 0) then value = ''
  name = strmid(tag_names(ev, /structure_name), 7, 1000)

  case (name) of

  "BUTTON": begin

        case (value) of

	  "DCHANGE" : begin

	    temp1 = displays.dshape
	    temp2 = displays.dswitch
	    temp3 = displays.dstring
	    temp4 = displays.dchars
	    temp5 = displays.ddate
	    change_circle, nb, nb, temp1, temp2, temp3, temp4, temp5
	    displays.dshape = temp1
	    displays.dswitch = temp2
	    displays.dstring = temp3
	    displays.dchars = temp4
	    displays.ddate = temp5

	  end

	  "CRB"					: begin

	    if crb eq 1 then begin
              widget_control, ev.id, set_value='Put Reversal Line On Plot'
	      crb = 0
	    endif else begin
              widget_control, ev.id, set_value='Take Reversal Line Off'
	      crb = 1
            endelse

	  end

	  "DRAW"				: begin

            if !d.window ne -1 then begin

	      draw_curve, final_x, final_y, ext.opos
	      draw = 1

	    endif

	  end

	  "IN_CURVE"				: begin

	      io_curve, final_x, final_y, /in
	      draw = 1

	  end

	  "OUT_CURVE"				: begin

	      io_curve, final_x, final_y, /out
	      draw = 1

	  end

	  "IN_LINE"				: begin

	      io_curve, l_x, l_y, /in
	      ext.line_x(0:n_elements(l_x)-1) = l_x
	      ext.line_y(0:n_elements(l_y)-1) = l_y
	      draw = 1

	  end

	  "OUT_LINE"				: begin

	      l_x = ext.line_x
	      l_y = ext.line_y
	      io_curve, l_x, l_y, /out
	      draw = 1

	  end

	  "CLEAR_ALL"				: begin

	    final_x = [-1.0]
	    final_y = [-1.0]
	    draw = 0

	  end

	  "CLEAR_1"				: begin

	    final_x = [-1.0]
	    final_y = [-1.0]
	    draw = 0

	  end

	  "LA_CLEAR_1"				: begin

	    n_e = n_elements(str_s)
	    if n_e ge 2 then begin
	      str_x = str_x(0:n_e-2)
	      str_y = str_y(0:n_e-2)
	      str_s = str_s(0:n_e-2)
	      str_o = str_o(0:n_e-2)
	      str_l = str_l(0:n_e-2)
	    endif else begin
	      str_x = [-1.0]
	      str_y = [-1.0]
	      str_s = [-1.0]
	      str_o = [-1.0]
	      str_l = ['']
            endelse

	  end

	  "LA_CLEAR_ALL"			: begin

	    str_x = [-1.0]
	    str_y = [-1.0]
	    str_s = [-1.0]
	    str_o = [-1.0]
	    str_l = ['']

	  end

	  "ARROW"				: begin

            if !d.window ne -1 then begin

	      draw_arrow, arrow_x, arrow_y, ext.opos
	      arrow = 1

	    endif

	  end

	  "AR_CLEAR_1"				: begin

	    loc = where(arrow_x eq -1, count)
	    if count ge 2 then begin
	      arrow_x = arrow_x(0:loc(count-2))
	      arrow_y = arrow_y(0:loc(count-2))
	    endif else begin
	      arrow_x = [-1.0]
	      arrow_y = [-1.0]
	      arrow = 0
            endelse

	  end

	  "AR_CLEAR_ALL"			: begin

	    arrow_x = [-1.0]
	    arrow_y = [-1.0]
	    arrow = 0

	  end

	  "LINE"				: begin

            if !d.window ne -1 then begin

	      loc = where(ext.line_x eq -1, count)
	      if count ge 1 then begin
	        l_x = ext.line_x(0:loc(count-1))
	        l_y = ext.line_y(0:loc(count-1))
	      endif else begin
		l_x = [-1.0]
		l_y = [-1.0]
	      endelse
	      draw_line, l_x, l_y, ext.opos
	      loc = where(l_x eq -1, count)
              ext.line_x(0:loc(count-1)) = l_x(0:loc(count-1))
              ext.line_y(0:loc(count-1)) = l_y(0:loc(count-1))

	    endif

	  end

	  "LI_CLEAR_1"				: begin

	    loc = where(ext.line_x eq -1, count)
	    n_e = n_elements(ext.line_x)-1
	    if count ge 2 then begin
	      ext.line_x(loc(count-2)+1:n_e-1) = 0.0
	      ext.line_y(loc(count-2)+1:n_e-1) = 0.0
	    endif

	  end

	  "LI_CLEAR_ALL"			: begin

	    n_e = n_elements(ext.line_x)-1
            ext.line_x(1:n_e) = 0.0
            ext.line_y(1:n_e) = 0.0
            ext.line_x(0) = -1.0
            ext.line_y(0) = -1.0

	  end

	  "Courier" 				: psfont = 0
	  "Courier, Bold" 			: psfont = 1
	  "Courier, Oblique" 			: psfont = 2
	  "Courier, Bold, Oblique"		: psfont = 3
	  "Helvetica"				: psfont = 4
	  "Helvetica, Bold"			: psfont = 5
	  "Helvetica, Oblique"			: psfont = 6
	  "Helvetica, Bold, Oblique"		: psfont = 8
	  "Avantgarde, Book"			: psfont = 12
	  "Avantgarde, Book, Oblique"		: psfont = 13
	  "Avantgarde, Demi"			: psfont = 14
	  "Avantgarde, Demi, Oblique"		: psfont = 15
	  "Schoolbook"				: psfont = 20
	  "Schoolbook, Bold"			: psfont = 21
	  "Schoolbook, Italic"			: psfont = 22
          "Schoolbook, Bold, Italic"		: psfont = 23
	  "Times" 				: psfont = 28
	  "Times, Bold" 			: psfont = 29
	  "Times, Italic"                       : psfont = 30
	  "Times, Bold, Italic"                 : psfont = 31

          else : 
	endcase

      end

    "TEXT_CH" : begin

	    widget_control, ev.id, get_value=value
	    if ev.id eq instring then begin
	      str = value(0)
	      str_get, str_x, str_y, str_s, str_o,str_l, str, ext.opos
	    endif 

    end

  else:
  endcase

  widget_control, ev.top, /destroy
  return

end

