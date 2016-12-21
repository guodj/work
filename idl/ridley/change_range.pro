pro change_range

  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common range_type, mmq, cr
  common rep_data, change

  base = widget_base(title = 'Choose Range Style :', /column)

  t1 = widget_button(base, value = 'All Plots Same Range', 	$
	uvalue = 'APSR')
  t1 = widget_button(base, value = 'Some Plots Same Range', 	$
	uvalue = 'SPSR')
  t1 = widget_button(base, value = 'All Plots Same Scale', 	$
	uvalue = 'APSS')
  t1 = widget_button(base, value = 'Some Plots Same Scale', 	$
	uvalue = 'SPSS')
  t1 = widget_button(base, value = 'Default Range and Scale', 	$
	uvalue = 'NONE')

  widget_control, base, /realize
  xmanager, 'ranges', base, event_handler='range_event', /modal

  maxminque = mmq
  change.range = cr

  return

end

