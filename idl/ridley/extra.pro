pro extra

  common type, plottype, types
  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common comdis, displays
  common text_id, intext, inps, instring, inmaxx, inmaxy, inminx, inminy

  temp = setup_extra(plottype)
  extrabase = temp(0)
  instring = temp(1)
                 
  widget_control, extrabase, /realize
  xmanager, 'extras', extrabase, event_handler='extra_event', /modal

  return

end

