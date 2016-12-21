pro new_symbol, screen = screen

  common extras, unconpts, columns, symtype, maxminque, rangepl, numofrang,   $
	         titleofplot, plotwhere, psfile, placement, numofgra,	      $
		 psfont, offdiv, offmul, maxran, checkinvar, filename, nb,    $
		 publish, checkstat, statfile, xyzpos, stat, lat, lon, alpha, $
		 mln, remmean, velocity, beam, step, crb, clr, draw,	      $
		 final_x, final_y, arrow_x, arrow_y, arrow, str_x, str_y,     $
		 str_s, str_o,str_l, str, checkselect, percent, inarow, ext
  common subsym, symbol_type

  temp = setup_symbol(plottype, screen)
  symbolbase = temp.base
                 
  if n_elements(screen) eq 0 then begin
    widget_control, symbolbase, /realize
    xmanager, 'symbols', symbolbase, event_handler='symbol_event', /modal
  endif else begin
    print, 'Enter Number of symbol you would like : '
    num = 0
    read, num
    if (num ge 1) and (num le n_elements(temp.types)) then		$
      symbol_type = temp.types(num-1)
  endelse 
    
  symtype = symbol_type

  return

end

