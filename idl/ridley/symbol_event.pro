pro symbol_event, ev

  common subsym, symbol_type

  widget_control, ev.id, get_uvalue=value
  if (n_elements(value) eq 0) then value = ''
  name = strmid(tag_names(ev, /structure_name), 7, 1000)

  case (name) of

    "BUTTON": begin

      case (value) of

	'Solid Lines' : symbol_type = 0
	'Line With Crosses' : symbol_type = -1
	'Line With Stars' : symbol_type = -2
	'Line With Diamonds' : symbol_type = -4
	'Line With Triangles' : symbol_type = -5
	'Line With Squares' : symbol_type = -6
	'Line With X' : symbol_type = -7
	'Crosses' : symbol_type = 1
	'Stars' : symbol_type = 2
	'Diamonds' : symbol_type = 4
	'Triangles' : symbol_type = 5
	'Squares' : symbol_type = 6
	'Dots' : symbol_type = 3
	'Xs' : symbol_type = 7
	'Circles' : symbol_type = 8

    endcase

    end

  else:
  endcase

  widget_control, ev.top, /destroy
  return

end

