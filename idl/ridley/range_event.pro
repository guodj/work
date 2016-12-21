pro range_event, ev

  common range_type, mmq, cr

  widget_control, ev.id, get_uvalue=value
  if (n_elements(value) eq 0) then value = ''
  name = strmid(tag_names(ev, /structure_name), 7, 1000)

  case (name) of

    "BUTTON": begin

      case (value) of

	'APSR' : begin
	  mmq = 1
	  cr = 1
	end

	'SPSR' : begin
	  mmq = 2
	  cr = 1
	end

	'APSS' : begin
	  mmq = 3
	  cr = 1
	end

	'SPSS' : begin
	  mmq = 4
	  cr = 1
	end

	'NONE' : begin
	  mmq = 5
	  cr = 1
	end

    endcase

    end

  else:
  endcase

  widget_control, ev.top, /destroy
  return

end

