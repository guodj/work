pro percent_event, ev

  common pblock, per

  widget_control, ev.id, get_uvalue=value
  if (n_elements(value) eq 0) then value = ''
  name = strmid(tag_names(ev, /structure_name), 7, 1000)

  case (name) of

    "BUTTON": begin

      widget_control, /destroy, ev.top

      return

    end

    "SLIDER": begin

      value = ev.value

      per = value
	      
    end

  else:
  endcase

end

