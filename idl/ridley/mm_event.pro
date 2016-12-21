pro mm_event, ev

  common crange, ranges, baseid

  name = strmid(tag_names(ev,/structure_name),7,1000)

  if name eq 'BUTTON' then begin

    widget_control, /destroy, ev.top
    return

  endif else begin

    widget_control, ev.id, get_value = value

    loc = where(baseid.point eq ev.id, count)

    if count gt 0 then ranges(loc) = value 			$
    else print, 'Error in program, count = 0'

  endelse

end
