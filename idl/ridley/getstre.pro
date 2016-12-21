pro getstre, ev

  common gs_common, instring, r_string

  widget_control, ev.id, get_value=value

  if ev.id eq instring then r_string = value(0)

  widget_control, ev.top, /destroy

  return

end

