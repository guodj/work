function getstr, str_to_get

  common gs_common, instring, r_string

  temp = sugetstr(str_to_get)
  strbase = temp(0)
  instring = temp(1)
                 
  widget_control, strbase, /realize
  xmanager, 'getstr', strbase, event_handler='getstre', /modal

  return, r_string

end

