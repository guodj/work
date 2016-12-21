function get_integer, range, percent, title=title

  common pblock, per

  if n_elements(percent) eq 0 then percent = 1.0

  if n_elements(title) eq 0 then title = 'Choose Percentage'

  p_base = widget_base(/column)
  l = widget_label(p_base, value = title)
  s = widget_slider(p_base, min = range(0), max = range(1), 		$
	value = percent)
  d = widget_button(p_base, value = 'Done', uvalue = 'DONE')

  widget_control, p_base, /realize
  xmanager, 'percent',p_base, event_handler='percent_event', /modal

  return, per

end
