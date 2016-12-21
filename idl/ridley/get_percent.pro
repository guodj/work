function get_percent, range, percent

  common pblock, per

  if n_elements(percent) eq 0 then percent = 1.0

  p_base = widget_base(/column)
  l = widget_label(p_base, value = 'Choose Percentage')
  s = widget_slider(p_base, min = range(0), max = range(1), 		$
	value = fix(percent*100.0))
  d = widget_button(p_base, value = 'Done', uvalue = 'DONE')

  widget_control, p_base, /realize
  xmanager, 'percent',p_base, event_handler='percent_event', /modal

  return, float(per)/100.0

end
