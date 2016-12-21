function sugetstr, str_to_get

  base = widget_base(title='Get'+str_to_get, /column) 

  dum = widget_label(base, value = 'Enter '+str_to_get+' : ')
  instring = widget_text(base, /editable, xsize=30, ysize=1, /frame)

  return, [base, instring]

end

