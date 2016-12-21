pro get_var, fbase, vars, n, numofvar

  common titles, var, nvar

  nvar = numofvar

  var = vars(n,*)
  widget_control, fbase, /realize
  xmanager, 'getvar', fbase, event_handler='get_var_event', /modal
  numofvar = nvar
  vars(n,*) = var

  return

end
