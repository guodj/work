;-----------------------------------------------------------------------------
;
;  Select Plot type menu
;
;-----------------------------------------------------------------------------

pro gettypee, ev

  COMMON type, plottype, types

  widget_control, ev.id, get_uvalue=value

  plottype=ev.index + 1

  widget_control, /destroy, ev.top

  return

end

