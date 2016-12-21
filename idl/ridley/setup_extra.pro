function setup_extra, plottype

  extrabase = widget_base(title='Extra! Extra!', /column) 

  dum = widget_button(extrabase, value="Done", uvalue = "DONE")
  dum = widget_label(extrabase, value = 'Label for Plot : ')
  instring = widget_text(extrabase, /editable, xsize=30, ysize=1, /frame)
  dum = widget_button(extrabase, value="Remove Last Label",		$
	uvalue="LA_CLEAR_1")
  dum = widget_button(extrabase, value="Remove All Labels",		$
	uvalue="LA_CLEAR_ALL")

  dum = widget_label(extrabase, value = '-')

  dum = widget_button(extrabase, value="Draw Curve on Plot",		$
	uvalue="DRAW")
  dum = widget_button(extrabase, value="Remove Last Curve",		$
	uvalue="CLEAR_1")
  dum = widget_button(extrabase, value="Remove All Curves",		$
	uvalue="CLEAR_ALL")
  dum = widget_button(extrabase, value="Read Curves from File",		$
	uvalue="IN_CURVE")
  dum = widget_button(extrabase, value="Save Curves to File",		$
	uvalue="OUT_CURVE")

  dum = widget_label(extrabase, value = '-')

  dum = widget_button(extrabase, value="Draw Line on Plot",	 	$
	uvalue="LINE")
  dum = widget_button(extrabase, value="Remove Last Line",		$
	uvalue="LI_CLEAR_1")
  dum = widget_button(extrabase, value="Remove All Lines",		$
	uvalue="LI_CLEAR_ALL")
  dum = widget_button(extrabase, value="Read Lines from File",		$
	uvalue="IN_LINE")
  dum = widget_button(extrabase, value="Save Lines to File",		$
	uvalue="OUT_LINE")

  dum = widget_label(extrabase, value = '-')

  dum = widget_button(extrabase, value="Draw Arrows on Plot",	 	$
	uvalue="ARROW")
  dum = widget_button(extrabase, value="Remove Last Arrow",		$
	uvalue="AR_CLEAR_1")
  dum = widget_button(extrabase, value="Remove All Arrows",		$
	uvalue="AR_CLEAR_ALL")

  dum = widget_label(extrabase, value = '-')

;  type_8   Azimuth Scans
;  type_9
;  type_13

  if (plottype eq 8) or							$
     (plottype eq 9) or							$
     (plottype eq 13) then begin

    dum = widget_button(extrabase, value="Add Reversal Line On Plot",	$
	  uvalue="CRB")

  endif

;  type_8   Azimuth Scans
;  type_13 
;  type_21  Magnetometer Clock Dial

  if plottype eq 8 or							$
     plottype eq 13 or							$
     plottype eq 21 then begin

    dum = widget_button(extrabase, value='Change Plotting Display',	$
	  uvalue = 'DCHANGE')

  endif

;  xpdmenu, ['/PS Font/ {', 						$
;	    '/Courier/', 						$
;	    '/Courier, Bold/', 						$
;	    '/Courier, Oblique/', 					$
;	    '/Courier, Bold, Oblique/',					$
;	    '/Helvetica/',						$
;	    '/Helvetica, Bold/',					$
;	    '/Helvetica, Oblique/',					$
;	    '/Helvetica, Bold, Oblique/',				$
;	    '/Avantgarde, Book/',					$
;	    '/Avantgarde, Book, Oblique/',				$
;	    '/Avantgarde, Demi/',					$
;	    '/Avantgarde, Demi, Oblique/',				$
;	    '/Schoolbook/',						$
;	    '/Schoolbook, Bold/',					$
;	    '/Schoolbook, Italic/',					$
;            '/Schoolbook, Bold, Italic/',				$
;	    '/Times/',							$
;	    '/Times, Bold/',						$
;	    '/Times, Italic/',						$
;	    '/Times, Bold, Italic/',					$
;	    '}'],extrabase

  return, [extrabase, instring]

end

