PRO prompt_for_next

  if !d.name ne 'PS' then					$
    dum = wmenu(['There is another time interval...',		$
                 '     Press the mouse botton      '],		$
	        init = 1)

  RETURN

END


