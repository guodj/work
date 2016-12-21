function setup_ranges, rs, numofrang

  font = '-adobe-helvetica-medium-r-*-*-14-140*'

  b = widget_base(title='Determine '+rs, /column)

  out = {base : b, point : lonarr(numofrang,2)}

  dum = widget_button(b, value = 'Done', uvalue = 'DONE', font=font)
  
  if rs eq 'Scale' then begin

    dum = widget_label(b, value = 					$
	'Since you have selected scale, the maximum and minimum',	$
	font=font)
    dum = widget_label(b, value = 					$
	'values can float, but you can select the range which is',	$
	font=font)
    dum = widget_label(b, value = 					$
	'maximum - minimum', font=font)

    for i=0,numofrang-1 do begin

      subb = widget_base(b, /row)

      dum = widget_label(subb, value = 'Range # '+tostr(i+1)+' : ',	$
	font=font)
      out.point(i,0) = widget_text(subb, /editable, xsize = 7, 		$
	ysize = 1, /frame, font=font)

    endfor

  endif else begin

    for i=0,numofrang-1 do begin

      subb = widget_base(b, /row)

      dum = widget_label(subb, value = 'Min Y # '+tostr(i+1)+' : ',	$
	font=font)
      out.point(i,0) = widget_text(subb, /editable, xsize = 7, 		$
	ysize = 1, /frame, font=font)

      dum = widget_label(subb, value = 'Max Y # '+tostr(i+1)+' : ',	$
	font=font)
      out.point(i,1) = widget_text(subb, /editable, xsize = 7, 		$
	ysize = 1, /frame, font=font)

    endfor

  endelse

  return, out

end


