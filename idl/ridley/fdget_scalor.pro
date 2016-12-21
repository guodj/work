pro fdget_scalor, id, mrow, col, uttime, data

  common ffinfo, header
 
; position the file pointer at the right row:

  point_lun, header.unit(id), long(mrow)*header.rl(id)

; read in the time

  uttime = double(0.0)
  readu, header.unit(id), uttime

; position the file pointer at the right column:

  len = header.loc(id,col) - 8
  if len gt 0 then 							$
    point_lun, header.unit(id), long(mrow)*header.rl(id) + len + 8

; create scalor:

  if (header.type(id,col) eq 't') then indata = double(0.0)		$
  else if (header.type(id,col) eq 'r') then indata = 0.0		$
  else if (header.type(id,col) eq 'i') then indata = long(0)		$
  else if (header.type(id,col) eq 'l') then indata = long(0)		$
  else if (header.type(id,col) eq 'c') then indata = ' '

  if n_elements(indata) eq 0 then begin

    print, 'Data type not defined. Can not read file!'

    indata = -1

  endif else begin

    readu, header.unit(id), indata

  endelse

  data = indata

  return

end
