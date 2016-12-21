function chopn,value,n

  if value lt 0 then begin

    strng = '0'+strcompress(string(fix(-1.0*value)), /remove_all)
    strng = '-'+chopr(strng,n)

  endif else 					$
    strng = chopr('0'+strcompress(string(fix(value)), /remove_all),n)

  return, strng

end

