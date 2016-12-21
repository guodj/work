pro break_down, array

  npts = n_elements(array)
  maxi = max(array)
  mini = min(array)

  print, 'There are ',npts,' points in this array.'
  print, 'The maximum is : ',maxi
  print, 'The minimum is : ',mini

  nmini = mini
  nmaxi = maxi

  totalper = 100.0

  while totalper gt 15.0 do begin

    per = 100.0

    fac = 1.0

    while (per lt 5.0) or (per gt 15.0) do begin

      nmaxi = nmini + fac*(nmaxi-nmini)

      loc = where(array ge nmini and array le nmaxi, count)

      per = 100.0*float(count)/float(npts)

      if per gt 15.0 then fac = 0.5
      if per lt 5.0 then fac = 1.5

      print, nmini,nmaxi,per, fac, count

    endwhile

    print, per,' % of the array is between',nmini,' and',nmaxi 

    nmini = nmaxi
    nmaxi = maxi
    totalper = totalper - per

  endwhile

  print, totalper,' % of the array is between',nmini,' and',maxi 

  return

end