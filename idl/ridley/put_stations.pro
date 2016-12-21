PRO Put_Stations, disfac, elev, change

  if n_elements(change) eq 0 then change = 1.0

  x = findgen(18)*!pi/9.

  usersym, 0.75*cos(x), .75*sin(x), /fill

  oplot, [0.0], [0.0], psym = 8

  xyouts, 2.0, 2.0, 'STF',alignment=0.0, charsize = change 

  alpha = - 1 * 24.0 * !pi / 180.0

  elev = elev * !pi / 180.0

  sinalpha = sin (alpha)
  cosalpha = cos (alpha)

; UMQ

  oldx = 140.0
  oldy = 388.0

  xyang = atan(oldx/oldy)
  if oldx lt 0.0 then xyang = xyang + !pi
  oldr = ((oldx^2 + oldy^2)^0.5) * cos(elev)

  oldx = oldr * sin(xyang)
  oldy = oldr * cos(xyang)

  newx = (oldx * cosalpha + oldy * sinalpha) * disfac
  newy = (oldy * cosalpha - oldx * sinalpha) * disfac

  oplot, [newx], [newy], psym = 8

  xyouts, newx-2.0, newy+2.0, 'UMQ',alignment=1.0, charsize = change 

; GDH

  oldx = 15.0
  oldy = 273.0

  xyang = atan(oldx/oldy)
  if oldx lt 0.0 then xyang = xyang + !pi
  oldr = ((oldx^2 + oldy^2)^0.5) * cos(elev)

  oldx = oldr * sin(xyang)
  oldy = oldr * cos(xyang)

  newx = (oldx * cosalpha + oldy * sinalpha) * disfac
  newy = (oldy * cosalpha - oldx * sinalpha) * disfac

  oplot, [newx], [newy], psym = 8

  xyouts, newx, newy+3.0, 'GDH',alignment=0.0, charsize = change 

; ITI

  oldx = -118.0
  oldy = 0.0

  oldx = oldx * cos(elev)

  newx = (oldx * cosalpha + oldy * sinalpha) * disfac
  newy = (oldy * cosalpha - oldx * sinalpha) * disfac

;  oplot, [newx], [newy], psym = 8

;  xyouts, newx-2.0, newy+2.0, 'ITI',alignment=1.0, charsize = change 

; ATU

  oldx = -58.0
  oldy = 149.0

  xyang = atan(oldx/oldy)
  oldr = ((oldx^2 + oldy^2)^0.5) * cos(elev)

  oldx = oldr * sin(xyang)
  oldy = oldr * cos(xyang)

  newx = (oldx * cosalpha + oldy * sinalpha) * disfac
  newy = (oldy * cosalpha - oldx * sinalpha) * disfac

  oplot, [newx], [newy], psym = 8

  xyouts, newx-2.0, newy+2.0, 'ATU',alignment=1.0, charsize = change 

; SKT

  oldx = -110.0
  oldy = -30.0

  xyang = atan(oldx/oldy)
  if oldx lt 0.0 then xyang = xyang + !pi
  oldr = ((oldx^2 + oldy^2)^0.5) * cos(elev)

  oldx = oldr * sin(xyang)
  oldy = oldr * cos(xyang)

  newx = (oldx * cosalpha + oldy * sinalpha) * disfac
  newy = (oldy * cosalpha - oldx * sinalpha) * disfac

;  oplot, [newx], [newy], psym = 8

;  xyouts, newx-2.0, newy+2.0, 'SKT',alignment=1.0, charsize = change 

; DYE

  oldx = 144.0
  oldy = -129.0

  xyang = atan(oldx/oldy)
  if oldy lt 0.0 then xyang = xyang + !pi
  oldr = ((oldx^2 + oldy^2)^0.5) * cos(elev)

  oldx = oldr * sin(xyang)
  oldy = oldr * cos(xyang)

  newx = (oldx * cosalpha + oldy * sinalpha) * disfac
  newy = (oldy * cosalpha - oldx * sinalpha) * disfac

;  oplot, [newx], [newy], psym = 8

;  xyouts, newx-2.0, newy+2.0, 'DYE',alignment=1.0, charsize = change 

; GHB

  oldx = -186.0
  oldy = -258.0

  xyang = atan(oldx/oldy)
  if oldx lt 0.0 then xyang = xyang + !pi
  oldr = ((oldx^2 + oldy^2)^0.5) * cos(elev)

  oldx = oldr * sin(xyang)
  oldy = oldr * cos(xyang)

  newx = (oldx * cosalpha + oldy * sinalpha) * disfac
  newy = (oldy * cosalpha - oldx * sinalpha) * disfac

  oplot, [newx], [newy], psym = 8

  xyouts, newx-2.0, newy+2.0, 'GHB',alignment=1.0, charsize = change 

  RETURN

END


