
;-----------------------------------------------------------------------------
;
;  procedure set_up_color
;
;	creates a 200 point color table if color is selected
;

PRO set_up_color

  red = fltarr(201)
  green = fltarr(201)
  blue = fltarr(201)

  red(1) = 255.0
  green(1) = 255.0
  blue(1) = 255.0

  for i=2,39 do begin

    red(i) = 200.0-5.0*float(i)
    green(i) = 200.0 - 5.0*float(i)
    blue(i) = 255.0 - float(i-1)*55.0/38.0

  endfor

  for i=40,79 do begin

    red(i) = 0.0
    green(i) = float(i-40)*200.0/39.0
    blue(i) = 200.0

  endfor

  for i=80,105 do begin

    red(i) = 0.0
    green(i) = 200.0
    blue(i) = 200.0 - float(i-79)*200.0/26.0

  endfor

  for i=106,132 do begin

    red(i) = float(i-105)*255.0/27.0
    green(i) = 200.0 + float(i-105)*55.0/27.0
    blue(i) = 0.0

  endfor

  for i=133,159 do begin

    red(i) = 255.0
    green(i) = 255.0 - float(i-133)*255.0/26.0
    blue(i) = 0.0

  endfor

  for i=160,200 do begin

    red(i) = 255
    green(i) = float(i-160)*5.0
    blue(i) = float(i-160)*5.0

  endfor

  tvlct, red,green,blue

  RETURN

END

