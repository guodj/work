
filelist = findfile("-t *.jpg")
file = filelist(0)

file = ask("jpg file to convert to tecplot",file)

isfullmap = fix(ask('whether to map over full sphere (0 or 1)','1'))

read_jpeg, file, image

s = size(image)
if (s(0) eq 3) then color = 1 else color = 0


if (color) then begin
    nLons = n_elements(image(0,*,0))
    nLats = n_elements(image(0,0,*))
endif else begin
    nLons = s(1)
    nLats = s(2)
endelse

if (isfullmap) then begin
    lats = findgen(nLats) / (nLats-1) * !pi - !pi/2
    lons = findgen(nLons) / (nLons-1) * 2*!pi
    r = 1.0
endif else begin
    lats = findgen(nLats) / (nLats-1) * !pi/3 - !pi/6
    lons = findgen(nLons) / (nLons-1) * !pi*3.0/4.0 + .25*!pi
    r = 1.01
endelse

openw,1,file+".dat"

printf,1,"TITLE = ""Earth Surface Map"""
printf,1,"VARIABLES = ""X [R]"", ""Y [R]"", ""Z [R]"",""Earth-Red"", ""Earth-Green"", ""Earth-Blue"""
printf,1,"ZONE T=""Earth"", I="+tostr(nLats)+", J="+tostr(nLons)+", K=1, ZONETYPE=Ordered, DATAPACKING=POINT"

for iLon = 0, nLons-1 do for iLat = 0,nLats-1 do begin

    x = r*cos(lats(iLat))*cos(lons(iLon))
    y = r*cos(lats(iLat))*sin(lons(iLon))
    z = r*sin(lats(iLat))

  if (color) then begin
      printf,1,x,y,z, $
        float(image(0,iLon,iLat))/255.0,$
        float(image(1,iLon,iLat))/255.0,$
        float(image(2,iLon,iLat))/255.0
  endif else begin
      printf,1,x,y,z, $
        1.0-float(image(iLon,iLat))/255.0,$
        1.0-float(image(iLon,iLat))/255.0,$
        1.0-float(image(iLon,iLat))/255.0
  endelse

endfor
close,1

end
