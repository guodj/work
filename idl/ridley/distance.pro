;
; function distance
;
;  computes the distance (in Km) between any two points on a sphere, given
;  the latitude and longitude of the points
;

function distance, lat1, long1, lat2, long2

  b = (90.0 - lat1)*!dtor
  a = (90.0 - lat2)*!dtor
  cc = abs(long1-long2)*!dtor
  
;  if (a eq 0) or (b eq 0) or (cc eq 0) then begin
;    dis = abs(b-a)*110.0*180.0/!pi
;  endif else begin
    cotb = 1.0/tan(b)
    cotcb = (cotb*sin(a) - cos(a)*cos(cc))/sin(cc)
    cb = atan(1.0/cotcb)
    c = asin(sin(b)*sin(cc)/sin(cb))
    dis = abs(c)*111.0/!dtor
;  endelse

  return, dis

end
