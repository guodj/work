
pro figure_out_positions, n, xxyy, pos, ppp = ppp

   if !p.charsize eq 0.0 then !p.charsize=1.0

   if (n_elements(ppp) eq 0) then begin
     ppp   = 6
     ny = 3
   endif else begin
     ny = ppp
   endelse

   space = max([float(!d.y_ch_size)/float(!d.y_size),$
                float(!d.x_ch_size)/float(!d.x_size)])*3.0*!p.charsize

   pos_space, ppp, space, sizes, ny = 3

   if (n eq 1 or n eq 2) then begin

     if (ppp eq 6) then begin
       get_position, ppp, space, sizes, (n-1)*2, pos1, /rect
       get_position, ppp, space, sizes, (n-1)*2+1, pos2, /rect

       pos2(2) = pos2(2) - (pos2(2) - pos2(0))*0.1

       pos = [pos1(0), pos1(1), pos2(2), pos1(3)]
     endif else begin
       get_position, ppp, space, sizes, n-1, pos, /rect
       pos(2) = pos(2) - (pos(2) - pos(0))*0.1
     endelse

     width  =  max([xxyy(0),xxyy(2)]) - min([xxyy(0),xxyy(2)])
     height =  max([xxyy(1),xxyy(3)]) - min([xxyy(1),xxyy(3)])
     aspectx = width/height

  endif else begin

     get_position, ppp, space, sizes, n + 1, pos
     aspectx = 1.0

  endelse

  aspectpos = (pos(2)-pos(0))/(pos(3)-pos(1)) $
               *float(!d.x_size)/float(!d.y_size)

  aspectratio = aspectpos/aspectx

  if aspectratio gt 1 then begin
    posmid=(pos(2)+pos(0))/2.
    posdif=(pos(2)-pos(0))/2.
    pos(0)=posmid - posdif/aspectratio
    pos(2)=posmid + posdif/aspectratio
  endif else begin
    posmid=(pos(3)+pos(1))/2.
    posdif=(pos(3)-pos(1))/2.
    pos(1)=posmid - posdif*aspectratio
    pos(3)=posmid + posdif*aspectratio
  endelse

  if (n eq 2) then pos([1,3]) = pos([1,3]) ;+ space/2.0

  if (n eq 3) then pos([0,2]) = pos([0,2]) + space/3.0
  if (n eq 4) then pos([0,2]) = pos([0,2]) - space/3.0

  return

end

