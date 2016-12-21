pro move_bullets, bullet, mouse, maxbul, b_speed, cent, rand

	common positions, grid
	common objects, ob

	if (mouse eq 1) and (n_elements(bullet) lt maxbul*2) then begin

	  if bullet(0) eq -1 then begin
	    bullet = fltarr(1,2)
	    bullet(0,*) = [ob.oldman.x,ob.oldman.y+15]
	  endif else begin
	    bullet = [bullet, fltarr(1,2)]
	    bullet(n_elements(bullet(*,0))-1,*) = [ob.oldman.x,ob.oldman.y+15]
	  endelse

	  mouse = 0

	endif

	nb = n_elements(bullet)
	if nb gt 1 then begin

	  for i=0,nb/2-1 do begin
	    tv, ob.dark_bullet, bullet(i,0)-2,bullet(i,1)-4
	    bullet(i,1) = bullet(i,1) + b_speed
	    xpos = bullet(i,0)/20
	    ypos = bullet(i,1)/20
	    if (ypos lt n_elements(grid(0,*))) and (ypos gt 0) then begin
	      if grid(xpos,ypos) eq 0 then				$
	        tv, ob.bullet, bullet(i,0)-2,bullet(i,1)-4		$
	        else begin
	          case(grid(xpos,ypos))of
	            1:begin
	              tv, ob.half_shroom1, 20*xpos, 20*ypos
	              grid(xpos,ypos) = 3
	            endcase
	            2:begin
	              tv, ob.half_shroom2, 20*xpos, 20*ypos
	              grid(xpos,ypos) = 4
	            endcase
	            3:begin
	              tv, ob.dark_shroom, 20*xpos, 20*ypos
	              grid(xpos,ypos) = 0
	            endcase
	            4:begin
	              tv, ob.dark_shroom, 20*xpos, 20*ypos
	              grid(xpos,ypos) = 0
	            endcase
	            5:begin
	              if rand(0) lt 0.5 then begin
	                tv, ob.mushroom1, 20*xpos, 20*ypos
	                grid(xpos,ypos) = 1
	              endif else begin
	                tv, ob.mushroom2, 20*xpos, 20*ypos
	                grid(xpos,ypos) = 2
	              endelse
; break apart the centipede
		      loc = where(cent(0,*,0) ne 0,nc)
		      done1 = 0
		      ii = 0
		      while done1 eq 0 do begin
			li = loc(ii)
			lx = where(cent(*,li,0) ne 0,len)
			jj = 0
			done2 = 0
			while done2 eq 0 do begin
			  if (cent(jj,li,3) eq xpos) and 	$
			     (cent(jj,li,4) eq ypos) then begin
			    done2 = 1
			    done1 = 1
			    if len-1 gt jj then begin
			      newc = where(cent(0,*,0) eq 0)
			      nc = newc(0)
			      for k=0,len-jj-2 do 		$
			        cent(k,nc,*) = cent(jj+1+k,li,*)
			      if cent(0,nc,1) eq -1 then	$ 
			        cent(0,nc,0) = 1		$
			      else				$
			        cent(0,nc,0) = 3
			      cent(jj:len-1,li,*) = 0
			    endif else cent(jj,li,*) = 0
			  endif
			  jj = jj+1
			  if jj ge len then done2 = 1
			endwhile
			ii = ii +1
			if ii ge nc then done1 = 1
		      endwhile
	            endcase
	          endcase
	        bullet(i,0) = -1
	        bullet(i,1) = -1
	      endelse
	    endif else begin
	      bullet(i,0) = -1
	      bullet(i,1) = -1
	    endelse
	  endfor

	  loc = where(bullet gt 0, count)
	  if count gt 0 then begin
	    loc = where(bullet(*,0) gt 0, count)
	    x = intarr(count)
	    x = bullet(loc,0)
	    y = intarr(count)
	    y = bullet(loc,1)
	    bullet = intarr(count,2)
	    bullet(0:count-1,0) = x
	    bullet(0:count-1,1) = y
	  endif else bullet = -1

	endif

	return

end
