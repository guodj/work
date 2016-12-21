pro d_colortable, pos, mm

  plot, [0,1],[mm(0),mm(1)], /nodata, 			$
	xstyle=5, ystyle=1, pos = pos,			$
	/noerase, yticklen=-0.1

  oplot, [0,1],[mm(0),mm(0)]
  oplot, [0,1],[mm(1),mm(1)]

  pos2 = pos
  pos2(0) = pos2(0) + 0.001
  pos2(1) = pos2(1) + 0.001
  pos2(2) = pos2(2) - 0.001
  pos2(3) = pos2(3) - 0.001

  plot, [0,1],[0,255], pos=pos2, /nodata, /noerase, 	$
	xstyle=5,ystyle=5

  dy = 255.0/256.0

  x = [0.0,0.0,1.0,1.0,0.0]
  y = [0.0,1.0,1.0,0.0,0.0]*dy

  for i=0,255 do					$
    polyfill, x, y+dy*i, color = i
    
  return

end
