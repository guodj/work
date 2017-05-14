;
; set_position
;
; used in conjunction with set_space. Determines the position of the current
; plotting region, given the output parameters from set_space.
;
; Input parameters:
; nb, space, bs, nbx, nby, xoff, yoff, xf, yf - Outputs from set_space
; pos_num - the number of the plot, ranges from 0 : bs-1
;
; Output parameters:
;
; pos - the position of the plot, used in the plot command
;
; modified to make rectangles on Jan 2, 1998

pro set_position, sizes, xipos, yipos, pos, rect = rect,		$
		  xmargin = xmargin, ymargin = ymargin

  nb = sizes.ppp
  space = sizes.space

  yf2 = sizes.yf
  yf = sizes.yf*(1.0-space)
  xf2 = sizes.xf
  xf = sizes.xf*(1.0-space)

  if keyword_set(rect) then begin

    if keyword_set(xmargin) then xmar = xmargin(0) 			$
    else xmar = space/2.0

    if keyword_set(ymargin) then ymar = ymargin(0) 			$
    else ymar = space/2.0

    xbuffer = 3.0*float(!d.x_ch_size)/float(!d.x_size) * !p.charsize +space/4.0
    xtotal = 1.0 - (space*float(sizes.nbx-1) + xmar + xf2*space/2.0) - xbuffer
    xbs = xtotal/(float(sizes.nbx)*xf)

    xoff = xmar - xf2*space/2.0 + xbuffer - space/4.0

    ybuffer = 3.0*float(!d.y_ch_size)/float(!d.y_size) * !p.charsize
    ytotal = 1.0 - (space*float(sizes.nby-1) + ymar + yf2*space/2.0) - ybuffer
    ybs = ytotal/(float(sizes.nby)*yf)

    yoff = space/4.0

  endif else begin

    xbs  = sizes.bs
    xoff = sizes.xoff
    ybs  = sizes.bs
    yoff = sizes.yoff

  endelse

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0 + xbs*xf

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space + xbs*xf

  ypos0 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf + ybs*yf) - yoff
  ypos1 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf) - yoff

  pos= [xpos0,ypos0,xpos1,ypos1]

  RETURN

END

