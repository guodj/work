
pro thermo_compare_sat_plot, lats, times, gitm, data, localtime, $
                             stime, etime, $
                             psfile, title, satellite

  d = [data,gitm]
  maxi = max(abs(d))
  n = n_elements(d)
  l = where(abs(d) gt maxi,c)
  while c lt 0.05*n do begin
     maxi = maxi*0.99
     l = where(abs(d) gt maxi,c)
  endwhile

  maxi = maxi*1.2

  if (min(data) lt 0.0) then mini = -maxi else mini = 0.0

  yrange = [mini,maxi]

  ppp = 3
  space = 0.04
  pos_space, ppp, space, sizes, ny = ppp

  get_position, ppp, space, sizes, 0, pos0, /rect
  pos0(0) = pos0(0) + 0.05
  pos0(2) = pos0(2) - 0.05

  get_position, ppp, space, sizes, 1, pos, /rect
  pos(0) = pos(0) + 0.05
  pos(2) = pos(2) - 0.05

  get_position, ppp, space, sizes, 2, pos2, /rect
  pos2(0) = pos2(0) + 0.05
  pos2(2) = pos2(2) - 0.05

  time_axis, stime, etime,btr,etr, xtickname, xtitle, xtickv, xminor, xtickn

  c_r_to_a, itime, stime
  c_a_to_y_m_d, itime, ymd
  c_r_to_a, itime, etime
  c_a_to_y_m_d, itime, ymd2

  setdevice, psfile, 'p', 5

  makect,'mid'

  levels = 2*findgen(31)*maxi/30.0+mini

  l = where(data gt levels(29),c)
  if (c gt 0) then data(l) = levels(29)
  l = where(gitm gt levels(29),c)
  if (c gt 0) then gitm(l) = levels(29)

  l = where(data lt levels(1),c)
  if (c gt 0) then data(l) = levels(1)
  l = where(gitm lt levels(1),c)
  if (c gt 0) then gitm(l) = levels(1)

  cPos = pos0
  cPos(2) = cPos(2)-0.05
  contour, data, times, lats, pos = cPos, $
           xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
           xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
           ytitle = 'Latitude (deg)', thick = 3, $
           xrange = [btr,etr], levels = levels, /fill

  xyouts, cPos(2), cPos(3)+0.01, 'Local Time : ' + $
          string(localtime,format='(f5.2)'), /norm, align = 1.0

  ctpos = cPos
  ctpos(0) = cPos(2)+0.01
  ctpos(2) = ctpos(0)+0.03
  plotct,254,ctpos,mm(levels),Satellite+title,/right

  cPos = pos
  cPos(2) = cPos(2)-0.05
  contour, gitm, times, lats, pos = cPos, $
           xtickname = xtickname, xtickv = xtickv, ystyle = 1, $
           xminor = xminor, xticks = xtickn, xstyle = 1, charsize = 1.2, $
           ytitle = 'Latitude (deg)', thick = 3, /noerase, $
           xrange = [btr,etr], levels = levels, /fill, $
           xtitle = xtitle

  ctpos = cPos
  ctpos(0) = cPos(2)+0.01
  ctpos(2) = ctpos(0)+0.03

  plotct,254,ctpos,mm(levels),'GITM'+title,/right

  closedevice

end
