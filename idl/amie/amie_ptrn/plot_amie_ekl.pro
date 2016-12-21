
files = ask('files to plot','*_sum')

psfile = ask('ps file name','ratio.ps')

filelist = findfile(files)
nfiles = n_elements(filelist)
if (nfiles gt 0 and strlen(filelist(0)) gt 0) then begin

  itime = intarr(6)
  tmp   = fltarr(11)
  data  = fltarr(1440*nfiles,11)
  time  = dblarr(1440*nfiles)
  line  = ''

  for ifile=0L,nfiles-1 do begin

    print, "Reading file ",filelist(ifile)
    openr,1,filelist(ifile)
    readf,1,line
    for i=0L,1440-1 do begin
      readf,1,itime,tmp
      c_a_to_r, itime, rtime
      data(ifile*1440L+i,*) = tmp
      time(ifile*1440L+i)   = rtime
    endfor
    close,1

  endfor

  stime = min(time)
  etime = max(time)

  get_imf, stime, etime, imfdata, dt = 60.0

  ekl = imfdata.ekl

  time_axis, stime, etime, btr, etr, $
           xtickname, xtitle, xtickv, xminor, xtickn

  title = strmid(xtitle,0,strpos(xtitle,"Univ")-1)

  time = time - stime

  cpcp = reform(data(*,0))

  eklmax = 30
  npts_min = 5

  dekl = 0.25
  nekl = eklmax/dekl + 1

  eklbin  = findgen(nekl)*dekl
  cpcpbin = fltarr(nekl)
  cpcpstd = fltarr(nekl)

  swnbin = fltarr(nekl)
  swnstd = fltarr(nekl)

  npts    = fltarr(nekl)

  maxekl = 0.0

  for i=0,nekl-1 do begin
      l = where(ekl ge eklbin(i)-dekl/2 and ekl lt eklbin(i)+dekl/2 and $
                abs(imfdata.bz) gt 0.01 and abs(imfdata.by) gt 0.01, count)

      if (count ge npts_min) then begin
          npts(i) = count

          cpcpbin(i) = mean(cpcp(l))
          cpcpstd(i) = stddev(cpcp(l))

          swnbin(i) = mean(imfdata.dynp(l)/imfdata.magp(l))
          swnstd(i) = stddev(imfdata.dynp(l)/imfdata.magp(l))

          if (abs(eklbin(i)) gt maxekl) then maxekl=abs(eklbin(i))
      endif
  endfor

  maxekl = max([10.0,maxekl])

  setdevice, psfile,'p',4

  plotdumb

  ppp = 4
  space = 0.01
  pos_space, ppp, space, sizes, ny = ppp

  get_position, ppp, space, sizes, 0, pos, /rect
  pos(0) = pos(0) + 0.075

  loc = where(npts ge npts_min, c)

  plot_io, eklbin(loc), npts(loc), xrange = [0.0,maxekl], psym = 2, $
    pos = pos, /noerase, xtickname = strarr(10)+' ', ytitle = 'Points', $
    charsize = 1.2, thick = 2.0, xstyle = 1, title = title

  get_position, ppp, space, sizes, 1, pos, /rect
  pos(0) = pos(0) + 0.075

  ymax = max([30.0,max(swnbin+swnstd)])
;  ymax = max(swnbin+swnstd)

  plot_io, eklbin, swnbin, min_val=0.01, xrange = [0.0,maxekl], psym = 2, $
    pos = pos, xtickname = strarr(10)+' ', ytitle = 'Pdyn (P)', $
    /noerase, charsize = 1.2, thick = 2.0, xstyle = 1, yrange = [0,ymax]

  for l=0,c-1 do begin
      i = loc(l)
      oplot, [eklbin(i), eklbin(i)], $
        [swnbin(i)-swnstd(i),swnbin(i)+swnstd(i)], thick = 2.0
      oplot, eklbin(i)+[-dekl,dekl]/2, swnbin(i)+[ swnstd(i), swnstd(i)], $
        thick = 2.0
      oplot, eklbin(i)+[-dekl,dekl]/2, swnbin(i)+[-swnstd(i),-swnstd(i)], $
        thick = 2.0
  endfor

  get_position, ppp, space, sizes, 2, post, /rect
  get_position, ppp, space, sizes, 3, posb, /rect
  pos = posb
  pos(3) = post(3)
  pos(0) = pos(0) + 0.075
  pos(1) = pos(1) + 0.06

  ymax = max([200.0,max(cpcpbin+cpcpstd)])

  plot, eklbin, cpcpbin, min_val=0.01, xrange = [0.0,maxekl], psym = 2, $
    pos = pos, xtitle = 'Ekl (mV/m)', ytitle = 'Potential (kV)', $
    /noerase, charsize = 1.2, thick = 2.0, xstyle = 1, yrange = [0,ymax]

  for l=0,c-1 do begin
      i = loc(l)
      oplot, [eklbin(i), eklbin(i)], $
        [cpcpbin(i)-cpcpstd(i),cpcpbin(i)+cpcpstd(i)], thick = 2.0
      oplot, eklbin(i)+[-dekl,dekl]/2, cpcpbin(i)+[ cpcpstd(i), cpcpstd(i)], $
        thick = 2.0
      oplot, eklbin(i)+[-dekl,dekl]/2, cpcpbin(i)+[-cpcpstd(i),-cpcpstd(i)], $
        thick = 2.0
  endfor

  loc = where(npts ge npts_min and eklbin ge 0.0,c)

  slope_s = 15.0
  yint_s  = 40.0

  test = eklbin(loc)*slope_s + yint_s
  error_s = 0.0
  for l=0,c-1 do begin
      i = loc(l)
      w = alog10(npts(i))^3.0/cpcpstd(i)
      error_s = error_s + abs(cpcpbin(i) - test(l))*w
  endfor

  slope = slope_s
  yint  = yint_s

  for n=1,1000 do begin

    slope = slope_s + (randomu(seed)-0.5)*2.0 * slope_s/5.0
    yint  = yint_s  + (randomu(seed)-0.5)*2.0 *  yint_s/5.0

    test = eklbin(loc)*slope + yint
    error = 0.0
    for l=0,c-1 do begin
        i = loc(l)
        w = alog10(npts(i))^3.0/cpcpstd(i)
        error = error + abs(cpcpbin(i) - test(l))*w
    endfor

    if (error lt error_s) then begin
      error_s = error
      slope_s = slope
      yint_s  = yint
      print, n, slope_s, yint_s, error_s
    endif

  endfor

  oplot, eklbin(0:nekl-1), eklbin(0:nekl-1)*slope_s+yint_s, $
    thick = 2.0  

  xyouts, charsize = 1.2, 3.0*dekl, ymax-ymax*0.05, $
    'Slope : '+string(slope_s,format='(f5.1)')+' (kV/mV/m)'
  xyouts, charsize = 1.2, 3.0*dekl, ymax-ymax*0.1, $
    'Y-Int. : '+string(yint_s,format='(f5.1)')+' (kV)'
  xyouts, charsize = 1.0, 3.0*dekl, ymax-ymax*0.15, $
    'Error : '+string(error_s,format='(f7.1)')

  closedevice

endif

end
