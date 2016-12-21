files = ask('files to plot','*_sum')

psfile = ask('ps file name','pci.ps')

filelist = findfile(files)
nfiles = n_elements(filelist)
if (nfiles gt 0 and strlen(filelist(0)) gt 0) then begin

  itime = intarr(6)
  tmp   = fltarr(12)
  data  = fltarr(1440*nfiles,12)
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

  c_r_to_a, itime, stime
  syear = tostr(itime(0))

  pcifile = getenv('PCI_DATA')+'/p15n'+syear+'.wdc'
  read_pci_wdc, pcifile, pci, pcitime

  loc = where(pcitime ge stime and pcitime le etime,nPci)

  if (nPci gt 0) then begin

    pcitime = pcitime(loc) - stime
    pci = pci(loc)

    pcicpcp = 19.35*pci+8.78
    pciram = -0.12*pci^2 + 2.5*pci + 11.0
    pcirpm = -0.20*pci^2 + 3.0*pci + 12.5

    pciram = pciram * 113.0*1000.0
    pcirpm = pcirpm * 113.0*1000.0

    pciarea = !pi*(0.25*pciram^2 + 0.25*pcirpm^2 + 0.5*(0.5*(pciram+pcirpm))^2)
    pciarea = pciarea/1.0e12

    pcie = (-0.17*pci^2 + 4.5*pci + 9.0) * 1.7

    loc = where(pcicpcp lt 0.0, count)
    if count gt 0 then pcicpcp(loc) = 0.0

    time_axis, stime, etime, btr, etr, $
               xtickname, xtitle, xtickv, xminor, xtickn

    title = strmid(xtitle,0,strpos(xtitle,"Univ")-1)

    time = time - stime

    cpcp = reform(data(*,0))
    area = reform(data(*,10))/1.0e12
    e    = reform(data(*,11)) * 1000.0

    ; let's redistribute the cpcp onto the pci time

    print, "redistributing CPCP"
    redistribute, cpcp, time, pcitime, /average
    print, "redistributing Area"
    redistribute, area, time, pcitime, /average
    print, "redistributing E"
    redistribute, e,    time, pcitime, /average

stop

    setdevice, psfile,'p',4

    plotdumb

    ppp = 6
    space = 0.05
    ysize = (1.0 - space*(ppp-1) )/ppp
    xspace = space

    half = 1.0 - (ysize*3 + space*2)
    space_top = 0.01
    ysize_top = (1.0-half - space_top*2.0)/3

    pos = fltarr(4)
    pos(0) = 0.1
    pos(2) = 0.98

    top = 150.0

    pos(3) = 1.0
    pos(1) = pos(3) - ysize_top

    print, pos(3)-pos(1)

    plot, pcitime, cpcp, $
          xstyle = 1, xticks = xtickn, xtickv = xtickv, $
          /noerase, xtickname = strarr(10)+' ', $
          xminor = xminor, yrange = [0,top], pos = pos, $
          ytitle = 'CPCP (kV)'
    oplot, pcitime, pcicpcp, linestyle = 2

    areatop = 50.0

    pos(3) = pos(1) - space_top
    pos(1) = pos(3) - ysize_top

    print, pos(3)-pos(1)

    plot, pcitime, area, $
          xstyle = 1, xticks = xtickn, xtickv = xtickv, $
          /noerase, xtickname = strarr(10)+' ', $
          xminor = xminor, yrange = [0,areatop], pos = pos, $
          ytitle = 'Area (1e12 m2)'
    oplot, pcitime, pciarea, linestyle = 2

    etop = 100

    pos(3) = pos(1) - space_top
    pos(1) = pos(3) - ysize_top

    print, pos(3)-pos(1), pos(1)

    plot, pcitime, e, $
          xstyle = 1, xticks = xtickn, xtickv = xtickv, $
          /noerase, xtitle = xtitle, xtickname = xtickname, $
          xminor = xminor, yrange = [0,etop], pos = pos, $
          ytitle = 'E-Field (mV/m)'
    oplot, pcitime, pcie, linestyle = 2

    pos(2) = (pos(0)+pos(2))/2.0 - xspace/2.0

    pos(3) = pos(1) - space
    pos(1) = pos(3) - ysize

    xposl = pos

    plot, cpcp, pcicpcp, pos = pos, /noerase, $
          xtitle = 'AMIE CPCP (kV)', ytitle = 'PCI CPCP (kV)', $
          psym = 3, xrange = [0,top], yrange = [0,top]
    oplot, [0,1000], [0,1000]

    pcicpcp = 19.35*pci+8.78
    slope = 19.35
    int   =  8.78
    err   = mean(abs(pcicpcp - cpcp))
    xyouts, 5.0, top*0.85, 'Slope : '+string(slope,format='(f5.2)')
    xyouts, 5.0, top*0.70, 'Int : '+string(int,format='(f5.2)')
    xyouts, top*0.95, 5.0, 'Err : '+string(err,format='(f5.2)'), align=1.0

    fit, pci, cpcp, slope, int, err, iter = 5000
    pcicpcp_new = pci*slope + int

    print, slope, int, err

    dx = pos(2) - pos(0)
    pos(0) = pos(2) + xspace
    pos(2) = pos(0) + dx
    xposr = pos

    plot, cpcp, pcicpcp_new, pos = pos, /noerase, $
          xtitle = 'AMIE CPCP (kV)', $
          psym = 3, xrange = [0,top], yrange = [0,top]
    oplot, [0,1000], [0,1000]

    xyouts, 5.0, top*0.85, 'Slope : '+string(slope,format='(f5.2)')
    xyouts, 5.0, top*0.70, 'Int : '+string(int,format='(f5.2)')
    xyouts, top*0.95, 5.0, 'Err : '+string(err,format='(f5.2)'), align=1.0
    loc = where(cpcp lt 150.0 and pci lt 50)
    cc = c_correlate(cpcp(loc),pci(loc),0)
    xyouts, top*0.95, 25.0, 'CC : '+string(cc,format='(f4.2)'), align=1.0


    ; ------------------------------------------------------------

    pos(3) = pos(1) - space
    pos(1) = pos(3) - ysize
    pos(0) = xposl(0)
    pos(2) = xposl(2)

    plot, area, pciarea, pos = pos, /noerase, $
          xtitle = 'AMIE Area', ytitle = 'PCI Area', $
          psym = 3, xrange = [0,areatop], yrange = [0,areatop]
    oplot, [0,1000], [0,1000]

    top = areatop

    error  = mean(abs(pciarea-area))
    xyouts, top*0.95, 1.5, 'Err : '+string(error,format='(f5.2)'), align=1.0

    x2 = -1.0
    x1 = 10.0
    x0 = 2.0
    pfit, pci, area, x2, x1, x0, error, iter = 5000
    pciarea_new = x2*pci^2 + x1*pci + x0

    pos(0) = xposr(0)
    pos(2) = xposr(2)

    plot, area, pciarea_new, pos = pos, /noerase, $
          xtitle = 'AMIE Area', $
          psym = 3, xrange = [0,areatop], yrange = [0,areatop]
    oplot, [0,1000], [0,1000]

    xyouts, 1.5, top*0.85, 'S2 : '+string(x2,format='(e9.2)')
    xyouts, 1.5, top*0.70, 'S1 : '+string(x1,format='(f5.2)')
    xyouts, 1.5, top*0.55, 'Int : '+string(x0,format='(f5.2)')
    xyouts, top*0.95, 1.5, 'Err : '+string(error,format='(f5.2)'), align=1.0

    cc = c_correlate(area,pci,0)
    xyouts, top*0.95, 8.5, 'CC : '+string(cc,format='(f4.2)'), align=1.0

    ; ------------------------------------------------------------

    pos(3) = pos(1) - space
    pos(1) = pos(3) - ysize
    pos(0) = xposl(0)
    pos(2) = xposl(2)

    plot, e, pcie, pos = pos, /noerase, $
          xtitle = 'AMIE E (mV/m)', ytitle = 'PCI E (mV/m)', $
          psym = 3, xrange = [0,etop], yrange = [0,etop]
    oplot, [0,1000], [0,1000]

    x2 = -0.17 * 1.7
    x1 = 4.5   * 1.7
    x0 = 9.0   * 1.7
    pcie = x2*pci^2 + x1*pci + x0

    error  = mean(abs(pcie-e))

    top = etop

    xyouts, 5.0, top*0.85, 'S2 : '+string(x2,format='(e9.2)')
    xyouts, 5.0, top*0.70, 'S1 : '+string(x1,format='(f5.2)')
    xyouts, 5.0, top*0.55, 'Int : '+string(x0,format='(f5.2)')
    xyouts, top*0.95, 5.0, 'Err : '+string(error,format='(f5.2)'), align=1.0

    pfit, pci, e, x2, x1, x0, error, iter = 5000
    pcie_new = x2*pci^2 + x1*pci + x0

    pos(0) = xposr(0)
    pos(2) = xposr(2)

    plot, e, pcie_new, pos = pos, /noerase, $
          xtitle = 'AMIE E (mV/m)', $
          psym = 3, xrange = [0,etop], yrange = [0,etop]
    oplot, [0,1000], [0,1000]

    xyouts, 5.0, top*0.85, 'S2 : '+string(x2,format='(e9.2)')
    xyouts, 5.0, top*0.70, 'S1 : '+string(x1,format='(f5.2)')
    xyouts, 5.0, top*0.55, 'Int : '+string(x0,format='(f5.2)')
    xyouts, top*0.95, 5.0, 'Err : '+string(error,format='(f5.2)'), align=1.0

    loc = where(e lt 150.0 and pci lt 50)
    cc = c_correlate(e(loc),pci(loc),0)
    xyouts, top*0.95, 20.0, 'CC : '+string(cc,format='(f4.2)'), align=1.0

    closedevice

  endif

endif

end
