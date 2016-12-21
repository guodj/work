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

  c_a_to_s, itime, strtime
  date = strmid(strtime,3,3)+", "+syear

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

;    loc = where(pcicpcp lt 0.0, count)
;    if count gt 0 then pcicpcp(loc) = 0.0

    fakepci = findgen(31)-5.0
    fakepcicpcp = 19.35*fakepci+8.78

    fakepciram = -0.12*fakepci^2 + 2.5*fakepci + 11.0
    fakepcirpm = -0.20*fakepci^2 + 3.0*fakepci + 12.5

    fakepciram = fakepciram * 113.0*1000.0
    fakepcirpm = fakepcirpm * 113.0*1000.0

    fakepciarea = !pi*(0.25*fakepciram^2 + 0.25*fakepcirpm^2 + 0.5*(0.5*(fakepciram+fakepcirpm))^2)
    fakepciarea = fakepciarea/1.0e12

    fakepcie = (-0.17*fakepci^2 + 4.5*fakepci + 9.0) * 1.7

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

    setdevice, psfile,'p',4

    plotdumb

    ppp = 3
    space = 0.075

    top = 175.0

    pos_space, ppp, space, sizes, ny=ppp

    get_position, ppp, space, sizes, 0, pos, /rect
    pos(0) = pos(0)+0.075

    plot, pci, cpcp, psym = 1, yrange = [0,top], pos = pos, $
          ytitle = 'CPCP (kV)', xrange = [-2,10], /noerase, title = date

    oplot, fakepci, fakepcicpcp

    slope = 19.35
    int   =  8.78
    loc = where(cpcp lt top)
    err   = sqrt(mean((pcicpcp(loc) - cpcp(loc))^2))
    oplot, [9.1,9.75], [12,12]
    xyouts, 9.0, 10, "P="+string(slope,format='(f5.2)')+"*PCI+"+ $
      string(int,format='(f5.2)')+ " (err:"+string(err,format='(f5.2)')+")", $
      alignment = 1.0

    loc = where(cpcp lt top)
    fit, pci(loc), cpcp(loc), slope, int, err, iter = 5000

    fakepcicpcp = slope*fakepci+int
    oplot, fakepci, fakepcicpcp, linestyle = 2

    oplot, [9.1,9.75], [32,32], linestyle = 2
    xyouts, 9.0, 30, "P="+string(slope,format='(f5.2)')+"*PCI+"+ $
      string(int,format='(f5.2)')+ " (err:"+string(err,format='(f5.2)')+")", $
      alignment = 1.0

    loc = where(cpcp lt 150.0 and pci lt 50)
    cc = c_correlate(cpcp(loc),pci(loc),0)
    xyouts, -1.75, 175.0, 'CC : '+string(cc,format='(f5.2)')

;-----------------------------------------------

    top = 30.0

    get_position, ppp, space, sizes, 1, pos, /rect
    pos(0) = pos(0)+0.075
    pos(1) = pos(1)+space*0.33
    pos(3) = pos(3)+space*0.33

    plot, pci, area, psym = 1, yrange = [0,top], pos = pos, $
          ytitle = 'Area (x10!E12!N m!E2!N)', xrange = [-2,10], /noerase

    oplot, fakepci, fakepciarea

    x2 = 0.3
    x1 = 10.0
    x0 = 20.0
    efit, pci, area, x2, x1, x0, err, iter = 5000, /nosubzero

    fakepciarea = x0 - x1*exp(-x2*fakepci)
    oplot, fakepci, fakepciarea, linestyle = 2

    loc = where(area lt 50.0 and pci lt 50)
    cc = c_correlate(area(loc),pci(loc),0)
    xyouts, -1.75, 27.0, 'CC : '+string(cc,format='(f5.2)')

    oplot, [9.1,9.75], [5,5], linestyle = 2
    xyouts, 9.0, 4, "A="+string(x0,format='(f6.2)')+"-"+ $
      string(x1,format='(f6.2)')+ $
      "*exp(-"+string(x2,format='(f5.2)')+"*PCI)"+$
      " (err:"+string(err,format='(f5.2)')+")",alignment = 1.0

    err   = sqrt(mean((pciarea - area)^2))
    oplot, [9.1,9.75], [2,2]
    xyouts, 9.0, 1, "Old Fit (err:"+string(err,format='(f5.2)')+")", $
      alignment = 1.0

;-----------------------------------------------

    top = 200.0

    get_position, ppp, space, sizes, 2, pos, /rect
    pos(0) = pos(0)+0.075
    pos(1) = pos(1)+space*0.66
    pos(3) = pos(3)+space*0.66

    plot, pci, e, psym = 1, yrange = [0,top], pos = pos, $
          ytitle = 'E (mV/m)', xrange = [-2,10], /noerase, $
      xtitle = 'PCI'

    oplot, fakepci, fakepcie

    slope = 4.5   * 1.7
    int = 9.0   * 1.7
    err   = sqrt(mean((pcie - e)^2))
    loc = where(e lt top)
    fit, pci(loc), e(loc), slope,int, err, iter = 5000

    fakepcie = slope*fakepci+int
    oplot, fakepci, fakepcie, linestyle = 2

;    oplot, [5.5,6.0], [32,32], linestyle = 2
;    xyouts, 9.5, 30, "P="+string(slope,format='(f5.2)')+"*PCI+"+ $
;      string(int,format='(f5.2)')+ " (err:"+string(err,format='(f5.2)')+")", $
;      alignment = 1.0

    loc = where(e lt top and pci lt 50)
    cc = c_correlate(e(loc),pci(loc),0)
    xyouts, -1.75, top*0.90, 'CC : '+string(cc,format='(f5.2)')

    oplot, [9.1,9.75], [32,32], linestyle = 2
    xyouts, 9.0, 30, "E="+string(slope,format='(f5.2)')+"*PCI+"+ $
      string(int,format='(f5.2)')+ " (err:"+string(err,format='(f5.2)')+")", $
      alignment = 1.0

    x2    = -0.29
    slope =  7.65
    int   = 15.30
    loc = where(e lt top)
    err   = sqrt(mean((pcie(loc) - e(loc))^2))
    oplot, [9.1,9.75], [12,12]
    xyouts, 9.0, 10, "E="+string(x2,format='(f5.2)')+"*PCI!E2!N+"+$
      string(slope,format='(f5.2)')+"*PCI+"+ $
      string(int,format='(f5.2)')+ " (err:"+string(err,format='(f5.2)')+")", $
      alignment = 1.0

  endif

endif

end
