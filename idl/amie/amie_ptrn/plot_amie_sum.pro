;+
; NAME:
;       plot_amie_sum.pro
;
; PURPOSE:
;       Takes a list of AMIE summary files and creates postscript plots of the result. AMIE
; 		summary files are created by make_amie_summary.
; EXPLANATION:
;       Read in the AMIE _sum data files then produce comparison plots.
;
; CALLING SEQUENCE:
;
;       This program is meant to be .run and then interviews the user.
;
; OPTIONAL INPUTS:
;
;       None
;
; OPTIONAL INPUT KEYWORDS:
;
;		None
;
; OUTPUTS:
;       Postscript images of sumarizing the difference between AMIE_sum1 and AMIE_sum2.
;
; EXAMPLE:
;		IDL> .run plot_amie_sum
;		% Compiled module: $MAIN$.
;		Enter files to plot [*_sum] :
;		Enter ps file name [summary.ps] :
;		Reading file b19980501.all_sum
;		Reading file b19980502.all_sum
;		Reading file b19980504.all_sum
;		Enter how many hours to average [1] : 1
;
; AUTHOR AND MODIFICATIONS:
;
;      Steve Trameroriginal plot_amie_sum
;      E.A. Kihn modifications 02/20/04


; Number of times in each file
num_samples = 1440L

cd,current=pwd
targetDir = ask('Data dir:', pwd)

print, 'First File Type: '
print, '1. All'
print, '2. Missing Sect'
print, '3. Only Weimer'
print, '4, RT no Weimer'
print, '5, RT Only'
print, '6, Statbox'
print, '7, Single'

;Assumes the data( AMIE runs) has extensions as described.
datatype1= fix(ask('Choice ?','1'))

  if datatype1 eq 1 then begin
    extension1 = 'all'
  endif
  if datatype1 eq 2 then begin
    extension1 = 'missect'
  endif
  if datatype1 eq 3 then begin
	extension1 = 'onlyweimer'
  endif
  if datatype1 eq 4 then begin
	extension1 = 'rtnoweimer'
  endif
  if datatype1 eq 5 then begin
	extension1 = 'rtonly'
  endif
  if datatype1 eq 6 then begin
	extension1 = 'statbox'
  endif
  if datatype1 eq 7 then begin
	extension1 = 'single'
  endif

print, 'Second File Type: '
print, '1. All'
print, '2. Missing Sect'
print, '3. Only Weimer'
print, '4, RT no Weimer'
print, '5, RT Only'
print, '6, Statbox'
print, '7, Single'

datatype2= fix(ask('Choice ?','5'))

  if datatype2 eq 1 then begin
    extension2 = 'all'
  endif
  if datatype2 eq 2 then begin
    extension2 = 'missect'
  endif
  if datatype2 eq 3 then begin
	extension2 = 'onlyweimer'
  endif
  if datatype2 eq 4 then begin
	extension2 = 'rtnoweimer'
  endif
  if datatype2 eq 5 then begin
	extension2 = 'rtonly'
  endif
  if datatype2 eq 6 then begin
	extension2 = 'statbox'
  endif
  if datatype2 eq 7 then begin
	extension2 = 'single'
  endif



;Find all the files with the right extension, can be trouble if weird extensions like .all_sum
summaryfilelist_1 = findfile(targetDir + '\*.' + extension1 +'_sum')
nfiles1 = n_elements(amiefilelist_1)


;Get the second set of extensions
summaryfilelist_2 = findfile(targetDir + '\*.' + extension2 + '_sum')
nfiles2 = n_elements(amiefilelist_2)


; Need a routine to sort and match them.

summaryfilelist_1 = summaryfilelist_1(sort(summaryfilelist_1))
summaryfilelist_2 = summaryfilelist_2(sort(summaryfilelist_2))


psfile = ask('ps file name','summary.ps')


nfiles = n_elements(summaryfilelist_1)
if (nfiles gt 0 and strlen(summaryfilelist_1(0)) gt 0) then begin
; There are 13 data elements in the _sum file like:
; CPCP        Ae        Au        Al        Bx        By        Bz       Vsw  Comp Dst  Comp HPI  Comp SJH   A (m^2)      Emax
  itime = intarr(6)
  itime2= intarr(6)
  tmp   = fltarr(13)
  tmp2  = fltarr(13)
  data  = fltarr(num_samples*nfiles,13)
  data2 = fltarr(num_samples*nfiles,13)
  time  = dblarr(num_samples*nfiles)
  line  = ''

  for ifile=0L,nfiles-1 do begin

    print, "Reading file ",summaryfilelist_1(ifile)
    openr,1,summaryfilelist_1(ifile)
    readf,1,line
    print, "Reading file ",summaryfilelist_2(ifile)
    openr,2,summaryfilelist_2(ifile)
    readf,2,line ;takes of the one line header.

    for i=0,num_samples-1 do begin
      readf,1,itime,tmp, FORMAT='(I5,5I3,13E10.2)'
      readf,2,itime2,tmp2, FORMAT='(I5,5I3,13E10.2)'
      c_a_to_r, itime, rtime
      data(ifile*num_samples+i,*) = tmp
      data2(ifile*num_samples+i,*) = tmp2
      time(ifile*num_samples+i)   = rtime
    endfor
    close,1
    close,2

  endfor

  stime = min(time)
  etime = max(time)

  ;get_imf, stime, etime, imfdata, dt=60.0

  time_axis, stime, etime, btr, etr, $
           xtickname, xtitle, xtickv, xminor, xtickn

  n_min_avg = fix(ask('How many minutes to avg?','1'))

  if n_min_avg gt 0 then begin


    n_times = fix (size(time,/n_elements)/n_min_avg)
    timeh = dblarr(n_times)
    datah = fltarr(n_times,13)
    datah2 = fltarr(n_times,13)
    window = double(n_min_avg * 60) ; We're working in milisec's
		t = 0.0
    for i=0,n_times-1 do begin
      t = stime + i*window
      t_loc = where(time ge t and time lt t+window, count)

      if count gt 0 then begin
          timeh(i) = t
          datah(i,*) = (TOTAL (data(t_loc,*),1))/count
          datah2(i,*) = (TOTAL (data2(t_loc,*),1))/count
      endif else begin
          timeh(i) = t
          datah(i,*) = datah(i-1,*)
      endelse
    endfor

    timeh = timeh - stime

    cpcp = reform(datah(*,0))
    cpcp2 = reform (datah2(*,0))
    by = reform(datah(*,5))
    bz = reform(datah(*,6))
    al = reform(datah(*,3))
    al2 = reform(datah2(*,3))
    area = reform(datah(*,11))
    area2 = reform(datah2(*,11))

    bt = sqrt(by(*,0)^2+bz(*,0)^2)
    clock = acos(bz(*,0)/(bt+0.01))

    vx = reform(abs(datah(*,7,0)))
    loc = where(vx lt 200.0,count)
    if count gt 0 then vx(loc) = 400.0

    ekl = bt*vx*sin(clock/2.0)*1000.0*1.0e-9*1000.0

    setdevice, psfile,'p',4
    plotdumb
	;PPP= Plots per page
    ppp = 5
    space = 0.01
    pos_space, ppp, space, sizes, ny=ppp

    get_position, ppp, space, sizes, ppp-1, pos, /rect
    pos(0) = pos(0) + 0.075

    plot, timeh, cpcp, xtickn = xtickname, xtitle = xtitle, $
        xtickv = xtickv, xminor = xminor, xticks = xtickn, $
        xrange = [btr,etr], pos = pos, /noerase, xstyle =1, $
	ytitle = 'CPCP (kV)', yrange = [0.0,max(cpcp)]

    oplot, timeh, cpcp2, linestyle = 1
    ;oplot, timeh, cpcp(*,2), linestyle = 1
    ;oplot, timeh, cpcp(*,3), linestyle = 2

    get_position, ppp, space, sizes, 0, pos, /rect
    pos(0) = pos(0) + 0.075

    plot, timeh, bz, xtickn = strarr(10)+' ', $
        xtickv = xtickv, xminor = xminor, xticks = xtickn, $
        xrange = [btr,etr], pos = pos, /noerase, xstyle =1, $
	ytitle = 'Bz (nT)'

    ;oplot, time, data(*,4), linestyle = 1

;    oplot, timeh, by(*,1), linestyle = 1
;    oplot, timeh, by(*,2), linestyle = 1
;    oplot, [btr,etr],[0.0,0.0], linestyle =2
;    oplot, timeh, by(*,3), linestyle = 2

    get_position, ppp, space, sizes, 1, pos, /rect
    pos(0) = pos(0) + 0.075

    plot, timeh, al, xtickn = strarr(10)+' ', $
        xtickv = xtickv, xminor = xminor, xticks = xtickn, $
        xrange = [btr,etr], pos = pos, /noerase, xstyle =1, $
	ytitle = 'Al (nT)'

    oplot, timeh, al2, linestyle = 1

;    oplot, timeh, bz(*,1), linestyle = 1
;    oplot, timeh, bz(*,2), linestyle = 1
    oplot, [btr,etr],[0.0,0.0], linestyle =2

    get_position, ppp, space, sizes, 2, pos, /rect
    pos(0) = pos(0) + 0.075

    plot, timeh, area, xtickn = strarr(10)+' ', $
        xtickv = xtickv, xminor = xminor, xticks = xtickn, $
        xrange = [btr,etr], pos = pos, /noerase, xstyle =1, $
	ytitle = 'Area m^2'

    oplot, timeh, area2, linestyle = 1

    get_position, ppp, space, sizes, 3, pos, /rect
    pos(0) = pos(0) + 0.075

    ;plot, imfdata.time-stime, imfdata.den, xtickn = strarr(10)+' ', $
     ;   xtickv = xtickv, xminor = xminor, xticks = xtickn, $
      ;  xrange = [btr,etr], pos = pos, /noerase, xstyle =1, $
	;ytitle = 'N (/cm3)'

    closedevice

  endif else begin

    time = time - stime

    plot, time, data(*,0), xtickn = xtickname, xtitle = xtitle, $
        xtickv = xtickv, xminor = xminor, xticks = xtickn

  endelse

endif

end
