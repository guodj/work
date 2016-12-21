
CCMCoutput = 1

basedir = '/csem1/Data6/Data/cdaweb.gsfc.nasa.gov/pub/istp/'

filelist = findfile('sat*_n000000.sat')
display, filelist

if (n_elements(sat) eq 0) then sat = 'geotail'
sat = ask('sat name (e.g., cluster, goes10, polar)',sat)

psfile = ask("ps filename", sat+".ps")

setdevice, psfile, "p", 4, 0.95

filelist = findfile('sat_'+sat+'*.sat')

nfiles=n_elements(filelist)
for iFile = 0, nFiles-1 do begin
    logfilename = filelist(iFile)
    spawn, 'wc '+logfilename,result
    if (fix(result) gt 0) then begin
        getlogpro, logfilename, $
          nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2
        if (iFile eq 0) then begin
            wlogtotal = wlog
        endif else wlogtotal = [wlogtotal, wlog]
    endif
endfor

;------------------------------------------------------
; Get delay
;------------------------------------------------------

if (n_elements(delay) eq 0) then delay = " 0.0"
delay = float(ask("delay in seconds",string(delay)))

npts = n_elements(wlogtotal(*,0))

time = dblarr(npts)

for i=0L,npts-1 do begin
    itime = reform(fix(wlogtotal(i,1:6)))
    c_a_to_r, itime, rtime
    time(i) = rtime + delay
endfor

print, " 1. Magnetometer Data"
print, " 2. Velocity Data"
print, " 3. Number Density and Temperature Data"
datatype = fix(ask("data type to plot"," 1"))

;getlogpro, logfilename, nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2

c_r_to_a, itime, min(time)
c_a_to_ymd, itime, ymd
sday = itime(2)

c_r_to_a, itime, max(time)
eday = itime(2)
c_a_to_ymd, itime, ymd2

print, sday, eday

days = 0
if (eday-sday ne 0) then days = 1

if (strpos(sat,'goes6') eq 0) or (strpos(sat,'goes06') eq 0) then begin
    filehead = 'goes/6_mag/'
    filemid  = '/g6_k0_mag_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'goes7') eq 0) or (strpos(sat,'goes07') eq 0) then begin
    filehead = 'goes/7_mag/'
    filemid  = '/g7_k0_mag_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'goes8') eq 0) or (strpos(sat,'goes08') eq 0) then begin
    filehead = 'goes/8_mag/'
    filemid  = '/g8_k0_mag_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'goes9') eq 0) or (strpos(sat,'goes09') eq 0) then begin
    filehead = 'goes/9_mag/'
    filemid  = '/g9_k0_mag_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'goes10') eq 0) then begin
    filehead = 'goes/0_mag/'
    filemid  = '/g0_k0_mag_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'goes11') eq 0) then begin
    filehead = 'goes/11_mag/'
    filemid  = '/goes11_k0_mag_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'goes12') eq 0) then begin
    filehead = 'goes/12_mag/'
    filemid  = '/goes12_k0_mag_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'polar') eq 0) then begin
    filehead = 'polar/mfe/'
    filemid  = '/po_k0_mfe_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'wind') eq 0) then begin
    filehead = 'wind/mfi_h0/'
    filemid  = '/wi_h0_mfi_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'geotail') eq 0) then begin
    filehead = 'geotail/mgf/'
    filemid  = '/ge_k0_mgf_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'1994-084') eq 0) then begin
    filehead = 'lanl/94_mpa/'
    filemid  = '/l4_k0_mpa_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'LANL-97A') eq 0) then begin
    filehead = 'lanl/97_mpa/'
    filemid  = '/l7_k0_mpa_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'lanl01') eq 0) then begin
    filehead = 'lanl/01a_mpa/'
    filemid  = '/a1_k0_mpa_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'lanl02') eq 0) then begin
    filehead = 'lanl/02a_mpa/'
    filemid  = '/a2_k0_mpa_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'lanl90') eq 0) then begin
    filehead = 'lanl/90_mpa/'
    filemid  = '/l0_k0_mpa_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'lanl94') eq 0) then begin
    filehead = 'lanl/94_mpa/'
    filemid  = '/l4_k0_mpa_'
    fileend  = '_???.cdf'
endif

if (strpos(sat,'lanl97') eq 0) then begin
    filehead = 'lanl/97_mpa/'
    filemid  = '/l7_k0_mpa_'
    fileend  = '_???.cdf'
endif

file1 = basedir+filehead+strmid(ymd,0,4)+filemid+ymd+fileend
file2 = ''
if (days) then $
  file2 = basedir+filehead+strmid(ymd2,0,4)+filemid+ymd2+fileend

cdffilename = file1+' '+file2

if (n_elements(cdffilename) eq 0) then begin
    filelist = findfile("*.cdf")
    if (strlen(filelist(0)) gt 0) then cdffilename = filelist(0) $
    else cdffilename = "idontknow.cdf"
endif

cdffilename = ask("MHD satellite output file",cdffilename)

if (cdffilename ne "idontknow.cdf" and $
    strpos(cdffilename, "cdf") gt 0) then begin

    list = findfile(cdffilename)
    if (strlen(list(0)) eq 0) then begin
        print, "Can not file file : ",cdffilename
        stop
    endif

    nFiles = n_elements(list)
    nrows_total = 0

    for iFile = 0, nFiles-1 do begin

        cdffilename = list(iFile)

        print, "Reading file : ", cdffilename

        id = cdf_open(cdffilename)
        re = cdf_inquire(id)
        nvars = re.nvars  
        nrows = re.maxrec

        isZVar = 0
        if (nrows lt 1) then begin
            cdf_control, id, variable=0, /zvariable,get_var_info=v
            nrows = v.maxrec
            if (nrows eq 0) then nrows = v.maxrecs
            nVars = re.nzvars
            isZVar = 1
        endif

        bvar = "idontknow"

        for i=0,nvars-1 do begin
            result = cdf_varinq(id, i, zvar=isZvar)

            if (datatype eq 1) then begin
                if (result.name eq "B_GSM_c" or $
                    result.name eq "BGSMc" or $
                    result.name eq "B_xyz_gse__C1_PP_FGM" or $
                    result.name eq "B_xyz_gse__C2_PP_FGM" or $
                    result.name eq "B_xyz_gse__C3_PP_FGM" or $
                    result.name eq "B_xyz_gse__C4_PP_FGM" or $
                    result.name eq "BGSM" or $
                    result.name eq "B_GSM" or $
                    result.name eq "IB_vector") then bvar = result.name
            endif

            if (datatype eq 2) then begin
                if (strpos(result.name,'V') gt -1 and $
                     strpos(result.name,'GSM') gt -1 and $
                     strpos(result.name,'label') eq -1) then begin
                    bvar = result.name
                endif
            endif

            if (datatype eq 3) then begin
                if (strpos(result.name,'THERMAL_SPD') gt -1 or $
                    strpos(result.name,'temp_hip') gt -1) then begin
                   Tvar = result.name
                endif
                if (strpos(result.name,'Np') gt -1 or $
                    strpos(result.name,'dens_hip') gt -1) then begin
                    npvar = result.name
                endif
            endif

        endfor

        timevar = "Epoch"

        IsFound = 0
        for i=0,nvars-1 do begin
            result = cdf_varinq(id, i, zvar=isZvar)
            print, i,'. ',result.name
            if (strpos(result.name,"Epoch") gt -1 $
                and not IsFound) then begin
                timevar = result.name
                IsFound = 1
            endif
        endfor

        cdf_varget, id, timevar, cdftime_tmp, rec_count=nrows
        siz = size(cdftime_tmp)
        if (siz(0) eq 2) then cdftime_tmp = reform(cdftime_tmp(0,*))
        if (siz(0) eq 3) then cdftime_tmp = reform(cdftime_tmp(0,0,*))
        swap = 0
        print, max(cdftime_tmp)
        if (max(cdftime_tmp) lt 1.0e10) then swap = 1

        if (swap) then cdftime_tmp = swap_endian(cdftime_tmp)

        if (iFile eq 0) then cdftime = cdftime_tmp $
        else cdftime = [cdftime,cdftime_tmp]

        if (datatype eq 1 or datatype eq 2) then begin
            cdf_varget, id, bvar, cdfbfield_tmp, rec_count=nrows
            if (swap) then cdfbfield_tmp = swap_endian(cdfbfield_tmp)
            if (iFile eq 0) then cdfbfield = cdfbfield_tmp $
            else begin
                tmp = fltarr(3,nrows_total + nrows)
                tmp(*,0:nrows_total-1) = cdfbfield
                tmp(*,nrows_total:nrows_total + nrows-1) = cdfbfield_tmp
                cdfbfield = tmp
            endelse
        endif else begin

            cdf_varget, id, npvar, np_tmp, rec_count=nrows
            if (swap) then np_tmp = swap_endian(np_tmp)
            s = size(np_tmp)
            cdf_varget, id, Tvar, T_tmp, rec_count=nrows
            if (swap) then T_tmp = swap_endian(T_tmp)
            if (s(0) eq 3) then begin
               np_tmp = reform(np_tmp(0,0,*))
               t_tmp = reform(t_tmp(0,0,*))
            endif
; Change units of t_tmp here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (iFile eq 0) then begin
                cdfbfield = fltarr(2,nrows)
                cdfbfield(0,*) = np_tmp
                cdfbfield(1,*) = T_tmp
            endif else begin
                tmp = fltarr(2,nrows_total + nrows)
                tmp(*,0:nrows_total-1) = cdfbfield
                tmp(0,nrows_total:nrows_total + nrows-1) = np_tmp(*)
                tmp(1,nrows_total:nrows_total + nrows-1) = T_tmp(*)
                cdfbfield = tmp
            endelse


        endelse

        nrows_total = nrows_total + nrows

    endfor

    nrows = nrows_total

    if (bvar eq "IB_vector") then cdfbfield = cdfbfield/10.0

    if (strpos(mklower(bvar),"gse") gt -1) then gse_gsm, cdfbfield, cdftime

    for i=0L,nrows-1 do begin
        cdf_epoch, cdftime(i), year, month, day, hour, minute, second, milli,/break
        itime = [year, month, day, hour, minute, second]
        c_a_to_r, itime, rtime
        cdftime(i) = rtime
    endfor

endif else begin

    if (strpos(cdffilename,"cdf") lt 0) then begin

        getlogpro, cdffilename, nlogfile, nVarsLog, wlognamesCdf, $
          wlogCdf, wlogCdf0, wlogCdf1, wlogCdf2

        nptsCdf = n_elements(wlogCdf(*,0))

        cdftime = dblarr(nptsCdf)

        loc = where(strpos(mklower(wlognamesCdf),"year") eq 0)
        if (loc(0) eq -1) then begin
            print, "Can not find year variable in the log file!!!"
            print, wlognamesCdf
            stop
        endif

        for i=0L,nptsCdf-1 do begin
            itime = reform(fix(wlogCdf(i,loc(0):loc(0)+5)))
            c_a_to_r, itime, rtime
            cdftime(i) = rtime
        endfor

        if (datatype eq 1) then begin
            loc = where(strpos(mklower(wlognamesCdf),"bx") eq 0)
            if (loc(0) eq -1) then begin
                print, "Can not find bx variable in the log file!!!"
                print, wlognamesCdf
                stop
            endif
        endif

        if (datatype eq 2) then begin
            loc = where(strpos(mklower(wlognamesCdf),"vx") eq 0)
            if (loc(0) eq -1) then begin
                print, "Can not find vx variable in the log file!!!"
                print, wlognamesCdf
                stop
            endif
        endif

        if (datatype eq 3) then begin
            loc = where(strpos(mklower(wlognamesCdf),"np") eq 0)
            if (loc(0) eq -1) then begin
                print, "Can not find Np variable in the log file!!!"
                print, wlognamesCdf
                stop
            endif
            loc2 = where(strpos(mklower(wlognamesCdf),"tp") eq 0)
            if (loc2(0) eq -1) then begin
                print, "Can not find tp variable in the log file!!!"
                print, wlognamesCdf
                stop
            endif
        endif

        if (datatype eq 1 or datatype eq 2) then begin
            cdfbfield = fltarr(3,nptsCdf)
            for i=0,2 do cdfbfield(i,*) = reform(wlogCdf(*,loc(0)+i))
        endif else begin
            cdfbfield = fltarr(2,nptsCdf)
            cdfbfield(0,*) = reform(wlogCdf(*,loc(0)))
            cdfbfield(1,*) = reform(wlogCdf(*,loc2(0)))
        endelse

    endif

endelse

;-----------------------------------------------------------
; Determine whether to plot all real data or all MHD output
;-----------------------------------------------------------

if cdffilename ne "idontknow.cdf" then begin

    alldata = $
      ask("whether you want all times or just the MHD times plotted (a/m)","a")
    if (strpos(mklower(alldata),"a") eq 0) then alldata = 1 else alldata = 0

endif else alldata = 0

x = reform(wlogtotal(*,8))
y = reform(wlogtotal(*,9))
z = reform(wlogtotal(*,10))

xmini = round(min([min(x),-10.0])/10.0)*10.0
xmaxi = (floor(max([max(x), 10.0])/10.0)+1)*10.0

ymaxi = (floor(max([max(abs(y)),10.0])/10.0)+1)*10.0
zmaxi = (floor(max([max(abs(z)),10.0])/10.0)+1)*10.0

xrange = [ xmaxi, xmini]
yrange = [ ymaxi,-ymaxi]
zrange = [-zmaxi, zmaxi]

if alldata then begin
    stime = min([min(time),min(cdftime)])
    etime = max([max(time),max(cdftime)])
endif else begin
    stime = min(time)
    etime = max(time)
endelse

endtime = tostrf((etime-stime)/3600.0)
endtime = float(ask('end time to plot (hours)',endtime))
etime   = stime + endtime * 3600.0

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

tmptime = time - stime

if (datatype eq 1) then vars = [15,16,17]
if (datatype eq 2) then vars = [12,13,14]
if (datatype eq 3) then begin
    vars = [11,18]
    ; Convert to temperature
    wlogtotal(*,18) = (wlogtotal(*,18)*1.0e-9)/(wlogtotal(*,11)*1.0e6*1.3807e-23)/11000.0
    wlognames(18) = "T"
endif

if (CCMCOutput) then begin

   c_r_to_a, itime, stime
   c_a_to_ymd, itime, ymd
   if (strpos(sat,'*') ge 0 or strpos(sat,'?') ge 0) then begin
      p = strpos(psfile,'.ps')
      s = strmid(psfile,0,p)
      CCMCOutputFile  = s+'_'+ymd+'_CCMC.dat'
      CCMCOutputFile2 = s+'_'+ymd+'_CCMC_mp.dat'
   endif else begin
      CCMCOutputFile  = sat+'_'+ymd+'_CCMC.dat'
      CCMCOutputFile2 = sat+'_'+ymd+'_CCMC_mp.dat'
   endelse
   close,2,3
   openw,2, CCMCOutputFile
   openw,3, CCMCOutputFile2

   itime(4) = 0
   itime(5) = 0
   c_a_to_r, itime, stime_ccmc
   c_r_to_a, itime, etime
   if (itime(4) gt 0) then itime(3) = itime(3) + 1
   itime(4) = 0
   itime(5) = 0
   c_a_to_r, itime, etime_ccmc

   dt = 60.0

   nTimesCCMC = (etime_ccmc - stime_ccmc)/dt + 1

   xyzCCMC = fltarr(3,nTimesCCMC)
   datCCMC = fltarr(3,nTimesCCMC)
   timCCMC = dblarr(nTimesCCMC)
   staCCMC = strarr(nTimesCCMC)+'I'

   for iT = 0, nTimesCCMC-1 do begin
      t = stime_ccmc + dt * iT
      ; iterpolate
      l = where(time gt t, c)

      if (c gt 0) then begin
         if (l(0) eq 0) then begin
            l = 1
            r = 1.0
         endif else begin
            l = l(0)
            r = (time(l) - t) / (time(l) - time(l-1))
         endelse
      endif else begin
         l = n_elements(time)-1
         r = 0.0
      endelse

      tOut = (1-r) * time(l) + r * time(l-1)
      c_r_to_a, itime, t
      epoch = 0.0
      cdf_epoch,epoch,$
                itime(0),itime(1),itime(2),itime(3),itime(4),itime(5),$
                /compute
      xyzCCMC(0,iT) = (1-r) * x(l) + r * x(l-1)
      xyzCCMC(1,iT) = (1-r) * y(l) + r * y(l-1)
      xyzCCMC(2,iT) = (1-r) * z(l) + r * z(l-1)
      timCCMC(iT) = epoch

      if (wlogtotal(l,24) ne 3 or wlogtotal(l-1,24) ne 3) then $
         staCCMC(iT) = 'O'

      if (datatype eq 1) then begin
         datCCMC(0,iT) = (1-r) * wlogtotal(l,15) + r * wlogtotal(l-1,15)
         datCCMC(1,iT) = (1-r) * wlogtotal(l,16) + r * wlogtotal(l-1,16)
         datCCMC(2,iT) = (1-r) * wlogtotal(l,17) + r * wlogtotal(l-1,17)
      endif

      if (datatype eq 3) then begin
         datCCMC(0,iT) = (1-r) * wlogtotal(l,11) + r * wlogtotal(l-1,11)
         datCCMC(1,iT) = (1-r) * wlogtotal(l,18) + r * wlogtotal(l-1,18)
      endif

   endfor

   gsm_gse, xyzCCMC, timCCMC

   if (datatype eq 1) then begin
      gsm_gse, datCCMC, timCCMC
   endif

   for iT = 0, nTimesCCMC-1 do begin
      t = stime_ccmc + dt * iT
      c_r_to_a, itime, t
      if (datatype eq 1) then $
         printf, 2, itime, xyzCCMC(*,iT), datCCMC(*,iT), $
                 format='(i5,5i3,6f10.2)'
      if (datatype eq 3) then $
         printf, 2, itime, xyzCCMC(*,iT), datCCMC(0:1,iT), $
                 format='(i5,5i3,5f10.2)'
      printf, 3, itime, xyzCCMC(*,iT), staCCMC(iT), format='(i5,5i3,3f10.2,a2)'
   endfor

   close,2,3

endif

nvars = n_elements(vars)
ppp = nvars+2
space = 0.01
pos_space, ppp, space, sizes, ny = ppp

makect, 'mid'

if strpos(!d.name,'X') gt -1 then thick = 1 else thick = 3

plotdumb

get_position, ppp, space, sizes, 0, posm

yspace = (posm(3) - posm(1))/2.0

posl = posm
posl(2) = posl(2)-posl(0)+space/2 + 0.075
posl(0) = space/2 + 0.075

prx = posl(2) - posl(0)
pmx = (posl(2)+posl(0))/2.0
pry = posl(3) - posl(1)
pmy = (posl(3)+posl(1))/2.0

if (xrange(0)-xrange(1) gt yrange(0)-yrange(1)) then begin
  pry = pry * (yrange(0)-yrange(1))/(xrange(0)-xrange(1))
  posl(3) = pmy + pry/2.0
  posl(1) = posl(3) - pry
endif else begin
  prx = prx * (xrange(0)-xrange(1))/(yrange(0)-yrange(1))
  posl(2) = pmx + prx/2.0
  posl(0) = posl(2) - prx
endelse

plot, xrange, yrange, /nodata, pos = posl, /noerase, $
  xtitle = 'X-GSM (Re)', ytitle = 'Y-GSM (Re)', $
  xstyle =1, ystyle = 1, $
  xtickname = ['',' ','',' ','',' ','',' ','',' ']

oplot, x, y, thick = thick, color = 10
oplot, [x(1*npts/4)], [y(1*npts/4)], psym = 4
oplot, [x(2*npts/4)], [y(2*npts/4)], psym = 2
oplot, [x(3*npts/4)], [y(3*npts/4)], psym = 5
oplot, [-1000.0,1000.0], [0.0,0.0], linestyle = 1
oplot, [0.0,0.0], [-1000.0,1000.0], linestyle = 1
theta = findgen(19)/18.0*2.0*!pi
polyfill, 3.0*cos(theta), 3.0*sin(theta)

posr = posm
posr(0) = 1.0 - space/2.0 - (posr(2)-posr(0))
posr(2) = 1.0 - space/2

prx = posr(2) - posr(0)
pmx = (posr(2)+posr(0))/2.0
pry = posr(3) - posr(1)
pmy = (posr(3)+posr(1))/2.0

if (xrange(0)-xrange(1) gt zrange(1)-zrange(0)) then begin
  pry = pry * (zrange(1)-zrange(0))/(xrange(0)-xrange(1))
  posr(3) = pmy + pry/2.0
  posr(1) = posr(3) - pry
endif else begin
  prx = prx * (xrange(0)-xrange(1))/(zrange(1)-zrange(0))
  posr(2) = pmx + prx/2.0
  posr(0) = posr(2) - prx
endelse

plot, xrange, zrange, /nodata, pos = posr, /noerase, $
  xtitle = 'X-GSM (Re)', ytitle = 'Z-GSM (Re)', $
  xstyle =1, ystyle = 1, $
  xtickname = ['',' ','',' ','',' ','',' ','',' ']

oplot, x, z, thick = thick, color = 10
oplot, [x(1*npts/4)], [z(1*npts/4)], psym = 4
oplot, [x(2*npts/4)], [z(2*npts/4)], psym = 2
oplot, [x(3*npts/4)], [z(3*npts/4)], psym = 5
oplot, [-1000.0,1000.0], [0.0,0.0], linestyle = 1
oplot, [0.0,0.0], [-1000.0,1000.0], linestyle = 1
polyfill, 3.0*cos(theta), 3.0*sin(theta)

for nv = 0,nvars-1 do begin

    if nv eq 0 then title = logfilename else title = ""
    ivar = vars(nv)

    get_position, ppp, space, sizes, nv+1, pos, /rect

    pos(0) = pos(0) + 0.075
    pos(1) = pos(1) - yspace
    pos(3) = pos(3) - yspace
    pos(2) = pos(2) - 0.03

    if nv lt nvars-1 then begin
      xtn = strarr(60)+' '
      xt  = ' '
    endif else begin
      xtn = xtickname
      xt  = xtitle
    endelse

    if (datatype eq 3) then begin
        wlogtotal(*,ivar) = alog10(wlogtotal(*,ivar))
        wlognames(ivar) = "log("+wlognames(ivar)+")"
        if cdffilename ne "idontknow.cdf" then $
          cdfbfield(nv,*) = alog10(cdfbfield(nv,*))
    endif

    if cdffilename ne "idontknow.cdf" then begin

        loc = where(cdftime ge stime and cdftime le etime and $
                    reform(abs(cdfbfield(nv,*))) lt 1000.0)
        
        yrange = [min([min(wlogtotal(*,ivar)),min(cdfbfield(nv,loc))]), $
                  max([max(wlogtotal(*,ivar)),max(cdfbfield(nv,loc))])]
        
    endif else begin

        yrange = [min(wlogtotal(*,ivar)),max(wlogtotal(*,ivar))]
        
    endelse

    if (datatype eq 1) then begin
        if (yrange(0) lt -300) then yrange(0) = -300.0
        if (yrange(1) gt  300) then yrange(1) =  300.0
    endif

    plot, tmptime, wlogtotal(*,ivar), xstyle=1,		$
      ytitle = wlognames(ivar),		$
      xtickname = xtn,			$
      xtitle = xt,			$
      xtickv = xtickv,			$
      xminor = xminor,			$
      xticks = xtickn,   $
      xrange = [btr, etr], $
      pos = pos, /noerase, $
      thick = thick, yrange = yrange, title = title

    oplot, [btr, etr], [0.0,0.0], linestyle = 1

   oplot, [tmptime(1*npts/4)], [wlogtotal(1*npts/4,ivar)], psym = 4, symsize=2
   oplot, [tmptime(2*npts/4)], [wlogtotal(2*npts/4,ivar)], psym = 2, symsize=2
   oplot, [tmptime(3*npts/4)], [wlogtotal(3*npts/4,ivar)], psym = 5, symsize=2

   if cdffilename ne "idontknow.cdf" then begin
       oplot, cdftime(loc)-stime, cdfbfield(nv,loc), $
         linestyle = 2, thick = thick, color = 10

       t1 = cdftime(loc)-stime
       t2 = tmptime
       d2 = reform(wlogtotal(*,ivar))

       data_mhd = interpolate_mine(t1, d2, t2)
       data_sat = reform(cdfbfield(nv,loc))

       rms  = sqrt(mean((data_sat-data_mhd)^2))
       rms2 = sqrt(mean((data_sat)^2))

       rmss = strcompress("RMS : "+string(rms,format = "(f8.2)")+ $
         " / "+string(rms2,format = "(f8.2)"))
       xyouts, pos(2)+0.01, pos(3), rmss, orient = -90, /norm

       rmss = strcompress("nRMS : "+string(rms/rms2,format = "(f8.3)"))
       xyouts, pos(2)+0.03, pos(3), rmss, orient = -90, /norm

   endif

endfor

closedevice


end

