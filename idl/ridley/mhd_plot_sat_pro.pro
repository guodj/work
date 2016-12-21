
pro mhd_plot_sat_pro, psfile, logfilename, datatype, $
                      cdffilename, delay, alldata, endtimein, datasets, $
                      rmstime = rmstime

swap = 1

if (n_elements(rmstime) eq 0) then rmstime = -1.0

; psfile = ask("ps filename", "idl.ps")

setdevice, psfile, "p", 4, 0.95

; if (n_elements(logfilename) eq 0) then begin
;     filelist = findfile("*.sat")
;     if (strlen(filelist(0)) gt 0) then logfilename = filelist(0) $
;     else logfilename = "idontknow.sat"
; endif 
; logfilename = ask("MHD satellite output file",logfilename)

filelist = findfile(logfilename)

if (strpos(logfilename," ") gt 0) then begin
    filelist(0) = strmid(logfilename,0,strpos(logfilename," "))
    filelist(1) = strmid(logfilename,strpos(logfilename," "), strlen(logfilename))
endif

; print, logfilename
; print, 'list: ',filelist
; print, n_elements(filelist)

; print, " 1. Magnetometer Data"
; print, " 2. Velocity Data"
; print, " 3. Number Density and Temperature Data"
; datatype = fix(ask("data type to plot"," 1"))

getlogpro, filelist(0), nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2

if (n_elements(filelist) gt 1) then begin
    wlogb = wlog
    getlogpro, filelist(1), nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2
    multi = 1
endif else multi = 0

if (n_elements(cdffilename) eq 0) then begin
    filelist = findfile("*.cdf")
    if (strlen(filelist(0)) gt 0) then cdffilename = filelist(0) $
    else cdffilename = "idontknow.cdf"
endif

; cdffilename = ask("MHD satellite output file",cdffilename)

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

;        print, nrows

        bvar = "idontknow"

        for i=0,nvars-1 do begin
            result = cdf_varinq(id, i, zvar=isZvar)
;            print, i,'. ',result.name


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
                if (strpos(result.name,'THERMAL_SPD') gt -1) then begin
                    Tvar = result.name
                endif
                if (strpos(result.name,'Np') gt -1) then begin
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
            cdf_varget, id, Tvar, T_tmp, rec_count=nrows

            if (swap) then np_tmp = swap_endian(np_tmp)
            if (swap) then T_tmp = swap_endian(T_tmp)

; Change units of t_tmp here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            if (iFile eq 0) then begin
                cdfbfield = fltarr(2,nrows)
                cdfbfield(0,*) = np_tmp(0,*) 
                cdfbfield(1,*) = T_tmp(0,*) 
            endif else begin
                tmp = fltarr(2,nrows_total + nrows)
                tmp(*,0:nrows_total-1) = cdfbfield
                tmp(0,nrows_total:nrows_total + nrows-1) = np_tmp(0,*)
                tmp(1,nrows_total:nrows_total + nrows-1) = T_tmp(0,*)
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

    if (strpos(cdffilename,"fgm") gt 0) then begin

        read_themis_fgm, cdffilename, nptsCdf, cdftime, cdfbfield

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

endelse

;------------------------------------------------------
; Get delay
;------------------------------------------------------

; if (n_elements(delay) eq 0) then delay = " 0.0"
; delay = float(ask("delay in seconds",string(delay)))

npts = n_elements(wlog(*,0))

time = dblarr(npts)

for i=0L,npts-1 do begin
    itime = reform(fix(wlog(i,1:6)))
    c_a_to_r, itime, rtime
    time(i) = rtime + delay
endfor

if (multi) then begin

    nptsb = n_elements(wlogb(*,0))

    timeb = dblarr(nptsb)

    for i=0L,nptsb-1 do begin
        itime = reform(fix(wlogb(i,1:6)))
        c_a_to_r, itime, rtime
        timeb(i) = rtime + delay
    endfor

endif

;-----------------------------------------------------------
; Determine whether to plot all real data or all MHD output
;-----------------------------------------------------------

; if cdffilename ne "idontknow.cdf" then begin
;
;    alldata = $
;      ask("whether you want all times or just the MHD times plotted (a/m)","a")
;    if (strpos(mklower(alldata),"a") eq 0) then alldata = 1 else alldata = 0
;
; endif else alldata = 0

x = reform(wlog(*,8))
y = reform(wlog(*,9))
z = reform(wlog(*,10))

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
    if (multi) then begin
        stime = min([min(timeb),stime])
        etime = max([max(timeb),etime])
    endif
endif else begin
    stime = min(time)
    etime = max(time)
    if (multi) then begin
        stime = min([min(timeb),stime])
        etime = max([max(timeb),etime])
    endif
endelse

if (endtimein gt 0.0) then endtime = endtimein else $
  endtime = tostrf((etime-stime)/3600.0)
; endtime = float(ask('end time to plot (hours)',endtime))
etime   = stime + endtime * 3600.0

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

tmptime = time - stime
if (multi) then tmptimeb = timeb - stime

if (datatype eq 1) then vars = [15,16,17]
if (datatype eq 2) then vars = [12,13,14]
if (datatype eq 3) then begin
    vars = [11,18]
    ; Convert to temperature
    wlog(*,18) = (wlog(*,18)*1.0e-9)/(wlog(*,11)*1.0e6*1.3807e-23)/11000.0
    wlognames(18) = "T"
    if (multi) then $
      wlogb(*,18) = (wlogb(*,18)*1.0e-9)/(wlogb(*,11)*1.0e6*1.3807e-23)/11000.0
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


get_position, ppp, space, sizes, 0, pos, /rect
pos(1) = pos(1) - 0.1
pos(3) = pos(3) - 0.1
pos(0) = pos(0) + 0.075

plot, [0,1], [0,1], pos = pos, /noerase, /nodata, xstyle = 5, ystyle = 5

if (multi) then oplot, [0.0, 0.1], [0.1, 0.1], linestyle = 2, color = 10
oplot, [0.30, 0.40], [0.1, 0.1]
if (multi) then oplot, [0.65, 0.75], [0.1, 0.1], color = 240

l = strpos(datasets,";")
if (l eq -1) then l = strlen(datasets)
ds1 = strmid(datasets,0,l)
datasets = strmid(datasets,l+1,strlen(datasets))

l = strpos(datasets,";")
ds2 = strmid(datasets,0,l)
datasets = strmid(datasets,l+1,strlen(datasets))

xyouts, 0.29, 0.07, ds1, align = 1.0
xyouts, 0.41, 0.07, ds2
if (multi) then xyouts, 0.76, 0.07, datasets

for nv = 0,nvars-1 do begin

;    if nv eq 0 then title = logfilename else title = ""
    title = ""
    ivar = vars(nv)

    get_position, ppp, space, sizes, nv+1, pos, /rect

    pos(0) = pos(0) + 0.075
    pos(1) = pos(1) - yspace
    pos(3) = pos(3) - yspace

    if nv lt nvars-1 then begin
      xtn = strarr(60)+' '
      xt  = ' '
    endif else begin
      xtn = xtickname
      xt  = xtitle
    endelse

    if (datatype eq 3) then begin
        wlog(*,ivar) = alog10(wlog(*,ivar))
        wlognames(ivar) = "log("+wlognames(ivar)+")"
        if cdffilename ne "idontknow.cdf" then $
          cdfbfield(nv,*) = alog10(cdfbfield(nv,*))
        if (multi) then wlogb(*,ivar) = alog10(wlogb(*,ivar))
    endif

    if cdffilename ne "idontknow.cdf" then begin

        loc = where(cdftime ge stime and cdftime le etime and $
                    reform(abs(cdfbfield(nv,*))) lt 1000.0)

print, mm(cdftime)
print, stime, etime
print, mm(reform(abs(cdfbfield(nv,*))))

        yrange = [min([min(wlog(*,ivar)),min(cdfbfield(nv,loc))]), $
                  max([max(wlog(*,ivar)),max(cdfbfield(nv,loc))])]

    endif else begin

        yrange = [min(wlog(*,ivar)),max(wlog(*,ivar))]
        
    endelse

    if (multi) then $
      yrange = [min([min(wlogb(*,ivar)),min(yrange)]), $
                max([max(wlogb(*,ivar)),max(yrange)])]
        
    if (datatype eq 1) then begin
        if yrange(0) lt -500 then yrange(0) = -500.0
        if yrange(1) gt  500 then yrange(1) =  500.0
    endif

    plot, tmptime, wlog(*,ivar), xstyle=1,		$
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

    if (multi) then $
      oplot,tmptimeb, wlogb(*,ivar), thick = thick, color = 240

   oplot, [tmptime(1*npts/4)], [wlog(1*npts/4,ivar)], psym = 4, symsize=2
   oplot, [tmptime(2*npts/4)], [wlog(2*npts/4,ivar)], psym = 2, symsize=2
   oplot, [tmptime(3*npts/4)], [wlog(3*npts/4,ivar)], psym = 5, symsize=2

   if cdffilename ne "idontknow.cdf" then begin
       oplot, cdftime(loc)-stime, cdfbfield(nv,loc), $
         linestyle = 2, thick = thick, color = 10

       t1 = cdftime(loc)-stime
       t2 = tmptime

       d1 = reform(cdfbfield(nv,loc))
       d2 = reform(wlog(*,ivar))

       if multi then begin
           t3 = tmptimeb
           d3 = reform(wlogb(*,ivar))
       endif

       if (rmstime gt 0) then etime = rmstime*3600.0 $
       else begin
           etimefinal = min([max(t1),max(t2)])
           if multi then etimefinal = min([etimefinal, max(t3)])
       endelse

       smallesttime = min([max(t1),max(t2)])
       if multi then smallesttime = min([smallesttime, max(t3)])

       l = where(t1 lt smallesttime)
       t1 = t1(l)
       d1 = d1(l)
       l = where(t2 lt smallesttime)
       t2 = t2(l)
       d2 = d2(l)
       if (multi) then begin
           t3 = t3(l)
           d3 = d3(l)
       endif

       data_mhd = interpol(d2, t2, t1)

       rms  = sqrt(mean((d1-data_mhd)^2))
       rms2 = sqrt(mean((d1)^2))
       rmss = strcompress("nRMS : "+string(rms/rms2,format = "(f5.3)"))
       xyouts, pos(2)+0.01, pos(3), rmss, orient = -90, /norm

       if (multi) then begin
           data_mhd = interpol(d3, t3, t1)
           rms  = sqrt(mean((d1-data_mhd)^2))
           rmss = strcompress(string(rms/rms2,format = "(f5.3)"))
           xyouts, pos(2)+0.01, pos(1), rmss, alignment = 1, orient = -90, /norm, color = 240
       endif
   endif

endfor

closedevice


end

