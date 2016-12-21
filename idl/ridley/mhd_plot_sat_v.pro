
psfile = ask("ps filename", "idl.ps")

setdevice, psfile, "p", 4, 0.95

if (n_elements(logfilename) eq 0) then begin
    filelist = findfile("*.sat")
    if (strlen(filelist(0)) gt 0) then logfilename = filelist(0) $
    else logfilename = "idontknow.sat"
endif 
logfilename = ask("MHD satellite output file",logfilename)

getlogpro, logfilename, nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2

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

;        print, nrows

        bvar = "idontknow"

        for i=0,nvars-1 do begin
            result = cdf_varinq(id, i, zvar=isZvar)
;            print, i,'. ',result.name
            tmp = mklower(result)
            if ((strpos(tmp,'v') gt -1 and $
                 strpos(tmp,'gsm') gt -1) and $ 
                Vp eq -1) then begin
                bvar = result.name
            endif
        endfor

;        print, bvar

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

        if (iFile eq 0) then cdftime = cdftime_tmp $
        else cdftime = [cdftime,cdftime_tmp]

        cdf_varget, id, bvar, cdfbfield_tmp, rec_count=nrows

        if (iFile eq 0) then cdfbfield = cdfbfield_tmp $
        else begin
            tmp = fltarr(3,nrows_total + nrows)
            tmp(*,0:nrows_total-1) = cdfbfield
            tmp(*,nrows_total:nrows_total + nrows-1) = cdfbfield_tmp
            cdfbfield = tmp
        endelse

        nrows_total = nrows_total + nrows

    endfor

    nrows = nrows_total

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

        loc = where(strpos(mklower(wlognamesCdf),"vx") eq 0)
        if (loc(0) eq -1) then begin
            print, "Can not find vx variable in the log file!!!"
            print, wlognamesCdf
            stop
        endif

        cdfbfield = fltarr(3,nptsCdf)
        for i=0,2 do cdfbfield(i,*) = reform(wlogCdf(*,loc(0)+i))

    endif

endelse

;------------------------------------------------------
; Get delay
;------------------------------------------------------

if (n_elements(delay) eq 0) then delay = " 0.0"
delay = float(ask("delay in seconds",string(delay)))

npts = n_elements(wlog(*,0))

time = dblarr(npts)

for i=0L,npts-1 do begin
    itime = reform(fix(wlog(i,1:6)))
    c_a_to_r, itime, rtime
    time(i) = rtime + delay
endfor

;-----------------------------------------------------------
; Determine whether to plot all real data or all MHD output
;-----------------------------------------------------------

if cdffilename ne "idontknow.cdf" then begin

    alldata = $
      ask("whether you want all times or just the MHD times plotted (a/m)","a")
    if (strpos(mklower(alldata),"a") eq 0) then alldata = 1 else alldata = 0

endif else alldata = 0

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
endif else begin
    stime = min(time)
    etime = max(time)
endelse

endtime = tostrf((etime-stime)/3600.0)
endtime = float(ask('end time to plot (hours)',endtime))
etime   = stime + endtime * 3600.0

time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

tmptime = time - stime

vars = [12,13,14]

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
oplot, cos(theta), sin(theta)

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
oplot, cos(theta), sin(theta)

for nv = 0,nvars-1 do begin

    if nv eq 0 then title = logfilename else title = ""
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

    if cdffilename ne "idontknow.cdf" then begin

        loc = where(cdftime ge stime and cdftime le etime and $
                    reform(abs(cdfbfield(nv,*))) lt 1000.0)
        
        yrange = [min([min(wlog(*,ivar)),min(cdfbfield(nv,loc))]), $
                  max([max(wlog(*,ivar)),max(cdfbfield(nv,loc))])]
        
    endif else begin

        yrange = [min(wlog(*,ivar)),max(wlog(*,ivar))]
        
    endelse

    if yrange(0) lt -500 then yrange(0) = -500.0
    if yrange(1) gt  500 then yrange(1) =  500.0

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

   oplot, [tmptime(1*npts/4)], [wlog(1*npts/4,ivar)], psym = 4, symsize=2
   oplot, [tmptime(2*npts/4)], [wlog(2*npts/4,ivar)], psym = 2, symsize=2
   oplot, [tmptime(3*npts/4)], [wlog(3*npts/4,ivar)], psym = 5, symsize=2

   if cdffilename ne "idontknow.cdf" then begin
       oplot, cdftime(loc)-stime, cdfbfield(nv,loc), $
         linestyle = 2, thick = thick, color = 10

       t1 = cdftime(loc)-stime
       t2 = tmptime
       d2 = reform(wlog(*,ivar))
       data_mhd = interpol(d2, t2, t1)
       data_sat = reform(cdfbfield(nv,loc))

       rms  = sqrt(mean((data_sat-data_mhd)^2))
       rms2 = sqrt(mean((data_sat)^2))
       rmss = strcompress("RMS : "+string(rms,format = "(f8.2)")+ $
         " / "+string(rms2,format = "(f8.2)"))
       xyouts, pos(2)+0.01, pos(3), rmss, orient = -90, /norm

   endif

endfor

closedevice


end

