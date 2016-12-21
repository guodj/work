
GetNewData = 1

if (n_elements(start) eq 0) then start='????'

old_start = start

start = ask('starting characters of the satellite',start)

filelist_new = findfile(start+"*.*ALL")
nfiles_new = n_elements(filelist_new)
if (nfiles_new eq 1) then begin
    filelist_new = findfile(start+"*.dat")
    nfiles_new = n_elements(filelist_new)
endif

if (nfiles_new eq 1) then begin
    filelist_new = findfile(start+"*.bin")
    nfiles_new = n_elements(filelist_new)
endif

if n_elements(nfiles) gt 0 then begin
    if (nfiles_new eq nfiles) then default = 'n' else default='y'
    if (strpos(old_start,start) eq -1) then default = 'n'
    GetNewData = mklower(strmid(ask('whether to reread data',default),0,1))
    if (GetNewData eq 'n') then GetNewData = 0 else GetNewData = 1
endif

if (GetNewData) then begin

    thermo_readsat, filelist_new, data, time, nTimes, Vars, nAlts, nSats, Files
    nFiles = n_elements(filelist_new)

endif

if (nSats eq 1) then begin

    nPts = nTimes

    Alts = reform(data(0,0:nPts-1,2,0:nalts-1))/1000.0
    Lons = reform(data(0,0:nPts-1,0,0)) * 180.0 / !pi
    Lats = reform(data(0,0:nPts-1,1,0)) * 180.0 / !pi

    c_r_to_a, itime, time(0)
    itime(3:5) = 0
    c_a_to_r, itime, basetime
    hour = (time/3600.0 mod 24.0) + fix((time-basetime)/(24.0*3600.0))*24.0
    localtime = (Lons/15.0 + hour) mod 24.0
    
    angle = 23.0 * !dtor * $
      sin((jday(itime(0),itime(1),itime(2)) - jday(itime(0),3,21))*2*!pi/365.0)
    angle = 0
    sza =  acos(sin(angle)*sin(Lats*!dtor) + $
                cos(angle)*cos(Lats*!dtor) * $ 
                cos(!pi*(LocalTime-12.0)/12.0))

    rho = reform(data(0,0:nPts-1,3,0:nalts-1))
    t   = reform(data(0,0:nPts-1,14,0:nalts-1))
    o   = reform(data(0,0:nPts-1,4,0:nalts-1))
    n   = reform(data(0,0:nPts-1,6,0:nalts-1)) +  $
          reform(data(0,0:nPts-1,9,0:nalts-1)) + $
          reform(data(0,0:nPts-1,10,0:nalts-1))
    no  = reform(data(0,0:nPts-1,8,0:nalts-1))

    nop  = reform(data(0,0:nPts-1,28,0:nalts-1))
    op   = reform(data(0,0:nPts-1,24,0:nalts-1))
    np   = reform(data(0,0:nPts-1,27,0:nalts-1))
    o2p  = reform(data(0,0:nPts-1,25,0:nalts-1))

endif

mint = min(time)
c_r_to_a, imint, mint
c_a_to_ymd, imint, ymd
c_a_to_s, imint, date
date = strmid(date,0,9)

time0 = imint
time = time - time(0)

save, rho, t, o, n, no, nop, op, np, o2p, alts, lats, lons, $
      time, time0, file='qb50_gitm_'+ymd+'.save'

setdevice, 'neutral_densities_'+ymd+'.ps','p',5
plotdumb

plot_oi, rho, alts, psym = 3, yrange = [80,400], xrange = [1e-13,1e-5], $
         xtitle = 'Rho (kg/m3)', ytitle = 'Alt (km)', thick = 4, $
         pos = [0.05, 0.6, 0.45, 1.0], /noerase, ystyle = 1

plot_oi, no, alts, psym = 3, yrange = [80,400], xrange = [1e10,1e20], $
         xtitle = '[NO] (/m3)', ytitle = 'Alt (km)', thick = 4, $
         pos = [0.05, 0.1, 0.45, 0.5], /noerase, ystyle = 1

plot_oi, o, alts, psym = 3, yrange = [80,400], xrange = [1e10,1e20], $
         xtitle = '[O] (kg/m3)', ytitle = ' ', thick = 4, $
         pos = [0.51, 0.6, 0.91, 1.0], /noerase, ystyle = 1

plot_oi, n, alts, psym = 3, yrange = [80,400], xrange = [1e10,1e20], $
         xtitle = '[N] (kg/m3)', ytitle = ' ', thick = 4, $
         pos = [0.51, 0.1, 0.91, 0.5], /noerase, ystyle = 1

xyouts, 0.5, 1.01, date, charsize = 1.3, alignment = 0.5, /norm

closedevice

setdevice, 'ion_densities_'+ymd+'.ps','p',5
plotdumb

plot_oi, op, alts, psym = 3, yrange = [80,400], xrange = [1e5,1e13], $
         xtitle = '[O+] (/m3)', ytitle = 'Alt (km)', thick = 4, $
         pos = [0.05, 0.6, 0.45, 1.0], /noerase, ystyle =1

plot_oi, np, alts, psym = 3, yrange = [80,400], xrange = [1e5,1e13], $
         xtitle = '[N+] (/m3)', ytitle = 'Alt (km)', thick = 4, $
         pos = [0.05, 0.1, 0.45, 0.5], /noerase, ystyle =1

plot_oi, nop, alts, psym = 3, yrange = [80,400], xrange = [1e5,1e13], $
         xtitle = '[NO+] (/m3)', ytitle = ' ', thick = 4, $
         pos = [0.51, 0.6, 0.91, 1.0], /noerase, ystyle =1

plot_oi, o2p, alts, psym = 3, yrange = [80,400], xrange = [1e5,1e13], $
         xtitle = '[O2+] (/m3)', ytitle = ' ', thick = 4, $
         pos = [0.51, 0.1, 0.91, 0.5], /noerase, ystyle =1

xyouts, 0.5, 1.01, date, charsize = 1.3, alignment = 0.5, /norm

closedevice

end
