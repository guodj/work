;Program to read in Auroral Electrojet data to determine onset times
;of substorms

file = 'ae_2012_november.dat'            ;filename goes here

lun = 1
close, lun
openr, lun, file

spawn, "wc "+file, wc
nHeaders = 15
nRows = long(wc(0))-nHeaders

line = ''
for i=0, nHeaders-1 do readf, lun, line

date   = strarr(nRows)
time   = strarr(nRows)
doy    = intarr(nRows)
year   = intarr(nRows)
month  = intarr(nRows)
day    = intarr(nRows)
hour   = intarr(nRows)
minute = intarr(nRows)
second = intarr(nRows)
;ms    = strarr(nRows)
;doy   = strarr(nRows)
ae = fltarr(nRows)
au = fltarr(nRows)
al = fltarr(nRows)
ao = fltarr(nRows)

line = ''

for iRow=0L,nRows-1 do begin
   readf,1,line
   year(iRow)   = fix(strmid(line,0,4))
   month(iRow)  = fix(strmid(line,5,2))
   day(iRow)    = fix(strmid(line,8,2))
   hour(iRow)   = fix(strmid(line,11,2))
   minute(iRow) = fix(strmid(line,14,2))
   second(iRow) = fix(strmid(line,17,2))
   ae(iRow) = float(strmid(line,32,9))
   au(iRow) = float(strmid(line,41,9))
   al(iRow) = float(strmid(line,51,9))
   ao(iRow) = float(strmid(line,61,9))
endfor

close, /all

;change this path/filename to print the file where needed:
fname = '/raid2/Gitm/Runs/201211/ae_2012_november_reformat.dat'
openw, lun, fname

printf, lun, '    YEAR     ', 'MONTH    ', 'DAY     ', 'HOUR  ', 'MINUTE      ', 'AE           ', 'MLat         ' , 'MLT  '

sumtest = 0
addup = 0

print, 'writing data...'

;these values for mlat and mlt are only placeholders to be used with GITM
mlat =  65.0
mlt  =  0.0

ae = fix(ae)
al = fix(al)
au = fix(au)

for t = 0L, nRows-31 do begin
   for i = 4, 30 do begin
      addup = al(t+i)
      sumtest = sumtest + addup
   end
   
   sumtest=sumtest/26
   
   if al(t+1) - al(t) lt -15 and al(t+2) - al(t) lt -30 and al(t+3) - al(t) lt -45 and sumtest - al(t) lt -100 then begin
      printf, lun, year(t), month(t), day(t), hour(t), minute(t), second(t), ae(t), al(t), au(t)
      t=t+20
   end
   
end



my_new_time = intarr(6, n_elements(year))
itime = intarr(size(year, /n_elements))
AE_time = ulon64arr(size(itime, /n_elements))
for i = 0, size(year, /n_elements)-1 do begin
   new_itime = [year(i), month(i), day(i), hour(i), minute(i), second(i)]
   my_new_time(*,i) = new_itime
   c_a_to_r, new_itime, rtime
   AE_time(i) = rtime
endfor

stime = min(AE_time)
etime = max(AE_time)

determine_min_max, [ae, au, al], mini, maxi
time_axis, stime, etime, btr, etr, xtickname, xtitle, xtickv, xminor, xtickn

makect, 'bry'

setdevice, 'ae.ps', 'l', 5

plot, AE_time - stime, [mini*1.1, maxi*1.5], $
      xstyle = 1, xtickname = xtickname, $
      xtitle = xtitle, xtickv = xtickv, $
      xminor = xminor, xticks = xtickn, $
      /noerase, /nodata


plots, [btr, etr], [0, 0], thick = 1, color = 0
oplot, AE_time-stime, ae, color = 0
oplot, AE_time-stime, au, color = 10
oplot, AE_time-stime, al, color = 100
xyouts, 1100, 9000, 'nT', orient = 90, charsize = 1.1, /device

;legend, ['AE', 'AU', 'AL'], textcolors = [0, 10, 100], color = [0, 10, 100], $
;        /horizontal, box=0, position = [0.05, 0.03]
xyouts, etr+1000, 300, 'AE', color = 0, orient = 270
xyouts, etr+1000, 0, 'AU', color = 10, orient = 270
xyouts, etr+1000, -300, 'AL', color = 100, orient = 270
plots, [etr+1500, etr+1500], [200, 100], thick = 4, color = 0
plots, [etr+1500, etr+1500], [-100, -200], thick = 4, color = 10
plots, [etr+1500, etr+1500], [-400, -500], thick = 4, color = 100

closedevice





end

