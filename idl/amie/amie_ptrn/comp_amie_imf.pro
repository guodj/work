
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

amie_file = ask('AMIE binary file name','b910603.all')
psfile = ask('ps file',amie_file+'.ps')

read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp

;--------------------------------------------------------------
; Need to figure out what to plot, so list to fields to the user
;--------------------------------------------------------------

nfields = n_elements(fields)
for i=0,nfields-1 do print, tostr(i+1)+'. '+fields(i)

;--------------------------------------------------------------
; Get field to be cross correlated
;--------------------------------------------------------------

type = fix(ask('field to plot','1'))-1
if (type lt 0) or (type gt nfields-1) then type = 0

n_start = 0
n_end = n_elements(time)-1

;--------------------------------------------------------------
; Put the contour data into data_1 array and get field name
;--------------------------------------------------------------

data_1 = reform(data(*,type,*,*))
field = strcompress(fields(type))

nlats = n_elements(lats)
nmlts = n_elements(mlts)

;--------------------------------------------------------------
; Set up device
;--------------------------------------------------------------

;setdevice, psfile, 'l', 4, 0.95

plot, data_1(*,12,5)

cc = c_correlate(data_1(*,12,5), imf(*,2), indgen(21)-10)

loc = where(abs(cc) eq max(abs(cc)))
print, 'Time shift of ',tostr((loc(0)-10)*5),' minutes is needed'
print, 'for a cross correlation of ',cc(loc(0))

plot, shift(reform(imf(*,2)),loc(0)-10), data_1(*,12,5)*1000.0, psym = 3

closedevice

end




