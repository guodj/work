
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

initial_guess = findfile('-t b*')
initial_guess = initial_guess(0)
if strlen(initial_guess) eq 0 then initial_guess='b970101'

amie_file = ask('AMIE binary file name',initial_guess)

read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp

;--------------------------------------------------------------
; figure out grid:
;--------------------------------------------------------------

nt = n_elements(time)
nlats = n_elements(lats)
nmlts = n_elements(mlts)

;--------------------------------------------------------------
; Need to figure out what to plot, so list to fields to the user
;--------------------------------------------------------------

nfields = n_elements(fields)
for i=0,nfields-1 do print, tostr(i+1)+'. '+fields(i)

;--------------------------------------------------------------
; Get field to be contoured
;--------------------------------------------------------------

type_1 = fix(ask('field to output','1'))-1
if (type_1 lt 0) or (type_1 gt nfields-1) then type_1 = 0

outfile = ask('output file name',amie_file+'.'+$
                                 strmid(fields(type_1),0,4)+'.dat')

data_1 = reform(data(*,type_1,*,*))
field_1 = strcompress(fields(type_1))

if (strpos(mklower(fields(type_1)),'potential') gt -1) then 		$
    units = '(kV)' 
if (strpos(mklower(fields(type_1)),'cond')  gt -1) then 		$
    units = '(mhos)' 
if (strpos(mklower(fields(type_1)),'electric field') gt -1) then 	$
    units = '(mV/m)' 
if (strpos(mklower(fields(type_1)),'current') gt -1) then 		$
    units = '(A/m!E2!N)' 
if (strpos(mklower(fields(type_1)),'joule') gt -1) then 		$
    units = '(mJ/m!E2!N)' 
if (strpos(mklower(fields(type_1)),'energy flux') gt -1) then 		$
    units = '(W/m2)' 
if (strpos(mklower(fields(type_1)),'mean energy') gt -1) then 		$
    units = '(keV)' 

quantity = strcompress(fields(type_1))+' in '+units
data = data_1

save, file = outfile, quantity, nlats, nmlts, lats, mlts, time, data

end




