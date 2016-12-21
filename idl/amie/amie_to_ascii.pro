
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

initial_guess = findfile('-t b*')
initial_guess = initial_guess(0)
if strlen(initial_guess) eq 0 then initial_guess='b970101'

amie_file = ask('AMIE binary file name',initial_guess)
amie_file = findfile(amie_file)

nfiles = n_elements(amie_file)

for iFile = 0, nFiles-1 do begin

  cAmieFile = amie_file(iFile)

  print, "Reading file ",cAmieFile

  read_amie_binary, cAmieFile, data, lats, mlts, time, fields, 		$
  		    imf, ae, dst, hp, cpcp
  
  ;--------------------------------------------------------------
  ; figure out grid:
  ;--------------------------------------------------------------
  
  nt = n_elements(time)
  nlats = n_elements(lats)
  nmlts = n_elements(mlts)
  
  if (iFile eq 0) then begin
  
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
  
  endif
  
  data_1 = reform(data(*,type_1,*,*))
  field_1 = strcompress(fields(type_1))
  
  if (strpos(mklower(fields(type_1)),'potential') gt -1) then 	$
    units = '(kV)' 
  if (strpos(mklower(fields(type_1)),'cond')  gt -1) then 		$
    units = '(mhos)' 
  if (strpos(mklower(fields(type_1)),'electric field') gt -1) then 	$
    units = '(mV/m)' 
  if (strpos(mklower(fields(type_1)),'current') gt -1) then 		$
    units = '(A/m!E2!N)' 
  if (strpos(mklower(fields(type_1)),'joule') gt -1) then 		$
    units = '(mJ/m!E2!N)' 
  if (strpos(mklower(fields(type_1)),'energy flux') gt -1) then 	$
    units = '(W/m2)' 
  if (strpos(mklower(fields(type_1)),'mean energy') gt -1) then 	$
    units = '(keV)' 
  
  outfile = cAmieFile+'.'+strmid(fields(type_1),0,4)+'.dat'
  
  openw,1,outfile
  
  printf,1,''
  printf,1,'#FIELD'
  printf,1,strcompress(fields(type_1)),' in ',units
  
  printf,1,''
  printf,1,'#NUMERICAL VALUES'
  printf,1,format='(i5,a)',nt,' Number of Times in File'
  printf,1,format='(i5,a)',nlats,' Number of Latitudes'
  printf,1,format='(i5,a)',nmlts,' Number of MLTS'
  
  printf,1,''
  printf,1,'#LATITUDES'
  for i=0,nlats-1 do printf,1,format = '(f8.2)',lats(i)
  
  printf,1,''
  printf,1,'#MLTS'
  for i=0,nmlts-1 do printf,1,format = '(f8.2)',mlts(i)
  
  if mean(abs(data_1(0,*,*))) lt 1 then begin
    format = '('+tostr(nmlts)+'e10.2)'
  endif else begin
    format = '('+tostr(nmlts)+'f8.2)'
  endelse
  
  printf,1,''
  printf,1,'#FORMAT'
  printf,1,format
  
  for i=0,nt-1 do begin
  
    rtime = time(i)
    c_r_to_a, itime, rtime
    printf,1,''
    printf,1,'#TIME'
    printf,1,format='(6i5)',itime
  
    printf,1,''
    printf,1,'#DATA'
  
    for j = 0, nlats-1 do $
      printf,1,format = format, data_1(i,*,j)
  
  endfor
  
  close,1

endfor

end




