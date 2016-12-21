pro orbitmaker_lanl_cdf, cdf_file, sat_name=sat_name, file_name=file_name

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;orbitmaker_lanl_cdf                                  DTW, Jan 2008
;
;This procedure generates an orbit file from the input cdf file.  It
;should only be used on CDAweb's cdf files of LANL spacecraft.  An
;orbit file is generated covering the time of the entire cdf file.
;
;This program uses the CXFORM library.
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

; Check arguments.
if n_elements(cdf_file) eq 0 then message, $
  'ERROR: no input file specified.  Aborting.'

cdf_file = findfile(cdf_file)
cdf_file = cdf_file(0)

if strlen(cdf_file) eq 0 then message, $
  'ERROR: cant find input file.'

if n_elements(sat_name) eq 0 then read, sat_name, $
    prompt =  "Please enter name of satellite (no spaces):"
    

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Open, parse, read and close CDF.
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

; Open cdf, collect file metadata.
id = cdf_open(cdf_file)
info = cdf_inquire(id)

; Look for variables "Epoch" and "sc_pos_geo".  
; If they are not found, we have problems.
pos_var = -1
timevar = -1
for i=0, info.nvars-1 do begin
    var_info = cdf_varinq(id, i)
    if (var_info.name eq 'Epoch')      then timevar = i
    if (var_info.name eq 'sc_pos_geo') then pos_var = i
endfor
if pos_var eq -1 then message, 'ERROR: Spacecraft position varible not found!'
if timevar eq -1 then message, 'Error: Epoch variable not found!'

; Collect data from file.
cdf_varget, id, pos_var, sc_pos, rec_count=info.maxrec
cdf_varget, id, timevar, epoch,  rec_count=info.maxrec

swap = 0
if (max(epoch) < 1.0e10) then swap = 1
if (swap) then begin
    epoch = swap_endian(epoch)
    sc_pos  = swap_endian(sc_pos)
endif

; Close CDF.
cdf_close, id

; Make sc_pos double precision.
;sc_pos = double(sc_pos)

; Fix bad geosynch radius locations.
loc = where(sc_pos[0,0,*] gt 10.0, count)
if count gt 0 then sc_pos[0,0,loc] = 6.61400

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Open orbit file, write header.
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

; Build output file name.
; The file_name keyword allows custom file names.
if n_elements(file_name) eq 0 then begin
    NameFmt = '(a8,"_",i4.4,2i2.2,a4)'
    cdf_epoch, epoch[0,0,i], yy, mm, dd, hh, mn, ss, ms, /BREAKDOWN_EPOCH
    outfile = string(sat_name,yy[0],mm[0],dd[0],'.dat',format=NameFmt)
endif else outfile = file_name

fmt = '(i5, 5(1x,i2.2), " 000", 3(1x,f7.2))'

; Open file and write header.
openw, lun, outfile, /get_lun

printf, lun, ''
printf, lun, 'Position read from file ', cdf_file
printf, lun, ''
printf, lun, '#COOR'
printf, lun, 'GSM'
printf, lun, ''
printf, lun, '#START'


;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
; Convert variables, write to file.
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2r = double(!pi)/180.0D
dog = 'bark' ;Dummy var for debugging

; FOR loop because some of these commands don't support vectors.
; Of particular note are the time routines, which, for CDFs, are
; retarded.  Also, the spacecraft location coordinates are retarded.
; Lots of down syndrome stuff going on here.
for i=0, info.maxrec-1 do begin
;for i=0, 2 do begin
    ;print, "-----"
    ;print, "sc_pos = ", sc_pos[0,*,i]
    xgeo = sc_pos[0,0,i] * cos(d2r*sc_pos[0,1,i]) * cos(d2r*sc_pos[0,2,i])
    ygeo = sc_pos[0,0,i] * cos(d2r*sc_pos[0,1,i]) * sin(d2r*sc_pos[0,2,i])
    zgeo = sc_pos[0,0,i] * sin(d2r*sc_pos[0,1,i])
    ;print, 'Geo XYZ = ', x_geo, y_geo, z_geo
    cdf_epoch, epoch[0,0,i], yy, mm, dd, hh, mn, ss, ms, /BREAKDOWN_EPOCH
;    temp_epoch = date2es(mm, dd, yy, hh, mn, ss)

    iday = jday(yy,mm,dd)
    recalc,yy,iday,hh,mn,ss
    GEOGSM,XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,1.0

;    x_gsm = CXFORM([x_geo, y_geo, z_geo],'GEO','GSM',temp_epoch)
    ;print, "GSM XYZ = ",x_gsm
    ; Pause for easy debugging.
    ;read, dog, prompt = 'Ready?  Hit enter.'

    printf, lun, yy,mm,dd,hh,mn,ss,xgsm,ygsm,zgsm, format=fmt
    ;printf, lun, yy,mm,dd,hh,mn,ss,x_geo,y_geo,z_geo, format=fmt
endfor

close, lun
free_lun, lun

end
