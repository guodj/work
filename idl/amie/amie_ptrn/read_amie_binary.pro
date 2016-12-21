pro read_amie_binary, amiefile, data, lats, mlts, ut, fields, imf,	$
                      ae, dst, hp, pot, date = date,			$
                      ltpos = ltpos, lnpos = lnpos, 			$
                      plotapot = plotapot, speed = speed, by = by, 	$
                      bz = bz, field = field

  if n_elements(plotapot) eq 1 then plotapot = 1 else plotapot = 0
  if n_elements(field) eq 0 then field = -1

  openr,1,amiefile, /f77

  nlats = 0L
  nmlts = 0L
  ntimes = 0L

  n = 0L
  iyr_tmp = 0L
  imo_tmp = 0L
  ida_tmp = 0L
  ihr_tmp = 0L
  imi_tmp = 0L
  nfields = 0L

  swv_tmp = 0.0
  bx_tmp  = 0.0
  by_tmp  = 0.0
  bz_tmp  = 0.0
  aei_tmp = 0.0
  ae_tmp  = 0.0
  au_tmp  = 0.0
  al_tmp  = 0.0
  dsti_tmp= 0.0
  dst_tmp = 0.0
  hpi_tmp = 0.0
  sjh_tmp = 0.0
  pot_tmp = 0.0

  readu,1,nlats,nmlts,ntimes

  clats = fltarr(nlats)
  mlts = fltarr(nmlts)

  readu,1,clats
  readu,1,mlts

  readu,1,nfields

  ut  = dblarr(ntimes)
  imf = fltarr(ntimes,4)
  ae  = fltarr(ntimes,4)
  dst = fltarr(ntimes,2)
  hp  = fltarr(ntimes,2)
  pot = fltarr(ntimes)
  data  = fltarr(ntimes,nfields,nmlts,nlats)
  dummy = fltarr(nlats,nmlts)
  dummy = fltarr(nmlts,nlats)
  fields = strarr(nfields)

  tmp = bytarr(30)

  for i=0,nfields-1 do begin
    readu,1,tmp
    fields(i) = string(tmp)
  endfor

  for i=0,ntimes-1 do begin

    readu,1,n,iyr_tmp,imo_tmp,ida_tmp,ihr_tmp,imi_tmp
    itime = [fix(iyr_tmp),imo_tmp,ida_tmp,ihr_tmp,imi_tmp,0]

    c_a_to_r, itime, rtime
    ut(i) = rtime

    readu,1,swv_tmp,bx_tmp,by_tmp,bz_tmp,aei_tmp,ae_tmp,au_tmp,al_tmp,    $
          dsti_tmp,dst_tmp,hpi_tmp,sjh_tmp,pot_tmp

    imf(i,0) = bx_tmp
    imf(i,1) = by_tmp
    imf(i,2) = bz_tmp
    imf(i,3) = swv_tmp

    ae(i,0)  = ae_tmp 
    ae(i,1)  = au_tmp 
    ae(i,2)  = al_tmp 
    ae(i,3)  = aei_tmp 

    dst(i,0) = dst_tmp
    dst(i,1) = dsti_tmp

    hp(i,0)  = hpi_tmp
    hp(i,1)  = sjh_tmp

    pot(i)   = pot_tmp

    for j=0,nfields-1 do begin
      readu,1,dummy
;      data(i,j,*,*) = rotate(dummy,1)
      data(i,j,*,*) = dummy
      if (strpos(mklower(fields(j)),'potential') gt -1) then 		$
        data(i,j,*,*) = dummy/1000.0
    endfor

  endfor

  close,1

  lats = 90.0 - clats

  if plotapot then begin

    if (field lt 0) then begin
      for i=0,nfields-1 do print, tostr(i+1)+'. '+fields(i)
      type = fix(ask('field to plot','1'))-1
      if (type lt 0) or (type gt nfields-1) then type = 0
    endif else type = field

    data = reform(data(*,type,*,*))
    fields = fields(type)

    ltpos = fltarr(nmlts,nlats)
    lnpos = fltarr(nmlts,nlats)
    for i=0,nmlts-1 do ltpos(i,*) = lats
    for j=0,nlats-1 do lnpos(*,j) = mlts*360.0/24.0

    date = strarr(ntimes)
    time = strarr(ntimes)

    for i=0,ntimes-1 do begin
      c_r_to_a, itime, ut(i)
      c_a_to_s, itime, stime
      time(i) = strmid(stime,10,5)
      syear = fix(strmid(stime,7,2))
      if syear lt 65 then syear = syear + 2000 else syear = syear + 1900
      date(i) = strmid(stime,3,1)+mklower(strmid(stime,4,2))+' '+	$
                strmid(stime,0,2)+', '+tostr(syear)
    endfor

    ut = time
    speed = reform(imf(*,3))
    by    = reform(imf(*,1))
    bz    = reform(imf(*,2))

    lats = nlats
    mlts = nmlts

  endif else begin

    if field gt -1 then data = reform(data(*,field,*,*))

  endelse

  return

end


