
; you have to do this before you can run the procedure:
; openu,lu,"cgm.1993"
; a = fstat(lu)
; fsize = a.size

pro read_mag_bin, fin, fsize, glat, glon, mlat, mmid, bh, bz, b_total

  len_of_row = 24

; 11
  glat11 = 2.0*float(floor(glat/2))
  if glon lt 0.0 then glon = glon + 360.0
  if glon ge 360.0 then glon = glon - 360.0
  glon11 = 2.0*float(floor(glon/2))
  dlat = long((90.0-glat11)/2)
  dlon = long(glon11/2)
  total = dlat*(360/2 + 1) + dlon
  total11 = long(len_of_row)*total

; 21
  glon21 = glon11
  glat21 = 2.0*float(ceil(glat/2))
  dlat = long((90.0-glat21)/2)
  total = dlat*(360/2 + 1) + dlon
  total21 = long(len_of_row)*total

; 22
  glat22 = glat21
  glon22 = 2.0*float(ceil(glon/2))
  dlon = long(glon22/2)
  total = dlat*(360/2 + 1) + dlon
  total22 = long(len_of_row)*total

; 12
  glon12 = glon22
  glat12 = 2.0*float(floor(glat/2))
  dlat = long((90.0-glat11)/2)
  total = dlat*(360/2 + 1) + dlon
  total12 = long(len_of_row)*total

  g11 = 0
  g12 = 0
  g21 = 0
  g22 = 0

  if (total11 le fsize-len_of_row) then begin
    point_lun,fin,total11
    readu,1,glat11,glon11,mlat11,bh11,bz11,mmid11
  endif else g11 = 1

  if (total12 le fsize-len_of_row) then begin
    point_lun,fin,total12
    readu,1,glat12,glon12,mlat12,bh12,bz12,mmid12
  endif else g12 = 1

  if (total21 le fsize-len_of_row) then begin
    point_lun,fin,total21
    readu,1,glat21,glon21,mlat21,bh21,bz21,mmid21
  endif else g21 = 1

  if (total22 le fsize-len_of_row) then begin
    point_lun,fin,total22
    readu,1,glat22,glon22,mlat22,bh22,bz22,mmid22
  endif else g22 = 1

  if (g11+g12+g21+g22 gt 0) then begin
    mlat = -999.0
    mmid = -999.0
    bh   = 0.0
    bz   = 0.0
    print, 'bad point'
  endif else begin

    w1a = 2.0 - (glat - glat11)
    w1o = 2.0 - (glon - glon11)
    w2a = 2.0 - w1a
    w2o = 2.0 - w1o

    mlat = (mlat11*w1a*w1o+mlat12*w1a*w2o+mlat21*w2a*w1o+mlat22*w2a*w2o)/ $
           (w1a*w1o+w1a*w2o+w2a*w1o+w2a*w2o)
    bh = (bh11*w1a*w1o+bh12*w1a*w2o+bh21*w2a*w1o+bh22*w2a*w2o)/ $
         (w1a*w1o+w1a*w2o+w2a*w1o+w2a*w2o)
    bz = (bz11*w1a*w1o+bz12*w1a*w2o+bz21*w2a*w1o+bz22*w2a*w2o)/ $
         (w1a*w1o+w1a*w2o+w2a*w1o+w2a*w2o)

    if (abs(mmid11) gt 24.0) then w1a = 0.0
    if (abs(mmid21) gt 24.0) then w2a = 0.0

    check = 0
    if (mmid11 lt 2.0) then begin
      if (mmid21 gt 22.0) then mmid21 = mmid21-24.0
      if (mmid12 gt 22.0) then mmid12 = mmid12-24.0
      if (mmid22 gt 22.0) then mmid22 = mmid22-24.0
      check = 1
    endif
    if (mmid21 lt 2.0) and (check ne 1) then begin
      if (mmid11 gt 22.0) then mmid11 = mmid11-24.0
      if (mmid12 gt 22.0) then mmid12 = mmid12-24.0
      if (mmid22 gt 22.0) then mmid22 = mmid22-24.0
      check = 1
    endif
    if (mmid12 lt 2.0) and (check ne 1) then begin
      if (mmid21 gt 22.0) then mmid21 = mmid21-24.0
      if (mmid11 gt 22.0) then mmid11 = mmid11-24.0
      if (mmid22 gt 22.0) then mmid22 = mmid22-24.0
      check = 1
    endif
    if (mmid22 lt 2.0) and (check ne 1) then begin
      if (mmid21 gt 22.0) then mmid21 = mmid21-24.0
      if (mmid11 gt 22.0) then mmid11 = mmid11-24.0
      if (mmid12 gt 22.0) then mmid12 = mmid12-24.0
      check = 1
    endif

    mmid = (mmid11*w1a*w1o+mmid12*w1a*w2o+mmid21*w2a*w1o+mmid22*w2a*w2o)/ $
           (w1a*w1o+w1a*w2o+w2a*w1o+w2a*w2o)

    mmid = (mmid + 24.0) mod 24.0

  endelse

  b_total = sqrt(bh^2.0 + bz^2.0)

  return

end
