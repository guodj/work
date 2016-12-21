
pro read_cdf_pos, filename, xyz, time

  id = CDF_OPEN(filename)
  
  info = cdf_inquire(id)
  
  natt = info.natts
  nvar = info.nvars
  if nvar eq 0 then begin
    nvar = info.nzvars
    zvar = 1
  endif else zvar = 0
  
  vars = strarr(nvar)
  
  nrows = info.maxrec
  
  if (nrows lt 1) then begin
    cdf_control, id, variable=0, /zvariable,get_var_info=v
    nrows = v.maxrec
    if (nrows eq 0) then nrows = v.maxrecs
  endif
  
  xyz  = fltarr(3, nrows) - 1.0e32
  time = dblarr(nrows)

  epoch = -1
  pos   = -1
  IsCluster = 0
  IsGSE = 0

  for i=0,nvar-1 do begin
  
      if (zvar eq 0) then r = cdf_varinq(id,i)		$
      else r = cdf_varinq(id,i,/zvariable)
      vars(i) = r.name
      tmp = mklower(vars(i))
      if (strpos(tmp,'epoch') gt -1 and epoch lt 0) then epoch = i
      if (strpos(tmp,'pos') gt -1 and $
          (strpos(tmp,'gsm') gt -1 or strpos(tmp,'sm') gt -1) and $
          pos eq -1) then pos = i
      if (strpos(tmp,'pgsm') gt -1 and strpos(tmp,'label') lt 0) then pos = i
      if (strpos(tmp,'sc_pos_mag') gt -1 and pos lt 0) then pos = i
      if (strpos(tmp,'sc_r_xyz_gse__cl_sp_aux') gt -1 and pos lt 0) then begin
          pos = i
          IsCluster = 1
          IsGSE = 1
      endif
      if (strpos(tmp,'sc_dr1_xyz_gse__cl_sp_aux') gt -1) then pos1 = i
      if (strpos(tmp,'sc_dr2_xyz_gse__cl_sp_aux') gt -1) then pos2 = i
      if (strpos(tmp,'sc_dr3_xyz_gse__cl_sp_aux') gt -1) then pos3 = i
      if (strpos(tmp,'sc_dr4_xyz_gse__cl_sp_aux') gt -1) then pos4 = i

  endfor

  if (pos eq -1) then begin

      for i=0,nvar-1 do begin
          if (zvar eq 0) then r = cdf_varinq(id,i)		$
          else r = cdf_varinq(id,i,/zvariable)
          vars(i) = r.name
          tmp = mklower(vars(i))
          if (strpos(tmp,'pos') gt -1 and strpos(tmp,'label') lt 0) then pos = i
      endfor

  endif

  swap = 0
  if (epoch gt -1) then begin
    cdf_varget, id, vars(epoch), times, rec_count=nrows
    time = reform(times(0,*))
    if (max(time) < 1.0e10) then swap = 1
    if (swap) then time = swap_endian(time)
  endif

  if (pos gt -1) then begin

      if (IsCluster) then begin
          mid  = fltarr(3, nrows) - 1.0e32
          tmp  = fltarr(3, nrows) - 1.0e32
          xyz  = fltarr(3, nrows, 4) - 1.0e32
          cdf_varget, id, vars(pos), mid, rec_count=nrows
          if (swap) then mid = swap_endian(mid)
          for i=0,3 do begin
              cdf_varget, id, vars(pos1+i), tmp, rec_count=nrows
              if (swap) then tmp = swap_endian(tmp)
              xyz(*,*,i) = tmp + mid
          endfor
          if (IsGSE) then begin
              print, "Data in GSE Coordinates - Converting to GSM"
              for i=0,nrows-1 do begin
                  for j=0,3 do begin
                      xgse = xyz(0,i,j)
                      ygse = xyz(1,i,j)
                      zgse = xyz(2,i,j)
                      GSMGSE,XGSM,YGSM,ZGSM,XGSE,YGSE,ZGSE,-1,epoch=times(i)
                      xyz(0,i,j) = xgsm
                      xyz(1,i,j) = ygsm
                      xyz(2,i,j) = zgsm
                  endfor
              endfor
          endif
      endif else begin
          cdf_varget, id, vars(pos), xyz, rec_count=nrows
          if (swap) then xyz = swap_endian(xyz)
          if (IsGSE) then begin
              print, "Data in GSE Coordinates - Converting to GSM"
              for i=0,nrows-1 do begin
                  xgse = xyz(0,i)
                  ygse = xyz(1,i)
                  zgse = xyz(2,i)
                  GSMGSE,XGSM,YGSM,ZGSM,XGSE,YGSE,ZGSE,-1,epoch=times(i)
                  xyz(0,i) = xgsm
                  xyz(1,i) = ygsm
                  xyz(2,i) = zgsm
              endfor
          endif
      endelse

  endif else begin
    print, 'Could not find position array. Please check CDF file!'
    print, 'Here are the variables within this file:'
    for i=0,nvar-1 do print, vars(i)
  endelse

  cdf_close, id 

end
