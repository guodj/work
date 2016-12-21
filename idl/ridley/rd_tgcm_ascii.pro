pro rd_tgcm_ascii, filein, data, vars, utdes, axdes, xax, yax, ntmax, run

  openr,1,filein

  nin = 0

  vt = ''
  ut = ''
  temp = ''
  run_temp = ''

  t = fltarr(6)

  while (not eof(1)) do begin

    readf,1, vt
    readf,1, ut
    print, ut
    readf,1, run_temp
    readf,1, temp

; read in x positions and information:

    readf,1, temp

    nx = fix(strmid(temp,3,4))
    ax = strmid(temp,12, strlen(temp)-12)
    ax = stripend(ax)

    n  = nx/6
    n2 = nx - n*6
    if n2 gt 0 then t2 = fltarr(n2)
    v1 = fltarr(nx)
    for i=0,n-1 do begin
      readf,1,t
      v1(i*6:(i+1)*6-1) = t
    endfor
    if n2 gt 0 then begin
      readf,1,t2
      v1(n*6:nx-1) = t2
    endif

; read in y positions and information:

    readf,1, temp

    ny = fix(strmid(temp,3,4))
    ay = strmid(temp,12, strlen(temp)-12)
    ay = stripend(ay)

    n  = ny/6
    n2 = ny - n*6
    if n2 gt 0 then t2 = fltarr(n2)
    v2 = fltarr(ny)
    for i=0,n-1 do begin
      readf,1,t
      v2(i*6:(i+1)*6-1) = t
    endfor
    if n2 gt 0 then begin
      readf,1,t2
      v2(n*6:ny-1) = t2
    endif

    if nin eq 0 then begin

      data  = fltarr(ntmax,nx,ny)
      vars  = strarr(ntmax)
      utdes = strarr(ntmax)
      run   = strarr(ntmax)
      axdes = strarr(2,ntmax)
      xax   = v1
      yax   = v2

    endif else begin

      if nin ge ntmax then begin
        nin = ntmax - 1
        print, 'Too many data points in the file. Erasing last elements'
      endif

    endelse

    vars(nin)    = vt
    utdes(nin)   = ut
    run(nin)     = run_temp
    axdes(0,nin) = ax
    axdes(1,nin) = ay

    nt = nx*ny
    n = nt/6
    n2 = ny - n*6
    if n2 gt 0 then t2 = fltarr(n2)
    v2 = fltarr(ny)
    for i=0,n-1 do begin
      readf,1,t
      for j=0,5 do begin
        ij = i*6+j
        x  = ij mod nx
        y  = ij/nx 
	data(nin,x,y) = t(j)
      endfor
    endfor
    if n2 gt 0 then begin
      readf,1,t2
      for j=0,n2-1 do begin
        ij = n*6+j
        x  = ij mod nx
        y  = ij/nx 
	data(nin,x,y) = t2(j)
      endfor
    endif

    nin = nin + 1

  endwhile

  close,1

  data  = data(0:nin-1,*,*)
  utdes = utdes(0:nin-1)
  vars  = vars(0:nin-1)
  run   = run(0:nin-1)
  axdes = axdes(*,0:nin-1)

  return

end
