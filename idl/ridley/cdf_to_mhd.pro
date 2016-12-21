
filein  = ask('file to convert','')
satfile = ask('whether this is a satellite file or an IMF file (s/i)','s')
fileout = ask('file to output','')

if (mklower(satfile) eq 's') then begin

    print, "File : ",filein
    filein = findfile(filein)
    filein = filein(0)
    if (strlen(filein) eq 0) then begin
        print, "File not found!!"
        stop
    endif

    read_cdf_pos, filein, xyz, time

    if (min(xyz) gt -1.0e31) then begin

        r = sqrt(xyz(0,*)^2 + xyz(1,*)^2 + xyz(2,*)^2)
        loc = where(r gt 6372.0, count)
        if (count gt 0) then xyz = xyz/6372.0
        
        ntimes = n_elements(time)
        
        s = size(xyz)

        if (s(0) eq 2) then begin

            openw,1,fileout
            
            printf,1,''
            printf,1,'Data taken from file : ', filein
            printf,1,''
            printf,1,'#START'

            for i=0,ntimes-1 do begin

                cdf_epoch, time(i), year, month, day, hour, minute, second, milli, /break

                printf,1, format = '(i5,5i3,i4,3f8.2)', $
                  year, month, day, hour, minute, second, milli, $
                  xyz(0,i), xyz(1,i), xyz(2,i)

            endfor

            close,1

        endif else begin

            for n=0,s(3)-1 do begin

                p = strpos(fileout,"cluster")
                if (p ge 0) then begin
                   fo = strmid(fileout,0,p+7)+tostr(n)+strmid(fileout,p+7,strlen(fileout))
                endif else begin
                   p = strpos(fileout,"dat")
                   fo = strmid(fileout,0,p-1)+"_"+tostr(n)+".dat"
                endelse
                openw,1,fo
            
                printf,1,''
                printf,1,'Data taken from file : ', filein
                printf,1,''
                printf,1,'Satellite Number : ',tostr(n)
                printf,1,''
                printf,1,'#START'

                for i=0,ntimes-1 do begin

                    cdf_epoch, time(i), year, month, day, hour, minute, second, milli, /break

                    printf,1, format = '(i5,5i3,i4,3f8.2)', $
                      year, month, day, hour, minute, second, milli, $
                      xyz(0,i,n), xyz(1,i,n), xyz(2,i,n)

                endfor

                close,1

            endfor

        endelse

    endif


endif else begin

;    ac_h0_mfi_20031023_v04.cdf
;    year 

    filein = findfile(filein)
    filein = filein(0)
    if (strlen(filein) eq 0) then begin
        print, "File not found!!"
        stop
    endif

  magfile = filein
  swfile  = ask('SOLAR WIND file to convert','')

  swfile = findfile(swfile)
  swfile = swfile(0)
  if (strlen(swfile) eq 0) then begin
      print, "File not found!!"
      stop
  endif

  print, "magfile : ", magfile
  print, "swfile  : ", swfile
  
  read_cdf_imf, magfile, swfile, magfield, sw, magtime, swtime, isgse, $
    delay, xyz

  ntimes = n_elements(magtime)
  ntsw   = n_elements(swtime)

  den = fltarr(ntimes)
  alp = fltarr(ntimes)
  tem = fltarr(ntimes)
  vel = fltarr(3,ntimes)

  if (n_elements(xyz(0,*)) ne ntimes) then xyz2 = fltarr(3,ntimes)

  if (max(sw(5,*)) gt 0.0) then begin
      llt0 = where(sw(5,*) lt 0.0,count)
      if (count gt 0) then begin
          lgt0 = where(sw(5,*) gt 0.0)
          sw(5,llt0) = mean(sw(5,lgt0))
      endif
      AlphaFac = 4.0
  endif else AlphaFac = 0.0

  ; The solar wind data is often screwed up.  Lets try to
  ; fix it a bit before proceeding.

  j = 0
  if (sw(0,j) lt -2500.0 or $
      sw(0,j) gt -100.0) then begin
     sw(0,j) = -400.0
     sw(1,j) = 0.0
     sw(2,j) = 0.0
  endif

  for j = 1, ntsw-1 do begin
     if (sw(0,j) lt -5000.0 or $
         sw(0,j) gt -100.0) then sw(*,j) = sw(*,j-1)
  endfor

  for i=0,ntimes-1 do begin

      if (abs(magfield(0,i)) lt 1000.0) then begin

          cdf_epoch, magtime(i), year, month, day, hour, $
            minute, second, milli,/break

          j = 0
          done = 0
          while (not done) do begin
              if (j eq ntsw-1) then done = 1
              if (swtime(j) gt magtime(i) and $
                  sw(0,j) gt -5000.0 and $
                  sw(0,j) lt -100.0) then done = 1
              if (not done) then j++
          endwhile

          done = 0
          j2 = j
          while (not done) do begin
              if (j2 eq ntsw-1) then done = 1
              if (swtime(j2) gt magtime(i) and sw(3,j2) gt -5000.0) then done = 1
              if (not done) then j2 = j2 + 1
          endwhile

          if (abs(sw(3,j2)) lt 150.0) then begin
              den(  i) = sw(  3,j2) * (1.0 + sw(5,j2)*AlphaFac)
          endif else den(i) = 5.0 * (1.0 + 0.1*AlphaFac)

          done = 0
          j3 = j
          while (not done) do begin
              if (j3 ge ntsw-2) then begin
                  done = 1
              endif else begin
                  if (swtime(j3) gt magtime(i) and sw(4,j3) gt -5000.0) then done = 1
                  j3 = j3 + 1
              endelse
          endwhile

          if (n_elements(xyz(0,*)) ne ntimes) then xyz2(*,i) = xyz(*,j)

          if (j eq ntsw-1 or j3 eq 1) then begin
              tem(  i) = sw(  4,j3)
              vel(*,i) = sw(0:2,j)
          endif else begin
              if (sw(  4,j3) le sw(  4,j3+1)) then begin
                  tem(  i) = sw(  4,j3)
              endif else begin
                  dt = (swtime(j3) - magtime(i)) / (swtime(j3+1) - swtime(j3))
                  tem(  i) = dt * sw(  4,j3) + (1.0-dt) * sw(4,j3+1)
              endelse
              v1 = sqrt(sw(0,  j)^2+sw(1,  j)^2+sw(2,  j)^2)
              v2 = sqrt(sw(0,j+1)^2+sw(1,j+1)^2+sw(2,j+1)^2)
              if (v1 le v2) then begin
                  vel(*,i) = sw(0:2,j)
              endif else begin
                  dt = (swtime(j) - magtime(i)) / (swtime(j+1) - swtime(j))
                  vel(*,i) = dt * sw(0:2,j) + (1.0-dt) * sw(0:2,j+1)
              endelse
          endelse

      endif

  endfor

  j = 0
  for i=0,ntimes-1 do begin
      if (abs(magfield(0,i)) lt 1000.0) then begin
          if (i ne j) then begin
              magfield(*,j) = magfield(*,i)
              magtime(j)    = magtime(i)
              vel(*,j)      = vel(*,i)
              xyz(*,j)      = xyz(*,i)
              den(j)        = den(i)
              tem(j)        = tem(i)
          endif
          j = j + 1
      endif
  endfor

  l = where(tem lt 0.0, count)
  if (count gt 0) then tem(l) = 150000.0

  if (j ne nTimes) then begin

      nTimes = j

      mf = fltarr(3,nTimes)
      mf(*,0:nTimes-1) = magfield(*,0:nTimes-1)
      magfield = mf

      v = fltarr(3,nTimes)
      v(*,0:nTimes-1) = vel(*,0:nTimes-1)
      vel = v

      x = fltarr(3,nTimes)
      x(*,0:nTimes-1) = xyz(*,0:nTimes-1)
      xyz = x

      magtime    = magtime(0:nTimes-1)
      den        = den(0:nTimes-1)
      tem        = tem(0:nTimes-1)

  endif

  if (n_elements(xyz(0,*)) ne ntimes) then xyz = xyz2

  dynamic = mklower( $
            ask('whether you want dynamic or static time delay (d/s/n)','d'))

  if (strpos(dynamic,'d') gt -1) then begin

      print, "Dynamically adjusting time delay"

      newmagtime = dblarr(nTimes)

      newmagfield = magfield
      newvel = vel
      newden = den
      newtem = tem

      magtime = magtime - (xyz(0,*)-32.0*6372.0)/vel(0,*) * 1000.0

      print, "mean x          : ", mean(xyz(0,*))/6372.0, $
             mean(xyz(1,*))/6372.0,mean(xyz(2,*))/6372.0
      print, "mean vx         : ", mean(vel(0,*))
      print, "Mean time delay : ",-mean((xyz(0,*)-32.0*6372.0)/vel(0,*))/60.0,' minutes'

      goodtimes = intarr(nTimes)

      IsDone = 0

      MaxTime = max(magtime)+1.0

      i = nTimes+1
      j = nTimes

      print, mm(magtime)
      magtime_tmp = reform(magtime-min(magtime))/3600000.0
      print, mm(magtime_tmp)

      while not IsDone do begin

          loc = where(magtime_tmp(0:j-1) eq max(magtime_tmp(0:j-1)))

          if (magtime_tmp(loc(0)) gt MaxTime or $
              magtime_tmp(loc(0)) lt 0) then IsDone = 1 else begin
              if MaxTime lt 0 then IsDone = 1 else begin
                  j = loc(0)
                  goodtimes(loc(0)) = 1
                  MaxTime = magtime_tmp(loc(0))
                  magtime_tmp(loc(0)) = -1.0
              endelse
          endelse
          if (j eq 0) then IsDone = 1
      endwhile

      j = 0
      for i = 0, nTimes-1 do begin
          if (goodtimes(i)) then begin
              newmagtime(j) = magtime(i)
              newmagfield(0,j) = newmagfield(0,i)
              newmagfield(1,j) = newmagfield(1,i)
              newmagfield(2,j) = newmagfield(2,i)
              newvel(0,j) = newvel(0,i)
              newvel(1,j) = newvel(1,i)
              newvel(2,j) = newvel(2,i)
              newden(  j) = newden(i)
              newtem(  j) = newtem(i)
              j = j + 1
          endif
      endfor

      nTimes = j

      magfield = newmagfield
      vel = newvel
      den = newden
      tem = newtem
      magtime = newmagtime

      delay = 0.0

  endif

  if (strpos(dynamic,'s') gt -1) then begin

      print, "Statically adjusting time delay - within the data"

      magtime = magtime - mean((xyz(0,*)-32.0*6372.0)/vel(0,*)) * 1000.0
      print, "Delay : ",delay, -mean((xyz(0,*)-32.0*6372.0)/vel(0,*))

      delay = 0.0

  endif

  openw,1,fileout

  printf,1,''
  printf,1,'Data taken from file : ', filein
  printf,1,'     and               ', swfile
  if (AlphaFac gt 1.0) then begin
      printf,1,''
      printf,1,'The number density was multiplied by (1.0 + 4.0*AlphaRatio)'
  endif else begin
      printf,1,''
      printf,1,'The Alpha Particles were not considered in the number density'
  endelse
  printf,1,''
  printf,1,'#TIMEDELAY'
  printf,1,delay
  printf,1,''
  if (isgse) then begin
    printf,1,'#COOR'
    printf,1,'GSE'
    printf,1,''
  endif
  printf,1,'#START'

  for i=0,ntimes-1 do begin

      if (abs(magfield(0,i)) lt 1000.0) then begin

          cdf_epoch, magtime(i), year, month, day, hour, minute, second, milli,/break

          printf,1, format = '(i5,5i3,i4,3f8.2,4f9.2,f13.1)', $
            year, month, day, hour, minute, second, milli, $
            magfield(0,i), magfield(1,i), magfield(2,i),   $
            vel(*,i), den(i), tem(i)

      endif

  endfor

  close,1

endelse

end
