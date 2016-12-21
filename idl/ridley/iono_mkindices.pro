
filein = 'it*.idl'
filein = ask('filename',filein)

list = findfile(filein)
nfiles = n_elements(list)
if nfiles eq 1 and strlen(list(0)) eq 0 then begin
  print, "I can't seem to find that file."
  stop
endif else begin

  if (nfiles gt 1) then begin
    print, "I found ",tostr(nfiles)," to plot."
    nskip = fix(ask('how many to skip between them','0'))+1
  endif else nskip = 1

endelse

line = ''

if (nfiles gt 1) then begin

  relative = mklower(ask('universal time (U), or relative time (R)','U'))

  dt = 0.0
  if (strpos(relative,'r') eq 0) then begin
    dt = float(ask('offset time on first file (in seconds)','60.0'))
    ut = ''
  endif else ut = ' UT'

  if (strpos(relative,'r') ne 0) then ut = ' UT'

endif else begin

  ut = ''
  relative = 'u'
  dt = 0.0

endelse

time = dblarr(nfiles/nskip + 1)
itime = intarr(6)

ntotal = 0

for n = 0, nfiles-1,nskip do begin

  swaptheta = 0
  filein = list(n)

  if strpos(filein,'save') gt -1 then begin

    restore, filein

    if (n eq 0) then begin
        for i=0,nvars-1 do begin
            if n eq 0 then print, chopr(' '+tostr(i),2),'. ',vars(i)
        endfor
    endif

    time(ntotal) = rtime

  endif else begin

      openr,1,filein

      done = 0

      while (not done) do begin

          readf,1, line

          if (strpos(mklower(line),"numerical") gt -1) then begin

              readf,1, nvars
              readf,1, nlats
              readf,1, nlons

              tmp = fltarr(nvars)
              vars = strarr(nvars)

          endif

          if (strpos(mklower(line),"variable") gt -1) then begin

              for i=0,nvars-1 do begin
                  readf,1,line
                  vars(i) = strmid(line,6,strlen(line)-6)
                  if n eq 0 then print, chopr(' '+tostr(i),2),'. ',vars(i)
              endfor

          endif

          if (strpos(mklower(line),"time") gt -1 and $
              strpos(mklower(line),"simulation") lt 0) then begin

              int_tmp = 0
              for i=0,5 do begin
                  readf, 1, int_tmp
                  itime(i) = int_tmp
              endfor

              c_a_to_r, itime, rtime
              time(ntotal) = rtime

          endif

          if (strpos(mklower(line),"northern") gt -1) then begin

              data = fltarr(2,nvars,nlons,nlats)
              for j=0,nlons-1 do for i=0,nlats-1 do begin
                  readf,1,tmp
                  data(0,*,j,i) = tmp
              endfor

          endif

          if (strpos(mklower(line),"all") gt -1) then begin

              nlons = nlons+1
              nlats = nlats/2
              data = fltarr(2,nvars,nlons,nlats)
              for j=0,nlons-2 do begin 
                  for i=nlats-1,0,-1 do begin
                      readf,1,tmp
                      data(1,*,j,i) = tmp
                  endfor
                  for i=nlats-1,0,-1 do begin
                      readf,1,tmp
                      data(0,*,j,i) = tmp
                  endfor
              endfor
              swaptheta = 1

              data(*,*,nlons-1,*) = data(*,*,0,*)

          endif

          if (strpos(mklower(line),"southern") gt -1) then begin

              for j=0,nlons-1 do for i=0,nlats-1 do begin
                  readf,1,tmp
                  data(1,*,j,i) = tmp
              endfor

          endif

          if eof(1) then done = 1
          
      endwhile

      close,1

      if (n eq 0) then begin
          nt = -1
          for i=0,nvars-1 do $
            if strpos(mklower(vars(i)),'theta') gt -1 then nt = i
          np = -1
          for i=0,5 do if strpos(mklower(vars(i)),'psi') gt -1 then np = i
          if (nt eq -1) or (np eq -1) then begin
              print, "Can't file Theta or Psi variable. Please check file."
              stop
          endif

          theta = reform(data(*,nt,*,*))
          phi   = reform(data(*,np,*,*))

          if (swaptheta) then theta = 90.0 - theta
      endif

  endelse

  if (n eq 0) then begin

      var = 0
      nvars_to_plot = 0
      while (var ge 0) do begin
        var = fix(ask('Variable Number to plot (-1 to exit)','-1'))
        if (var ge 0) and (var lt nvars) then begin
          if nvars_to_plot eq 0 then $
            varlist = [var]          $
          else varlist = [varlist,var]
          nvars_to_plot = nvars_to_plot + 1
        endif
      endwhile

      print, "You have selected "+tostr(nvars_to_plot)+" variables to plot." 

      if (nvars_to_plot eq 0) then begin
	print, "Can not continue!"
	stop
      endif

      data_to_plot = fltarr(2,nfiles/nskip+1,nvars_to_plot,nlons,nlats)

  endif

  data_to_plot(*,ntotal,*,*,*) = reform(data(*,varlist,*,*))

  print, 'Finished Reading File '+filein
  ntotal = ntotal + 1

endfor

toff = 0.0
if (strpos(relative,'r') eq 0) then toff = time(0) - dt - 24.0*3600.0 else toff = -dt
time = time(0:ntotal-1) - toff

;
; We want to write an indices file, so figure out maximum & minimum for all
; Selected variables
;

indices = fltarr(2,ntotal,nvars_to_plot,2)
indices_vars = strarr(nvars_to_plot)

for i = 0, nvars_to_plot-1 do indices_vars(i) = vars(varlist(i))

for hem = 0, 1 do for i = 0, nvars_to_plot-1 do for j = 0, ntotal-1 do begin
  indices(hem,j,i,*) = mm(data_to_plot(hem,j,i,*,*))
endfor

save, indices, indices_vars, time, file = list(0)+'.ind.save'

end    
