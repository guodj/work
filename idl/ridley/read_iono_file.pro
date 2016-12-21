
pro read_iono_file, list, data_to_plot, time, theta, phi, $
                    nskip = nskip, variables = variables, all = all

  if n_elements(variables) eq 0 then variables = -1
  if n_elements(nskip) eq 0 then nskip = 1
  if n_elements(all) eq 0 then all = 0

  nfiles = n_elements(list)

  ntimestotal = (nfiles + (nfiles mod nskip))/nskip

  time = dblarr(ntimestotal)
  itime = intarr(6)

  ntotal = 0

  line = ''

  for n = 0, nfiles-1,nskip do begin

    swaptheta = 0
    filein = list(n)

    if strpos(filein,'save') gt -1 then begin
      restore, filein
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
	    if (n eq 0) and (variables(0) eq -1) then $
              print, chopr(' '+tostr(i),2),'. ',vars(i)
          endfor

        endif

        if (strpos(mklower(line),"time") gt -1) then begin

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

        if (all) then begin
          varlist = indgen(nvars)
	  nvars_to_plot = nvars
        endif else begin
          if (variables(0) gt -1) then begin
            varlist = variables
	    nvars_to_plot = n_elements(varlist)
          endif else begin

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

            print, "You have selected "+tostr(nvars_to_plot)+$
                   " variables to plot." 

            if (nvars_to_plot eq 0) then begin
	      print, "Can not continue!"
	      stop
            endif

          endelse

        endelse

        data_to_plot = fltarr(2,ntimestotal,nvars_to_plot,nlons,nlats)
        theta = fltarr(2,nlons,nlats)
        phi   = fltarr(2,nlons,nlats)

        nt = -1
        for i=0,nvars-1 do $
           if strpos(mklower(vars(i)),'theta') gt -1 then nt = i
        np = -1
        for i=0,5 do if strpos(mklower(vars(i)),'psi') gt -1 then np = i
        if (nt eq -1) or (np eq -1) then begin
          print, "Can't file Theta or Psi variable. Please check file."
          stop
        endif

        theta(*,*,*) = reform(data(*,nt,*,*))
        phi(*,*,*)   = reform(data(*,np,*,*))

        if (swaptheta) then theta = 90.0 - theta

      endif

      data_to_plot(*,ntotal,*,*,*) = reform(data(*,varlist,*,*))

    endelse

    print, 'Finished Reading File '+filein
    ntotal = ntotal + 1

  endfor

  return

end
