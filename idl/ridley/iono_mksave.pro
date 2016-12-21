
list = findfile("-t i*.idl")
if strlen(list(0)) gt 0 then filein = list(0) $
else filein = 'in000000.idl'

filein = ask('filename to make a save file from',filein)

list = findfile(filein)
nfiles = n_elements(list)
if nfiles eq 1 and strlen(list(0)) eq 0 then begin
    print, "I can't seem to find that file."
    stop
endif else nskip = 1

line = ''

for n = 0, nfiles-1,nskip do begin

    swaptheta = 0
    filein = list(n)

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
            endfor

        endif

        if (strpos(mklower(line),"time") gt -1 and $
            strpos(mklower(line),"simulation") lt 0) then begin

            int_tmp = 0
            itime = intarr(6)
            for i=0,5 do begin
                readf, 1, int_tmp
                itime(i) = int_tmp
            endfor

            c_a_to_r, itime, rtime

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

    nt = -1
    for i=0,nvars-1 do if strpos(mklower(vars(i)),'theta') gt -1 then nt = i
    np = -1
    for i=0,5 do if strpos(mklower(vars(i)),'psi') gt -1 then np = i
    if (nt eq -1) or (np eq -1) then begin
        print, "Can't file Theta or Psi variable. Please check file."
        stop
    endif

    theta = reform(data(*,nt,*,*))
    phi   = reform(data(*,np,*,*))

    if (swaptheta) then theta = 90.0 - theta

    save, file=filein+".save",data,nlons,nlats,rtime,vars,theta,phi,nvars

    print, 'Finished Reading File '+filein

endfor

end   
 
