
pro read_mag_obs, file, nObs, stats, lats, uts

    spawn, 'mag_obs.txt'

    line = ''

    openr,1,'mag_obs.txt'

    readf,1,line

    if (eof(1)) then begin
        nObs = 0
        return
    endif

    if (strpos(line,'CGM') gt -1) then begin

        for i=1,5 do readf,1,line

        nObs = 14

        stats = strarr(nObs)
        lats = fltarr(nObs)
        uts = fltarr(nObs)

        for i=0,5 do begin
            readf,1,line
            stats(i) = strmid(line, 0, 2)
            lats(i) = -float(strmid(line,45,6))
            uts(i) = float(strmid(line,70,2))*60.0 + float(strmid(line,73,2))
        endfor

        readf,1,line

        for i=6,nObs-1 do begin
            readf,1,line
            stats(i) = strmid(line, 0, 20)
            if (strmid(line,44,1) eq 'N') then s=1 else s=-1
            lats(i) = s*float(strmid(line,45,6))
            uts(i) = float(strmid(line,70,2))*60.0 + float(strmid(line,73,2))
        endfor

        close,1

    endif else begin

        print, "don't understand mag file."
        nObs = 0

    endelse


end

