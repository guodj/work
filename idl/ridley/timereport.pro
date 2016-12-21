
filelist = findfile("/Users/ridley/.timekeeper_output.*")

nFiles = n_elements(filelist)

time   = dblarr(nFiles,10000)
items  = strarr(nFiles,10000)
nItems = lonarr(nFiles)
line   = ""

for iFile=0,nFiles-1 do begin

  filename = filelist(iFile)
  openr,1,filename

  nItems(iFile) = 0

  while not eof(1) do begin

      readf,1,line

      l = strpos(line,'<')
      items(iFile,nItems(iFile)) = strmid(line,0,l-1)
      stime = strmid(line,l+1,18)
      c_s_to_a, itime, stime
      c_a_to_r, itime, rtime
      time(iFile,nItems(iFile)) = rtime
      nItems(iFile) = nItems(iFile) + 1

  endwhile

  close,1

endfor

Tasks_Cum    = strarr(1000)
Times_Cum    = fltarr(1000)
Tasks_Cum(0) = items(0,0)
nTasks_Cum   = 1

Tasks_Cum_Cat    = strarr(1000)
Times_Cum_Cat    = fltarr(1000)
Tasks_Cum_Cat(0) = strmid(items(0,0), 0, strpos(items(0,0),':'))
nTasks_Cum_Cat   = 1

for iFile=0,nFiles-1 do begin

    Tasks_Daily = strarr(1000)
    Times_Daily = fltarr(1000)
    Tasks_Daily(0) = items(iFile,0)
    nTasks_Daily = 1

    for iItem = 1, nItems(iFile) - 1 do begin

        dt = time(iFile,iItem) - time(iFile,iItem-1)

        ; Daily

        iFound = -1
        for i = 0, nTasks_Daily-1 do begin
;            print, "-->",items(iFile, iItem),"<-->",Tasks_Daily(i),"<--"
            if (strpos(items(iFile,iItem),Tasks_Daily(i)) eq 0) then iFound = i
        endfor
;        print, "iFound : ",iFound, nTasks_Daily
        if (iFound gt -1) then begin
            Times_Daily(iFound) = Times_Daily(iFound) + dt
        endif else begin
            Tasks_Daily(nTasks_Daily) = items(iFile,iItem)
            Times_Daily(nTasks_Daily) = dt
            nTasks_Daily = nTasks_Daily + 1
        endelse

        ; Cummulative

        iFound = -1
        for i = 0, nTasks_Cum-1 do $
          if (strpos(items(iFile, iItem),Tasks_Cum(i)) eq 0) then iFound = i
        if (iFound gt -1) then begin
            Times_Cum(iFound) = Times_Cum(iFound) + dt
        endif else begin
            Tasks_Cum(nTasks_Cum) = items(iFile,iItem)
            Times_Cum(nTasks_Cum) = dt
            nTasks_Cum = nTasks_Cum + 1
        endelse

        ; Cummulative

        iFound = -1

        cat = strmid(items(iFile, iItem),0, strpos(items(iFile, iItem),':'))
        for i = 0, nTasks_Cum_Cat-1 do begin
            if (strpos(Tasks_Cum_Cat(i), cat) eq 0) then iFound = i
        endfor

        if (iFound gt -1) then begin
            Times_Cum_Cat(iFound) = Times_Cum_Cat(iFound) + dt
        endif else begin
            Tasks_Cum_Cat(nTasks_Cum_Cat) = cat
            Times_Cum_Cat(nTasks_Cum_Cat) = dt
            nTasks_Cum_Cat = nTasks_Cum_Cat + 1
        endelse

    endfor

    c_r_to_a, itime, time(iFile,0)
    c_a_to_s, itime, stime

    print, ""
    print, "Daily Report for : ",strmid(stime,0,9)

    for i=0,nTasks_Daily-1 do begin
        td = Times_Daily(i)
        if (td gt 120.0) then begin
            td = td/60.0
            if (td gt 120.0) then begin
                td = td/60.0
                units = " hours"
            endif else units = " minutes"
        endif else units = " seconds"
        print, format="(a40,f10.2, a)", Tasks_Daily(i), td, units
    endfor

endfor

c_r_to_a, itime, time(0,0)
c_a_to_s, itime, stime1

c_r_to_a, itime, time(nFiles-1,0)
c_a_to_s, itime, stime2

print, ""
print, "Categorical Cummulative Report for : ",strmid(stime1,0,9), $
  " thru ",strmid(stime2,0,9)

for i=0,nTasks_Cum_Cat-1 do begin
    td = Times_Cum_Cat(i)
    if (td gt 120.0) then begin
        td = td/60.0
        if (td gt 120.0) then begin
            td = td/60.0
            units = " hours"
        endif else units = " minutes"
    endif else units = " seconds"
    print, format="(a40,f10.2,a)", Tasks_Cum_Cat(i), td, units
endfor

print, ""
print, "Cummulative Report for : ",strmid(stime1,0,9), $
  " thru ",strmid(stime2,0,9)

for i=0,nTasks_Cum-1 do begin
    td = Times_Cum(i)
    if (td gt 120.0) then begin
        td = td/60.0
        if (td gt 120.0) then begin
            td = td/60.0
            units = " hours"
        endif else units = " minutes"
    endif else units = " seconds"
    print, format="(a40,f10.2,a)", Tasks_Cum(i), td, units
endfor

end
