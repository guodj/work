
pro fix_mags, data, time, LastDayOffset

  LastDayOffset = fltarr(3)*0.0

  nDays = n_elements(time(0,*))
  nPts  = n_elements(time(*,0))

  ; We really can't do too much with too little data...
  if (nDays lt 5) then return

  ; We may have more than 5 days, but all of the data may be missing,
  ; so check to see that we have at least 5 days of good data.
  loc = where(abs(data) lt 10000, count)
  if (count/nPts lt 5) then return

  ; The first thing to check is whether there is a clear linear trend in
  ; the data.

  DailyAverage = fltarr(3, nDays)-99999.0
  DailyStdDev  = fltarr(3, nDays)-99999.0
  DailyTime    = dblarr(nDays)

  for iComp = 0, 2 do begin

      IsBad = 0

      for iDay = 0, nDays-1 do begin
          loc = where(abs(reform(data(iComp,*,iDay))) lt 10000,count)
          if count gt 2 then begin
              DailyAverage(iComp,iDay) = mean(data(iComp,loc,iDay))
              DailyStdDev(iComp,iDay)  = stddev(data(iComp,loc,iDay))
              DailyTime(iDay)          = mean(time(      loc,iDay))
          endif

          ; First check to see if any of the daily values are quite large

          if (abs(DailyAverage(iComp,iDay)) gt $
              abs(2*DailyStdDev(iComp,iDay))) then IsBad = IsBad + 1

      endfor

     ; If we find 4 days that are "bad", then let us consider removing
     ; a line or seeing if there base-line shift

      if (IsBad gt 4) then begin

          print, 'Investigating Component ', iComp,' for base-line shifts'
          print, '   and trend removal.'

          ;----------------------------------------------
          ; First compute a trend line

          DailyTime = (DailyTime - min(DailyTime))/(24.0*3600.0)
          Weights = abs(1.0/reform(DailyStdDev(iComp,*)))

          loc = where(abs(reform(DailyAverage(iComp,*))) lt 10000.0, count)
          yInt      = DailyAverage(iComp,loc(0))
          slope     = (DailyAverage(iComp,loc(count-1)) - yInt) / $
            DailyTime(loc(count-1))

          TestDailyAverage = DailyTime * slope + yInt

          Error = total(Weights*(TestDailyAverage-DailyAverage(iComp,*))^2)

          yInt_save = yInt
          slope_save = slope
          error_save = Error

          for i=1,100 do begin

              yInt  = yInt_save  + 0.1 * (randomu(s,1)-0.5) * yInt_Save
              slope = slope_save + 0.1 * (randomu(s,1)-0.5) * slope_Save
              TestDailyAverage = DailyTime * slope(0) + yInt(0)
              Error = total(Weights*(TestDailyAverage-DailyAverage(iComp,*))^2)

              if (Error lt Error_Save) then begin
                  yInt_save = yInt(0)
                  Slope_Save = Slope(0)
                  Error_Save = Error
              endif

          endfor

          TestDailyAverage = DailyTime * Slope_Save(0) + yInt_Save(0)

          print, "Y-Intercept : ",yInt_Save
          print, "Slope (nT/day) : ", slope_save

          testdata_trend = reform(data(iComp,*,*))

          for iDay = 0, nDays-1 do begin
              loc = where(abs(reform(testdata_trend(*,iDay))) lt 10000,count)
              if count gt 0 then $
                testdata_trend(loc,iDay) = testdata_trend(loc,iDay) - $
                TestDailyAverage(iDay)
          endfor

          ;----------------------------------------------
          ; Second compute a base-line shifted pattern

          loc = where(abs(reform(DailyAverage(iComp,*))) lt 10000.0, count)
          baseline1 = DailyAverage(iComp,loc(0))
          baseline2 = DailyAverage(iComp,loc(count-1))

          iDay = 1
          Done = 0

          while (not done) do begin

              if (iDay ge nDays-1) then Done = 1 else begin

                  if (abs(DailyAverage(iComp,iDay)) lt 10000.0) then begin

                      if (abs(DailyAverage(iComp,iDay)-baseline1) gt $
                          abs(DailyAverage(iComp,iDay)-baseline2)) then $
                        Done=1 $
                      else iDay = iDay + 1

                  endif else iDay = iDay + 1

              endelse

          endwhile

          iDay = iDay - 1

          iDayShift = iDay

          baseline1 = 0.0
          nDays1 = 0
          baseline2 = 0.0
          nDays2 = 0

          for iDay = 0, nDays-1 do begin
              if (iDay le iDayShift) then begin
                  if (abs(DailyAverage(iComp,iDay)) lt 10000) then begin
                      baseline1 = baseline1 + DailyAverage(iComp,iDay)
                      nDays1 = nDays1 + 1
                  endif
              endif else begin
                  if (abs(DailyAverage(iComp,iDay)) lt 10000) then begin
                      baseline2 = baseline2 + DailyAverage(iComp,iDay)
                      nDays2 = nDays2 + 1
                  endif
              endelse
          endfor

          baseline1 =baseline1/nDays1
          baseline2 =baseline2/nDays2

          testdata_shift = reform(data(iComp,*,*))

          print, "Day Boundary : ", iDayShift
          print, "Baseline 1 :",baseline1
          print, "Baseline 2 :",baseline2

          ; One problem with this method is that we may have storms that
          ; cause these "shifts".  Have to do something about this...

          if (iDayShift lt nDays-1) then begin

              ; lets compare the last few points of the previous day
              ; and the first few points of the next day....

              l1 = where(abs(data(iComp,*,iDayShift)) lt 10000.0, c1)
              l2 = where(abs(data(iComp,*,iDayShift+1)) lt 10000.0, c2)

              if (c1 gt 10 and c2 gt 10) then begin

                  m1 = mean(data(iComp,l1(c1-11:c1-1),iDayShift))
                  m2 = mean(data(iComp,l1(0:10),iDayShift+1))

                  ; if the base line shift happens in the first few points,
                  ; then go ahead, otherwise, just skip it....

                  if (abs(m2-m1) gt 0.25*abs(baseline2-baseline1)) then begin
              
                      for iDay = 0, nDays-1 do begin

                          if (iDay le iDayShift) then bl = baseline1 else bl = baseline2

                          loc = where(abs(reform(testdata_shift(*,iDay))) lt 10000,count)
                          if count gt 0 then $
                            testdata_shift(loc,iDay) = testdata_shift(loc,iDay) - bl

                      endfor

                  endif else begin

                      print, "Skipping baseline shift correction : "
                      print, "m1, m2 : ", m1, m2
                      baseline1 = 0.0
                      baseline2 = 0.0
                      
                  endelse

              endif

          endif

          loc = where(abs(testdata_trend) lt 10000.0, count)
          if (count gt 0) then ms_trend = mean(testdata_trend(loc)^2)
          loc = where(abs(testdata_shift) lt 10000.0, count)
          if (count gt 0) then ms_shift = mean(testdata_shift(loc)^2)

          print, "ms_trend : ", ms_trend
          print, "ms_shift : ", ms_shift

          if (ms_trend lt ms_shift) then begin
              data(iComp,*,*) = testdata_trend 
              loc = where(abs(reform(DailyAverage(iComp,*))) lt 10000.0, count)
              LastDayOffset(iComp) = $
                DailyTime(loc(count-1)) * Slope_Save(0) + yInt_Save(0)
          endif else begin
              data(iComp,*,*) = testdata_shift
              LastDayOffset(iComp) = baseline2
          endelse

      endif

  endfor

end
