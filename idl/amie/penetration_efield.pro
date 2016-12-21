

year = fix(ask('year (-1 to enter file names manually)','1998'))

if (year gt 0) then begin

    if (year lt 65) then year = year + 2000 
    if (year gt 65 and year lt 100) then year = year + 1900 

    syear = tostr(year)
    smonth = chopr('0'+tostr(fix(ask('month','01'))),2)

    sday = ask('day (?? acceptable)','??')

    if (strpos(sday,'?') eq -1) then sday = chopr('0'+tostr(fix(sday)),2)

    amie_file_north = syear+'/'+smonth+'/North/b'+syear+smonth+sday+'n'
    amie_file_south = syear+'/'+smonth+'/South/b'+syear+smonth+sday+'s'

endif else begin

    amie_file_north = ask('north file','')
    amie_file_south = ask('south file','')

endelse

NorthFileList = findfile(amie_file_north, count = nFilesNorth)
SouthFileList = findfile(amie_file_south, count = nFilesSouth)

nFilesMin = 0

NorthFile = ''
SouthFile = ''
for i=0,nFilesNorth-1 do begin
    
    for j = 0,nFilesSouth-1 do begin

        if (strpos(SouthFileList(j),strmid(NorthFileList(i),15,8)) gt -1) then begin
            if (nFilesMin eq 0) then begin
                SouthFile = SouthFileList(j)
                NorthFile = NorthFileList(i)
            endif else begin
                SouthFile = [SouthFile,SouthFileList(j)]
                NorthFile = [NorthFile,NorthFileList(i)]
            endelse
            nFilesMin = nFilesMin + 1

        endif

    endfor

endfor

if (strlen(NorthFile(0)) gt 0) then begin
    NorthFileList = NorthFile
    SouthFileList = SouthFile
endif else begin
    print, "Can't match north and south files!"
    print, "NorthFiles : ", NorthFileList
    print, "SouthFiles : ", SouthFileList
    stop
endelse

nFilesTotal = n_elements(NorthFileList)

for iFile = 0, nFilesTotal-1 do begin

    amie_file_north = NorthFileList(iFile)
    amie_file_south = SouthFileList(iFile)

    print, "Reading north file : ", amie_file_north
    read_amie_binary, amie_file_north, data, lats, mlts, time, fields, 	$
      imf, ae, dst, hp, cpcp, version

    NorthPotentialHL = reform(data(*,0,*,*))
    NorthPedersenHL = reform(data(*,1,*,*))
    NorthHallHL = reform(data(*,3,*,*))
    NorthAuroralMeanEnergy = reform(data(*,5,*,*))
    NorthAuroralEnergyFlux = reform(data(*,6,*,*))
    
    print, "Reading south file : ", amie_file_south
    read_amie_binary, amie_file_south, data, lats, mlts, time, fields, 	$
      imf, ae, dst, hp, cpcp, version
    
    SouthPotentialHL = reform(data(*,0,*,*))
    SouthPedersenHL = reform(data(*,1,*,*))
    SouthHallHL = reform(data(*,3,*,*))
    SouthAuroralMeanEnergy = reform(data(*,5,*,*))
    SouthAuroralEnergyFlux = reform(data(*,6,*,*))

    dl = lats(0)-lats(1)

    lowlat = min(lats)

    nLatsHL = n_elements(lats)
    nLatsLL = 2*lowlat/dl + 1
    nMltsLL = n_elements(mlts)
    nTimes  = n_elements(time)

    LatsLL = lowlat - findgen(nLatsLL)*dl
    MltsLL = Mlts

    nLats = n_elements(lats)
    nMlts = n_elements(mlts)

    for iTime = 0, nTimes-1 do begin

        tmp = reform(SouthPedersenHL(iTime,*,*))
        for iLat = 0, nLats-1 do $
          SouthPedersenHL(iTime,*,iLat) = tmp(*,nLats-1-iLat)

        tmp = reform(SouthHallHL(iTime,*,*))
        for iLat = 0, nLats-1 do $
          SouthHallHL(iTime,*,iLat) = tmp(*,nLats-1-iLat)

    endfor

    PotentialLL = fltarr(nTimes,nMltsLL, nLatsLL)
    PedersenLL  = fltarr(nTimes,nMltsLL, nLatsLL)
    HallLL      = fltarr(nTimes,nMltsLL, nLatsLL)

    Theta = fltarr(nMltsLL, nLatsLL)
    for i = 0,nMltsLL-1 do Theta(i,*) = (90.0-LatsLL)*!dtor
    sn = sin(Theta)
    cs = cos(Theta)
    sn2= sn*sn
    cs2 = cs*cs
    cs3 = 1.00 + 3.00*cs2
    cs4 = sqrt(cs3)
    dPsi = (mlts(1) - mlts(0)) / 12.0 * !pi
    dTheta = dL * !dtor

    for iT = 0, nTimes-1 do begin

        c_r_to_a, itime, time(it)
        SeasonAngle = 23.0 * sin((float(jday(itime(0),itime(1),itime(2))) - $
                                  jday(itime(0),3,21))/366*!pi*2)
        utime = itime(3)*3600.0 + itime(4)*60.0 + itime(5)
        CenterTime = 12*3600.0 + 45*60.0
        DipoleAngle = 11.0 * sin((utime-CenterTime)/24.0*3600.0*!pi*2)
        
        Dec = SeasonAngle + DipoleAngle
        SinDec = sin(Dec*!dtor)
        CosDec = cos(Dec*!dtor)

        sza = fltarr(nMlts, nLatsLL)
        cossza = fltarr(nMlts, nLatsLL)

        for iMlt = 0, nMlts-1 do begin 
            cossza(iMlt, *) =  $
              (SinDec*sin(LatsLL*!dtor) + $
               CosDec*cos(LatsLL*!dtor)* $
               cos(!pi*(Mlts(iMlt)-12.0)/12.0))
        endfor

        for i= 0, nMlts-1 do begin

            l = where(reform(cossza(i,*)) lt -0.15, count)

            factor = findgen(nLatsLL)/(nLatsLL-1)

            if (count eq 0) then begin

                add = 0.0
                if (min(cossza(i,*)) lt 0.35) then add = 0.1
                if (min(cossza(i,*)) lt 0.2) then add = 0.15
                if (min(cossza(i,*)) lt 0.0) then add = 0.25

                diff = NorthPedersenHL(iT,i,nLats-1)/$
                  (cossza(i,0)+add)*(1.0-factor) + $
                  SouthPedersenHL(iT,i,0)/(cossza(i,nLatsLL-1)+add)*factor

                cond = (reform(cossza(i,*))+add)*diff
                PedersenLL(iT,i,*) = cond

                diff = NorthHallHL(iT,i,nLats-1)/$
                  (cossza(i,0)+add)*(1.0-factor) + $
                  SouthHallHL(iT,i,0)/(cossza(i,nLatsLL-1)+add)*factor

                cond = (reform(cossza(i,*))+add)*diff
                HallLL(iT,i,*) = cond

            endif else begin

                diff = NorthPedersenHL(iT,i,nLats-1)*(1.0-factor) + $
                  SouthPedersenHL(iT,i,0)*factor
                PedersenLL(iT,i,*) = diff

                diff = NorthHallHL(iT,i,nLats-1)*(1.0-factor) + $
                  SouthHallHL(iT,i,0)*factor
                HallLL(iT,i,*) = diff

            endelse

        endfor

        Sigma0 = 1000.0
        SigmaP = reform(PedersenLL(iT,*,*))
        SigmaH = reform(HallLL(iT,*,*))

        loc = where(SigmaP lt 1.25,count)
        if count gt 0 then sigmap(loc) = 1.25

        loc = where(SigmaH lt 2.5,count)
        if count gt 0 then sigmah(loc) = 2.5

        for i= 1, nMlts-2 do begin
            sigmap(i,*) = (sigmap(i-1,*) + sigmap(i,*) + sigmap(i+1,*))/3
            sigmah(i,*) = (sigmah(i-1,*) + sigmah(i,*) + sigmah(i+1,*))/3
        endfor

        C = 4.00*Sigma0*cs2 + SigmaP*sn2

        SigmaThTh = Sigma0*SigmaP*cs3/C
        SigmaThPs = 2.00*Sigma0*SigmaH*cs*cs4/C
        SigmaPsPs = SigmaP+SigmaH*SigmaH*sn2/C

        dSigmaThTh_dTheta = fltarr(nMltsLL, nLatsLL)
        dSigmaThTh_dPsi   = fltarr(nMltsLL, nLatsLL)
        dSigmaThPs_dTheta = fltarr(nMltsLL, nLatsLL)
        dSigmaThPs_dPsi   = fltarr(nMltsLL, nLatsLL)
        dSigmaPsPs_dTheta = fltarr(nMltsLL, nLatsLL)
        dSigmaPsPs_dPsi   = fltarr(nMltsLL, nLatsLL)

        for j = 0, nMltsLL-1 do begin
            if (j gt 0 and j lt nMlts-1 ) then begin
                for i = 1, nLats-2 do begin
                    dSigmaThTh_dTheta(i,j) = $
                      (SigmaThTh(i+1,j)-SigmaThTh(i-1,j)) / (2*dTheta)
                    dSigmaThTh_dPsi(i,j) = $
                      (SigmaThTh(i,j+1)-SigmaThTh(i,j-1))/  (2*dPsi)

                    dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ $
                      (2*dTheta)
                    dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,j-1))/ $
                      (2*dPsi)

                    dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ $
                      (2*dTheta)
                    dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,j-1))/ $
                      (2*dPsi)
                endfor
            endif else if (j eq 0) then begin
                for i = 1, nLats-2 do begin
                    dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ $
                      (2*dTheta)
                    dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,nMlts-2))/$
                      (2*dPsi)
                    
                    dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/$
                      (2*dTheta)
                    dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,nMlts-2))/$
                      (2*dPsi)

                    dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/$
                      (2*dTheta)
                    dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,nMlts-2))/$
                      (2*dPsi)
                endfor
            endif else begin
                for i = 1, nLats-2 do begin
                    dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/$
                      (2*dTheta)
                    dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,1)-SigmaThTh(i,j-1))/$
                      (2*dPsi)

                    dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/$
                      (2*dTheta)
                    dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,1)-SigmaThPs(i,j-1))/$
                      (2*dPsi)
                    
                    dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/$
                      (2*dTheta)
                    dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,1)-SigmaPsPs(i,j-1))/$
                      (2*dPsi)
                endfor
            endelse
        endfor

        sn  = sin(Theta)
        cs  = cos(Theta) 
        sn2 = sn*sn 
        cs2 = cs*cs 
  
        kappa_Theta2 = SigmaThTh*sn2 
        kappa_Theta1 = SigmaThTh*sn*cs   $
          + dSigmaThTh_dTheta*sn2 - dSigmaThPs_dPsi*sn 
        kappa_Psi2   = SigmaPsPs
        kappa_Psi1   = dSigmaThPs_dTheta*sn + dSigmaPsPs_dPsi
           
        C_A = -2.0 * (kappa_Theta2/dTheta^2 + kappa_Psi2/dPsi^2)
        
        C_B = kappa_Theta2/dTheta^2 - kappa_Theta1/dTheta

        C_C = kappa_Theta2/dTheta^2 + kappa_Theta1/dTheta

        C_D = kappa_Psi2/dPsi^2 - kappa_Psi1/dPsi

        C_E = kappa_Psi2/dPsi^2 + kappa_Psi1/dPsi

; We basically are solving the equation:
;  
;  C_A*phi(i,j) + C_B*phi(i-1,j) + C_C*phi(i+1,j) +
;                 C_D*phi(i,j-1) + C_E*phi(i,j-1) = 0

        ; high lat BC
        PotentialLL(iT, *, 0) = NorthPotentialHL(iT, *, nLatsHL-1)

  ; low lat BC
        PotentialLL(iT, *, nLatsLL-1) = SouthPotentialHL(iT, *, nLatsHL-1)

        if (iT gt 0) then $
          PotentialLL(iT,*,1:nLatsLL-2) = $
          PotentialLL(iT-1,*,1:nLatsLL-2) / $
          max(abs([PotentialLL(iT-1,*,0),PotentialLL(iT-1,*,nLatsLL-1)])) * $
          max(abs([PotentialLL(iT,*,0),PotentialLL(iT,*,nLatsLL-1)]))

        Error = 1.0e6
        PotentialOld = 0.0
        iIters = 0

        Done = 0

        SaveMax = max([abs(PotentialLL(iT,*,0)),abs(PotentialLL(iT,*,nLatsLL-1))])

        while (not done) do begin

            for iLat = 1, nLatsLL-2 do begin

                iMlt = 0
                PotentialLL(iT, iMlt, iLat) = -( $
                       c_c(iMlt,iLat) * PotentialLL(iT, iMlt, iLat+1) + $
                       c_b(iMlt,iLat) * PotentialLL(iT, iMlt, iLat-1) + $
                       c_e(iMlt,iLat) * PotentialLL(iT, iMlt+1, iLat) + $
                       c_d(iMlt,iLat) * PotentialLL(iT, nMltsLL-2, iLat)) / ( $
                       c_a(iMlt,iLat))

                PotentialLL(iT, nMltsLL-1, iLat) = $
                  PotentialLL(iT, iMlt, iLat)

                for iMlt = 1, nMltsLL-2 do begin
                    PotentialLL(iT, iMlt, iLat) = -( $
                      c_c(iMlt,iLat) * PotentialLL(iT, iMlt, iLat+1) + $
                      c_b(iMlt,iLat) * PotentialLL(iT, iMlt, iLat-1) + $
                      c_e(iMlt,iLat) * PotentialLL(iT, iMlt+1, iLat) + $
                      c_d(iMlt,iLat) * PotentialLL(iT, iMlt-1, iLat)) / ( $
                      c_a(iMlt,iLat))

                endfor

            endfor

            PotTemp = reform(PotentialLL(iT,*,*))

            Error = mean((PotentialOld - PotTemp)^2)

            if (error lt 1.0e-6 or max(PotTemp) gt SaveMax) then done = 1

            PotentialOld = PotTemp

            iIters = iIters + 1

        endwhile

;  contour, sigmap, /fill, nlevels = 61, levels = findgen(61)
;  contour, potentialll(iT,*,*), /fill, nlevels = 30, levels = findgen(61)-30
;  plot, potentialll(iT,12,*)
        print, iT, iIters, Error, mm(potentialll(it,*,*))

;stop

    endfor

    nLatsTotal = 90/dl + 1
    nMltsTotal = nMlts

    colats = findgen(nLatsTotal)*dl

    Potential = fltarr(nMltsTotal, nLatsTotal)
    AuroralMeanEnergy = fltarr(nMltsTotal, nLatsTotal)
    AuroralEnergyFlux = fltarr(nMltsTotal, nLatsTotal)

;------------------------------------------------------------------
; Write Equatorial File - North
;------------------------------------------------------------------

    print,'writing file : ',amie_file_north+'_eq'
    openw,1,amie_file_north+'_eq', /f77

    writeu,1, long(nLatsTotal), long(nMltsTotal), long(nTimes)

    writeu,1, colats
    writeu,1, mlts

    writeu,1, 3L

    tmp = bytarr(30)
    writeu,1,'Electric Potential            '
    writeu,1,'Auroral Mean Energy (keV)     '
    writeu,1,'Auroral Energy Flux (W/m2)    '

    for i=0,ntimes-1 do begin

        c_r_to_a, itime, time(i)

        writeu, 1, long(i), long(itime(0:4))

        writeu, 1, imf(i,*), ae(i,*), dst(i,*), hp(i,*), cpcp(i)
        
        for iLat = 0, nLatsHL-1 do $
          Potential(*,iLat) = reform(NorthPotentialHL(i,*,iLat))

        for iLat = nLatsHL, nLatsHL+NLatsLL/2-1 do $
          Potential(*,iLat) = reform(PotentialLL(i,*,iLat-nLatsHL))

        for iLat = 0, nLatsHL-1 do $
          AuroralMeanEnergy(*,iLat) = reform(NorthAuroralMeanEnergy(i,*,iLat))

        for iLat = 0, nLatsHL-1 do $
          AuroralEnergyFlux(*,iLat) = reform(NorthAuroralEnergyFlux(i,*,iLat))

        writeu, 1, Potential*1000.0
        writeu, 1, AuroralMeanEnergy
        writeu, 1, AuroralEnergyFlux

    endfor

    writeu,1,version

    close,1

;------------------------------------------------------------------
; Write Equatorial File - South
;------------------------------------------------------------------

    print,'writing file : ',amie_file_south+'_eq'
    openw,1,amie_file_south+'_eq', /f77

    writeu,1, long(nLatsTotal), long(nMltsTotal), long(nTimes)

    writeu,1, colats
    writeu,1, mlts

    writeu,1, 3L

    tmp = bytarr(30)
    writeu,1,'Electric Potential            '
    writeu,1,'Auroral Mean Energy (keV)     '
    writeu,1,'Auroral Energy Flux (W/m2)    '

    for i=0,ntimes-1 do begin

        c_r_to_a, itime, time(i)

        writeu, 1, long(i), long(itime(0:4))

        writeu, 1, imf(i,*), ae(i,*), dst(i,*), hp(i,*), cpcp(i)

        for iLat = 0, nLatsHL-1 do $
          Potential(*,iLat) = reform(SouthPotentialHL(i,*,iLat))

        for iLat = nLatsHL, nLatsHL+NLatsLL/2-1 do $
          Potential(*,iLat) = reform(PotentialLL(i,*,(nLatsLL-1) - (iLat-nLatsHL)))

        for iLat = 0, nLatsHL-1 do $
          AuroralMeanEnergy(*,iLat) = reform(SouthAuroralMeanEnergy(i,*,iLat))

        for iLat = 0, nLatsHL-1 do $
          AuroralEnergyFlux(*,iLat) = reform(SouthAuroralEnergyFlux(i,*,iLat))

        writeu, 1, Potential*1000.0
        writeu, 1, AuroralMeanEnergy
        writeu, 1, AuroralEnergyFlux

    endfor

    writeu,1,version

    close,1

endfor

end
