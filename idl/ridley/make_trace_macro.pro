
openw, 1, "field_line.mcr"

printf,1,"#!MC 1000"
printf,1,"$!VarSet |MFBD| = './'"

;printf,1,"$!READDATASET  '""|MFBD|/tecplot_in.plt"" ' "
;printf,1,"  READDATAOPTION = NEW"
;printf,1,"  RESETSTYLE = YES"
;printf,1,"  INCLUDETEXT = NO"
;printf,1,"  INCLUDEGEOM = NO"
;printf,1,"  INCLUDECUSTOMLABELS = NO"
;printf,1,"  VARLOADMODE = BYNAME"
;printf,1,"  VARNAMELIST = '""X [R]"" ""Y [R]"" ""Z [R]"" ""`r [amu/cm^3]"" ""U_x [km/s]"" ""U_y [km/s]"" ""U_z [km/s]"" ""B_x [nT]"" ""B_y [nT]"" ""B_z [nT]"" ""p [nPa]"" ""J_x [`mA/m^2]"" ""J_y [`mA/m^2]"" ""J_z [`mA/m^2]""' " 
printf,1,"$!GLOBALTHREEDVECTOR UVAR = 8"
printf,1,"$!GLOBALTHREEDVECTOR VVAR = 9"
printf,1,"$!GLOBALTHREEDVECTOR WVAR = 10"
printf,1,"$!GLOBALSTREAM RODRIBBON{NUMRODPOINTS = 6}"
printf,1,"$!GLOBALSTREAM RODRIBBON{WIDTH = 0.075}"
printf,1,"$!GLOBALSTREAM CELLFRACTION = 0.5"
printf,1,"$!GLOBALSTREAM MINCELLFRACTION = 1E-02"

r = 10.0
nTheta = 2

;Theta = (findgen(nTheta)+1.0) * !pi/2 / (nTheta+1)
Theta = fltarr(nTheta)
Theta(0) = -20.0 * !pi/180.0
Theta(1) = 20.0 * !pi/180.0

nPhi = 36
Phi = findgen(nPhi) * 2.0*!pi / nPhi

for iTheta = 0, nTheta-1 do begin
    for iPhi = 0, nPhi-1 do begin

        x = r * cos(Theta(iTheta)) * cos(Phi(iPhi))
        y = r * cos(Theta(iTheta)) * sin(Phi(iPhi))
        z = r * sin(Theta(iTheta))

        printf,1,"$!STREAMTRACE ADD"
        printf,1,"  STREAMTYPE = VOLUMEROD"
        printf,1,"  DIRECTION = BOTH"
        printf,1,"  STARTPOS"
        printf,1,"    {"
        printf,1,"    X = ", x
        printf,1,"    y = ", y
        printf,1,"    Z = ", z
        printf,1,"    }"

    endfor
endfor

printf,1,"$!CREATESTREAMZONES" 
printf,1,"  CONCATENATE = NO"

iTotal = 1

for iTheta = 0, nTheta-1 do begin
    for iPhi = 0, nPhi-1 do begin

        sTotal = chopr('00'+tostr(iTotal),3)

        printf,1,"$!WRITEDATASET  ""|MFBD|/trace_"+sTotal+".dat""" 
        printf,1,"  INCLUDETEXT = NO"
        printf,1,"  INCLUDEGEOM = NO"
        printf,1,"  INCLUDECUSTOMLABELS = NO"
        printf,1,"  ASSOCIATELAYOUTWITHDATAFILE = NO"
        printf,1,"  ZONELIST =  ["+tostr(iTotal+1)+"]"
        printf,1,"  VARPOSITIONLIST =  [1-3]"
        printf,1,"  BINARY = NO"
        printf,1,"  USEPOINTFORMAT = YES"
        printf,1,"  PRECISION = 9"

        iTotal = iTotal + 1

    endfor
endfor

printf,1,"$!RemoveVar |MFBD|"

close,1

end
