
pro plot_data, adata, alon, alat, idata, ilon, ilat, $
               pn, ppp, space, sizes, units

  aloc = where(alat gt 50.0)
  iloc = where(ilat gt 50.0)
  mini = min([min(adata(*,aloc)),min(idata(*,iloc))])
  maxi = max([max(adata(*,aloc)),max(idata(*,iloc))])

  get_position, ppp, space, sizes, pn, pos
  contour_circle, adata, alon, alat, $
    mini = mini, maxi = maxi, /sym, $
    nLevels = 31, $
    pos = pos, $
    maxrange = 40.0

  if (pn eq 0) then begin
      xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.03, 'AMIE', alignment = 0.5, /norm
  endif else begin
      xyouts, (pos(0)+pos(2))/2.0, pos(1)-0.03, 'AMIE', alignment = 0.5, /norm
  endelse

  get_position, ppp, space, sizes, pn+1, pos
  contour_circle, idata, ilon, ilat, $
    mini = mini, maxi = maxi, /sym, $
    nLevels = 31, $
    pos = pos, $
    maxrange = 40.0, /no06

  ctpos = [pos(2)+0.01, pos(1), pos(2)+0.03, pos(3)]
  plotct, 255, ctpos, [mini,maxi], title, /right

  if (pn eq 0) then begin
      xyouts, (pos(0)+pos(2))/2.0, pos(3)+0.03, 'SWMF', alignment = 0.5, /norm
  endif else begin
      xyouts, (pos(0)+pos(2))/2.0, pos(1)-0.03, 'SWMF', alignment = 0.5, /norm
  endelse

end

;--------------------------------------------------------------------------
; Main Code
;--------------------------------------------------------------------------

Amie_File = 'b20020110n.save'

read_amie_binary, Amie_File, Amie_Data, Amie_Lats, Amie_Mlts, $
  Amie_Times, Amie_Vars, imf, ae, dst, hp, pot, version

filelist = findfile("it*.idl")
nfiles = n_elements(filelist)

for iFile = 0, nFiles-1, 10 do begin

    Iono_File = filelist(iFile)

    iono_read_file, Iono_File, Iono_nVars, Iono_nLats, Iono_nLons, $
      Iono_Vars, Iono_Time, Iono_Data

    c_r_to_a, itime, Iono_Time
    c_a_to_s, itime, stime
    stime = strmid(stime,0,15)+' UT'

    index = itime(3)*60 + itime(4)

    print, 'Seeking index number : ', index

    Amie_Pot = reform(Amie_Data(index,  0, *, *))
    Amie_FAC = reform(Amie_Data(index, 18, *, *))
    Amie_Ped = reform(Amie_Data(index,  1, *, *))
    Amie_Hal = reform(Amie_Data(index,  3, *, *))
    Amie_Lats = Amie_Lats
    Amie_Lons = Amie_Mlts

    Iono_Pot = reform(Iono_Data(0, 5, *, *))
    Iono_FAC = reform(Iono_Data(0, 4, *, *))
    Iono_Ped = reform(Iono_Data(0, 3, *, *))
    Iono_Hal = reform(Iono_Data(0, 2, *, *))
    Iono_Lats = 90.0-reform(iono_data(0,0,0,*))
    Iono_Lons = (reform(iono_data(0,1,*,0))+180.0) mod 360

    shr = chopr('0'+tostr(itime(3)),2)
    smi = chopr('0'+tostr(itime(4)),2)
    con_file = 'con_'+shr+smi+'.ps'
    pot_file = 'pot_'+shr+smi+'.ps'

    setdevice, con_file, 'p', 5

    ppp = 4
    space = 0.05
    pos_space, ppp, space, sizes

    makect, 'mid'

    plotdumb
    plot_data, Amie_Ped, Amie_Lons, Amie_Lats, $
      Iono_Ped, Iono_Lons, Iono_Lats, $
      0, ppp, space, sizes, $
      Amie_Vars(1)

    plot_data, Amie_Hal, Amie_Lons, Amie_Lats, $
      Iono_Hal, Iono_Lons, Iono_Lats, $
      2, ppp, space, sizes, $
      Amie_Vars(3)

    get_position, ppp, space, sizes, 0, pos
    xyouts, pos(0)-0.03, (pos(1)+pos(3))/2.0, 'Pedersen Conductance', $
      alignment = 0.5, orient = 90, /norm
    xyouts, 0.5, pos(3)+0.05, stime, alignment = 0.5, charsize = 2.0, /norm
    get_position, ppp, space, sizes, 2, pos
    xyouts, pos(0)-0.03, (pos(1)+pos(3))/2.0, 'Hall Conductance', $
      alignment = 0.5, orient = 90, /norm

    closedevice

    setdevice, pot_file, 'p', 5

    plotdumb
    plot_data, Amie_Pot, Amie_Lons, Amie_Lats, $
      Iono_Pot, Iono_Lons, Iono_Lats, $
      0, ppp, space, sizes, $
      Amie_Vars(0)

    plot_data, Amie_FAC, Amie_Lons, Amie_Lats, $
      Iono_FAC, Iono_Lons, Iono_Lats, $
      2, ppp, space, sizes, $
      Amie_Vars(18)

    get_position, ppp, space, sizes, 0, pos
    xyouts, pos(0)-0.03, (pos(1)+pos(3))/2.0, 'Potential', $
      alignment = 0.5, orient = 90, /norm
    xyouts, 0.5, pos(3)+0.05, stime, alignment = 0.5, charsize = 2.0, /norm
    get_position, ppp, space, sizes, 2, pos
    xyouts, pos(0)-0.03, (pos(1)+pos(3))/2.0, 'Field-Aligned Current', $
      alignment = 0.5, orient = 90, /norm

    closedevice

endfor

end
