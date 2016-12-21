
pro readapex, filename, apex

  openr,1,filename, /f77

  nLats = 179
  nLons = 73

  glats = fltarr(nLons, nLats)
  glons = fltarr(nLons, nLats)
  alats = fltarr(nLons, nLats)
  alons = fltarr(nLons, nLats)
  bm = fltarr(nLons, nLats)
  bx = fltarr(nLons, nLats)
  by = fltarr(nLons, nLats)
  bz = fltarr(nLons, nLats)

  GeoLat=0.0
  GeoLon=0.0
  alat=0.0
  alon=0.0
  bmag=0.0
  xmag=0.0
  ymag=0.0
  zmag=0.0

  for iLat = -89, 89 do begin
      for iLon = 0, 360, 5 do begin
          
          readu,1,GeoLat,GeoLon,alat,alon,bmag,xmag,ymag,zmag

          glats(iLon/5,iLat+89) =   GeoLat
          glons(iLon/5,iLat+89) =   GeoLon
          alats(iLon/5,iLat+89) =   alat
          alons(iLon/5,iLat+89) =   alon
          bm   (iLon/5,iLat+89) =   bmag
          bx   (iLon/5,iLat+89) =   xmag
          by   (iLon/5,iLat+89) =   ymag
          bz   (iLon/5,iLat+89) =   zmag

      endfor
  endfor

  close,1

  apex = create_struct(name = 'apexdata', $
                      'glats',   glats, $
                      'glons',   glons, $
                      'alats', alats, $
                      'alons', alons, $
                      'bmag', bm, $
                      'bx', bx, $
                      'by', by, $
                      'bz', bz)

end

;filename = 'apex2002.dat'
;readapex, filename, apex
;
;contour, apex.alats, levels = [-88,88]
;
;end
