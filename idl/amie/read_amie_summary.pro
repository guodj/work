
pro read_amie_summary, file, data

  line = ''
  close,1
  openr,1,file

  readf,1,line
  if (strpos(line,'nMag') gt 0) then begin
      mags = 1
      tmp = fltarr(20)
  endif else begin
      mags = 0
      tmp = fltarr(18)
  endelse

  nPtsMax = 1440
  time    = dblarr(nPtsMax)
  cpcp    = fltarr(nPtsMax)
  ae      = fltarr(3,nPtsMax)
  imfsw   = fltarr(4,nPtsMax)
  dst     = fltarr(nPtsMax)
  hpi     = fltarr(nPtsMax)
  sjh     = fltarr(nPtsMax)
  area    = fltarr(nPtsMax)
  emax    = fltarr(nPtsMax)
  nmags   = intarr(nPtsMax)

  itime = intarr(6)
  nPts = 0
  while (not eof(1)) do begin

      readf,1,tmp
      itime(*) = tmp(0:5)
      c_a_to_r, itime, rtime

      time(nPts) = rtime

      cpcp(nPts)    = tmp(6)
      ae(*,nPts)    = tmp(7:9)
      imfsw(*,nPts) = tmp(10:13)
      dst(nPts)     = tmp(14)
      hpi(nPts)     = tmp(15)
      sjh(nPts)     = tmp(16)
      area(nPts)    = tmp(17)
      if (mags) then begin
          nmags(nPts) = tmp(19)
          emax(nPts)    = tmp(18)
      endif

      nPts = nPts + 1

      if (nPts gt nPtsMax) then begin
          print, 'Reading in to many points....'
          print, 'Reset nPtsMax in read_amie_summary.pro'
          stop
      endif

  endwhile

  data = create_struct(name = 'AMIESummary', $
                       'time' , time, $
                       'cpcp' , cpcp, $
                       'ae'   , ae, $
                       'imfsw', imfsw, $
                       'dst'  , dst, $
                       'hpi'  , hpi, $
                       'sjh'  , sjh, $
                       'area' , area, $
                       'emax' , emax, $
                       'nMags', nMags)

  close,1

end
