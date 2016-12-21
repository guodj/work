
pro read_themis_fgm, cdffilename, nptsCdf, cdftime, cdfbfield

  line = ''
  openr,1,cdffilename
  readf,1,line
  readf,1,line

  nptsCdf = 0L
  cdftime = dblarr(100000L)
  cdfbfield = fltarr(3,100000L)

  itime = intarr(6)

  while not (eof(1)) do begin

      readf,1,line

      itime(0) = fix(strmid(line, 0,4))
      itime(1) = fix(strmid(line, 5,2))
      itime(2) = fix(strmid(line, 8,2))
      itime(3) = fix(strmid(line,11,2))
      itime(4) = fix(strmid(line,14,2))
      itime(5) = fix(strmid(line,17,2))
      c_a_to_r, itime, rtime
      cdftime(nptsCdf) = rtime

      cdfbfield(0,nptsCdf) = float(strmid(line, 23,14))
      cdfbfield(1,nptsCdf) = float(strmid(line, 37,14))
      cdfbfield(2,nptsCdf) = float(strmid(line, 51,14))

      nptsCdf = nptsCdf + 1L

  endwhile

  cdftime = cdftime(0:nptsCdf-1)
  cdfbfield = cdfbfield(*,0:nptsCdf-1)

end
