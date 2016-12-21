
filelist = findfile('*.cdf')
nFiles = n_elements(filelist)

nFiles = 10

time = dblarr(nFiles)

for iFile=0,nFiles-1 do begin

   print, iFile

   file = filelist(iFile)

   id = cdf_open(file)
   r  = cdf_inquire(id)

   cdf_varget, id, 'theta', theta, /zvariable
   cdf_varget, id, 'phi', phi, /zvariable
   cdf_varget, id, 'ep', potential, /zvariable
   cdf_varget, id, 'eflux', eflux, /zvariable
   cdf_varget, id, 'eave', avee, /zvariable

   cdf_close, id

   nTheta = n_elements(theta)
   nPhi   = n_elements(phi)

   theta2d = fltarr(nPhi, nTheta)
   phi2d = fltarr(nPhi, nTheta)

   pot2d = fltarr(nPhi, nTheta)
   eflux2d = fltarr(nPhi, nTheta)
   avee2d = fltarr(nPhi, nTheta)

   i = 0
   for iP = 0, nPhi-1 do for iT = 0, nTheta-1 do begin
      pot2d(iP,iT) = potential(i)*1000.0
      eflux2d(iP,iT) = eflux(i)*1000.0
      avee2d(iP,iT) = avee(i)
      i++
   endfor

   if (iFile eq 0) then begin
      nThetaHalf = (nTheta-1)/2
      dataNorth = fltarr(nFiles, 3, nPhi, nThetaHalf)
      dataSouth = fltarr(nFiles, 3, nPhi, nThetaHalf)
      ThetaSouth = -theta(0:nThetaHalf-1)
      ThetaNorth = reverse(theta(nThetaHalf+1:nTheta-1))
      cpcpn = fltarr(nFiles)
      cpcps = fltarr(nFiles)
   endif

   l = strpos(file,'.cdf')
   itime = [ $
           fix(strmid(file,l-19,4)), $
           fix(strmid(file,l-15,2)), $
           fix(strmid(file,l-13,2)), $
           fix(strmid(file,l-10,2)), $
           fix(strmid(file,l-8,2)), $
           fix(strmid(file,l-6,2))]
   c_a_to_r, itime, rtime

   time(iFile) = rtime

   dataSouth(iFile,0,*,*) = pot2d(*,0:nThetaHalf-1)
   dataSouth(iFile,1,*,*) = eflux2d(*,0:nThetaHalf-1)
   dataSouth(iFile,2,*,*) = avee2d(*,0:nThetaHalf-1)

   for i=0, nPhi-1 do begin
      dataNorth(iFile,0,i,*) = reverse(reform(pot2d(i,nThetaHalf+1:nTheta-1)))
      dataNorth(iFile,1,i,*) = reverse(reform(eflux2d(i,nThetaHalf+1:nTheta-1)))
      dataNorth(iFile,2,i,*) = reverse(reform(avee2d(i,nThetaHalf+1:nTheta-1)))
   endfor

   cpcps(iFile) = max(pot2d(*,0:nThetaHalf-1)) - min(pot2d(*,0:nThetaHalf-1))
   cpcpn(iFile) = max(pot2d(*,nThetaHalf+1:nTheta-1)) - $
                  min(pot2d(*,nThetaHalf+1:nTheta-1))

endfor

Vars = ['Potential (V)', $
        'Energy Flux (ergs/cm2)','Mean Energy (eV)']
Version = 0.1
mlts = (phi+180.0)/360.0*24.0
l = where(mlts gt 24.01,c)
if (c gt 0) then mlts(l) = mlts(l) mod 24.0

fileout = 'swmf_north.bin'
lats = ThetaNorth
amie_write_binary, fileout, Vars, lats, mlts, time, dataNorth, $
                   cpcp = cpcpn, Version = Version

fileout = 'swmf_south.bin'
lats = ThetaSouth
amie_write_binary, fileout, Vars, lats, mlts, time, dataSouth, $
                   cpcp = cpcps, Version = Version


end

