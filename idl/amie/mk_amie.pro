
variability = 0.0
cc = 0.0

ByBase =   0.0
BzBase =  -2.0
rBz = 0.0
rBy = 0.0
rH  = 0.0

filename = 'nightside_'+$
           tostr(BzBase)+'_'+$
           chopr('0'+tostr(rBz),2)+'_'+$
           chopr('0'+tostr(ByBase),2)+'_'+$
           chopr('0'+tostr(rBy),2)+'_'+$
           chopr('0'+tostr(rH),2)+'_'+$
           chopr('0'+tostr(100*variability),2)+'_'+$
           chopr('0'+tostr(100*cc),2)+'.bin'

filename = 'substorm_withoutpbi'+'.bin'

nTimes = 24*6*10*2
;nTimes = 96*4

iTime = [2000,3,20,0,0,0]
c_a_to_r, iTime, basetime

imf = fltarr(nTimes, 4)
ae = fltarr(nTimes, 4)
dst = fltarr(nTimes, 2)
hpi = fltarr(nTimes, 2)
cpcp = fltarr(nTimes)

; bx has no meaning, really.  Can be ignored

bx = fltarr(nTimes)+0.0 + 0.0*randomn(s,nTimes)

; by sets a circular potential at the pole.  You could set this to 
; zero if you want a dawn-dusk symmetric potential pattern.

by = fltarr(nTimes) + ByBase + rBy * randomn(s,nTimes)

; This sets the strength of the background two cell convection pattern

bz = fltarr(nTimes)+BzBase + rBz*randomu(s,nTimes)

; This sets the resolution:

lats = 90.0-findgen(81)/2
mlts = findgen(97)/4.0

; This sets where the potential and aurora maximizes (potential has
; max/min at the ocflb, while the aurora is equatorward of this.)

ocflbBase = 75.0 - 3.0*cos(mlts*!pi/12.0)

nLats = n_elements(lats)
nMlts = n_elements(mlts)

nVars = 3

data = fltarr(nTimes, nVars, nMlts, nLats)

potential = fltarr(nMlts, nLats)
eflux     = fltarr(nMlts, nLats)
avee      = fltarr(nMlts, nLats)
; area is used to calculate the hemispheric power
area      = fltarr(nMlts, nLats)
for iLat = 0,nLats-1 do $
   area(*,iLat) = 111.0*(111.0*360.0/nMlts) * cos(lats(iLat)*!dtor)*1000.0*1000.0

time = dindgen(nTimes)/nTimes * 2.0D*24.0D*3600.0D + basetime

; I have a substorm going off on the second day, starting at 11 UT.
; You can "turn this off" by setting the time after the run is
; finished (i.e. > 48 hours)

SubstormGrowth   = (24.0+11.0)*3600.0D + basetime
SubstormOnset    = (24.0+12.0)*3600.0D + basetime
SubstormRecovery = (24.0+12.5)*3600.0D + basetime
DtRecovery = 3600.0D

; This is stuff for the substorm

DeltaOCFLBOnset    =  -5.0
DeltaOCFLBRecovery =  10.0

for i=0,nTimes-1 do begin

   ocflb = ocflbBase

   fac_pt_base = 5.0
   fac_ef_base = 1.0
   fac_ae_base = 6.0

   ; This is the amplitude of the potential in the harang
   ; discontinuity.  I would set this to zero, or you could play with
   ; it for a strong jet....

   amp_h = 10.0 + randomn(s,1) * rH
   amp_h = amp_h(0)

   ; This is the amplitude of the potential in the PBIs
   ; I would set this to zero, or you could play with
   ; it for a strong jet....

   amp_pbi = 30.0 + randomn(s,1) * rH
   amp_pbi = amp_pbi(0)*1000.0

   ; This is the amplitude, multiplied by Bz for the strength of the CPCP
   amp    = 10.0*1000.0
   ; Amplitude of the E-flux
   amp_ef = 1.5                 ; /1000.0
   ; Amplitude of the average energy
   amp_ae = 2.5

   pbi = 0

;   ; Substorm growth - expand the potential pattern
;   if (time(i) gt SubstormOnset-3600.0*6 and $
;       time(i) le SubstormGrowth) then begin
;      Dt = (time(i)-SubstormGrowth)/(3600.0*6)
;
;      saw = dt*18.0 - floor(dt*18.0)
;
;      loctim = 24.0 - 2.0*cos(float(floor(dt*18.0))/4 * !pi)
;      ocflb_pbi = ocflb + 2.0 - 2.0*saw
;      print, loctim
;      pbi = 1
;
;      intensity = exp(-2*(saw-0.2)^2.0)
;
;; 20 minute lifecyle
;; 22 MLT and 02 MLT
;; pot ~ 10ish kV  (going for 100mV/m)
;; 3-5 keV
;; possibly 10 times the bckground e-flux
;; moves equatorward + eastward
;
;
;   endif

   pbi = 0

   ; Substorm growth - expand the potential pattern
   if (time(i) gt SubstormGrowth and $
       time(i) le SubstormOnset) then begin
      Dt = (time(i)-SubstormGrowth)/(SubstormOnset-SubstormGrowth)
      ocflb = ocflbBase + Dt*DeltaOCFLBOnset
      ; Change Bz over 10 minute
      Dt = (time(i)-SubstormGrowth)/600.0
      if (dt gt 1.0) then Dt = 1.0
      bz(i) = bz(i) - 3.0*Dt
   endif

   ; Onset - retract the auroral oval, and broaden it greatly and
   ; strengthen the precipitation

   if (time(i) gt SubstormOnset and $
       time(i) le SubstormRecovery) then begin
      Dt = (time(i)-SubstormOnset)/(SubstormRecovery-SubstormOnset)
      ocflb = ocflbBase + DeltaOCFLBOnset + $
              Dt*DeltaOCFLBRecovery*(2.0+cos(mlts*!pi/12.0))/3.0
      fac_ef_base = 1.0 + 2.0 * dt
      amp_h = amp_h + 20.0 * dt
      amp_ef = amp_ef + (10.0-amp_ef)*dt
      amp_ae = amp_ae + (5.0-amp_ae)*dt
      fac_pt_base = 5.0 + 10.0*dt
      ; Change Bz over 10 minutes
      Dt = (time(i)-SubstormOnset)/600.0
      if (dt gt 1.0) then Dt = 1.0
      bz(i) = bz(i) - 3.0*(1.0-Dt)
   endif

   ; recovery - go back to the original state

   if (time(i) gt SubstormRecovery and $
       time(i) le SubstormRecovery+DtRecovery) then begin
      Dt = (time(i)-SubstormRecovery)/(DtRecovery)
      ocflb = ocflbBase + (1.0-Dt)*(DeltaOCFLBOnset + $
              DeltaOCFLBRecovery*(2.0+cos(mlts*!pi/12.0))/3.0)
      fac_ef_base = 3.0 - 2.0 * dt
      amp_h = amp_h + 20.0 * (1.0-dt)
      amp_ef = amp_ef + (10.0-amp_ef)*(1.0-dt)
      amp_ae = amp_ae + (5.0-amp_ae)*(1.0-dt)
      fac_pt_base = 5.0 + 10.0*(1.0-dt)
   endif

   amp_h = amp_h * 1000.0

   m = (90.0-mean(ocflb))/2

   for iMlt = 0,nMlts-1 do begin

      fac = fltarr(nLats) + fac_pt_base/2
      fac_ef = fltarr(nLats) + 0.6*fac_ef_base * (0.5+1.5*(cos(mlts(iMlt)*!pi/12.0)+1.0))
      fac_ae = fltarr(nLats) + fac_ae_base

      d = lats-ocflb(iMlt)

      ; build the potential pattern

      l = where(d ge 0)
      fac(l) = fac_pt_base
      line = exp(-abs(lats-ocflb(iMlt))/fac)
      potential(iMlt,*) = -bz(i) * amp * line * sin(mlts(iMlt)*!pi/12.0)

      ; Some By
      line = -amp/2 * by(i) * exp(-(90.0-lats)/m)
      potential(iMlt,*) = potential(iMlt,*) + line

      ; Harang Discontinuity - centered at midnight
      line = -amp_h * exp(-2*(lats-ocflb(iMlt))^2/fac^2) * $
             ((cos((mlts(iMlt)-1)*!pi/12.0)+1.0)/2.0)^3.0
      l = where(line gt 0,c)
      if (c gt 0) then line(l) = 0.0
      potential(iMlt,*) = potential(iMlt,*) + line


      if (pbi) then begin
         ; PBIs - centered at midnight
         line = -intensity * amp_pbi * exp(-(lats-ocflb_pbi(iMlt))^2/1.5^2) * $
                exp(50.0*((cos((mlts(iMlt)-loctim)*!pi/12.0)-1.0)/2.0))
         l = where(line gt 0,c)
         if (c gt 0) then line(l) = 0.0
        potential(iMlt,*) = potential(iMlt,*) + line
;         potential(iMlt,*) = line
      endif

      ; Add some variability
;      fac(l) = 5.0
      vari_line  = $
         variability*randomn(s,nLats) * exp(-abs(lats-ocflb(iMlt))/fac)
      vari_line2 = $
         variability*randomn(s,nLats) * exp(-abs(lats-ocflb(iMlt))/fac)

      potential(iMlt,*) = potential(iMlt,*) + $
                          vari_line*potential(iMlt,*)

      ; build the energy flux pattern

      l = where(d lt 0)
      fac_ef(l) = fac_ef(l)*2
      line = exp(-abs(lats-ocflb(iMlt))^2/fac_ef^2)
      l = where(line gt 0.5,c)
      if (c gt 0) then line(l) = 1.0
      eflux(iMlt,*) = amp_ef * line * (cos(mlts(iMlt)*!pi/12.0)+2.0)
      eflux(iMlt,*) = eflux(iMlt,*) + $
                      cc * vari_line * eflux(iMlt,*) + $
                      sign(cc) * (1.0 - abs(cc)) * vari_line2 * eflux(iMlt,*)

      if (pbi) then begin
         ; PBIs
         line = intensity * 10.0 * exp(-(lats-ocflb_pbi(iMlt))^2/1.5^2) * $
                exp(50.0*((cos((mlts(iMlt)-loctim)*!pi/12.0)-1.0)/2.0))
         l = where(line lt 0,c)
         if (c gt 0) then line(l) = 0.0
         eflux(iMlt,*) = eflux(iMlt,*) + line
      endif

      ; build the averaged energy pattern

      line = exp(-(lats-ocflb(iMlt))^2/fac_ae^2)
      avee(iMlt,*) = amp_ae * line * (cos(mlts(iMlt)*!pi/12.0)+5.0)/6
      avee(iMlt,*) = avee(iMlt,*) + $
                     cc * vari_line * avee(iMlt,*) + $
                     sign(cc) * (1.0 - abs(cc)) * vari_line2 * avee(iMlt,*)

      if (pbi) then begin
         ; PBIs
         line = intensity * 3.0 * exp(-(lats-ocflb_pbi(iMlt))^2/1.5^2) * $
                exp(50.0*((cos((mlts(iMlt)-loctim)*!pi/12.0)-1.0)/2.0))
         l = where(line lt 0,c)
         if (c gt 0) then line(l) = 0.0
         avee(iMlt,*) = avee(iMlt,*) + line
      endif

   endfor

   ; set base values for aurora

   l = where(avee lt 0.5)
   avee(l) = 0.5

   l = where(eflux lt 0.1)
   eflux(l) = 0.1

   ; fill the data array

   data(i,0,*,*) = potential
   data(i,1,*,*) = eflux
   data(i,2,*,*) = avee

   ; calculate hemispheric power

   hpi(i,0) = total(area*eflux*0.001)/1.0e9
   print, 'hp: ', i, hpi(i,0), amp_h

   ; calculate cpcp

   cpcp(i) = max(potential)-min(potential)

endfor

imf(*,0) = 400.0
imf(*,1) = bx
imf(*,2) = by
imf(*,3) = bz

Version = 1.0

vars = ['Potential (kV)', 'Total Energy Flux (ergs/cm2/s)', 'Mean Energy (ergs)']

; write amie file

amie_write_binary, filename, Vars, lats, mlts, time, data, $
                   imf = imf, ae = ae, dst = dst, hpi = hpi, cpcp = cpcp, $
                   Version = Version

end

