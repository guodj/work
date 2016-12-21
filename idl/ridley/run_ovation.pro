
;*****************************************************

function dF_bin, dF_call
  ndF = 12
  dFave = 4421.
  dFstep = dFave/8.
  dF = dF_call
  bin = fix(dF/dFstep)
  if( bin lt 0 )then bin = 0
  if( bin gt (ndF-1) )then bin = ndF-1
  return,bin
end

;****************************************************

function prob_estimate,b1,b2,Prob,dF,atype,i,j
;common prob_array,Prob
  ndF = 12

  p = b1 + b2*dF
  if( p lt 0. )then p = 0.
  if( p gt 1. )then p = 1.

  if( (b1 ne 0.) or (b2 ne 0.) )then begin
     return,p
  endif

  dFbin = dF_bin(dF)
  p0 = Prob[dFbin]
  if( p0 ne 0. )then begin
     return,p0
  endif

  dFbin_m1 = dFbin -1
  dFbin_p1 = dFbin + 1
  if( dFbin_m1 lt 0 )then dFbin_m1 = dFbin + 2
  if( dFbin_p1 gt (ndF-1) )then dFbin_p1 = dFbin - 2
  p1 = (Prob[dFbin_m1] + Prob[dFbin_p1])/2.
  return,p1

end

;***************************************************************

function j_plot_min,atype,jtype
  min = 0.
  if( jtype eq 3 )then min = 1.0e7
  if( jtype eq 4 )then min = 1.0e6
  if( jtype eq 5 )then min = 100.
  if( jtype eq 6 )then min = 100.
  if( jtype eq 7 )then min = -0.4
  if( jtype eq 8 )then min = 0.4
  return,min
end

;***************************************************************


function j_plot_max,atype,jtype
  if( jtype eq 1 )then begin    ;energy flux (e-)
     max = 1.5
     if( atype eq 0 )then max = 1.5 ;diff aurora e- flux ergs/cm2 s
     if( atype eq 1 )then max = 0.75 ;mono e- flux
     if( atype eq 2 )then max = 0.50 ;wave
     return,max
  endif
  if( jtype eq 2 )then begin    ;ion energy flux
     max = 0.5
     return,max
  endif
  if( jtype eq 3 )then begin    ;e- number flux
     max = 5.0e8
     if( atype eq 0 )then max = 5.0e8 ;diff aurora
     if( atype eq 1 )then max = 3.0e8 ;mono e- flux
     if( atype eq 2 )then max = 3.0e8 ;wave
     return,max
  endif
  if( jtype eq 4)then begin     ;ion number flux
     max = 3.0e7
     return,max
  endif
  if( (jtype eq 5) or (jtype eq 6) )then begin
     max = 30000.
     return,max
  endif
  if( jtype ge 7 )then begin
     max = 0.8
     return,max
  endif

  max = 0.
  return,max
end

;*******************************************************************

function mlt_bin, mlt_call
  mlt = mlt_call
  if(mlt lt 0.) then mlt = mlt + 24.
  if( mlt gt 24.) then mlt = mlt - 24.
  if( (mlt lt 0.) or (mlt gt 24.) )then begin
     return, -1
  endif

;bin = fix(mlt*2)
  bin = fix(mlt*4)
;if(bin gt 47) then bin =47
  if( bin gt 95 )then bin =95
  if(bin lt 0) then bin = 0
  return,bin
end

;*******************************************************************

function mlat_bin, mlat_call

  mlat = mlat_call
  if( mlat lt 0. )then begin
     mlat = -mlat
     if( (mlat lt 50.) or (mlat gt 90.) ) then begin
        return, -1
     endif
     bin = fix( (mlat-50.)*2 )
     if(bin lt 0) then bin = 0
     if(bin gt 79) then bin = 79
     return, bin
  endif

  if( mlat gt 0. )then begin
     if( (mlat lt 50.) or (mlat gt 90.) ) then begin
        return, -1
     endif
     bin = fix( (mlat-50.)*2 )
     if(bin lt 0) then bin = 0
     if(bin gt 79) then bin = 79
     bin = bin + 80
     return, bin
  endif
end

;*******************************************************************

function area_sphere,mlt1,mlt2,mlat1,mlat2
;area described by solid angle on Earth returned in km2
  Re = 6371.                    ;Earth radius in km
  delta_phi = (2.*!pi/24.)*abs(mlt2-mlt1)
  theta1 = (90. - mlat1)*!pi/180.
  theta2 = (90. - mlat2)*!pi/180.
  A = delta_phi*Re*Re*abs(cos(theta1)-cos(theta2))
  return,A
end

;*******************************************************************

pro season_weights,doy,winter_w,spring_w,summer_w,fall_w

  winter_w = 0.
  spring_w = 0.
  summer_w = 0.
  fall_w = 0.

  if( (doy ge 79) and (doy lt 171) )then begin
     summer_w = 1. - float(171-doy)/92.
     spring_w = 1. - summer_w
     return
  endif

  if( (doy ge 171) and (doy lt 263) )then begin
     fall_w = 1. - float(263-doy)/92.
     summer_w = 1. - fall_w
     return
  endif

  if( (doy ge 263) and (doy lt 354) )then begin
     winter_w = 1. - float(354-doy)/91.
     fall_w = 1. - winter_w
     return
  endif

; must be in the range 354 to 78.  Convert 354-365 to negative numbers
  doy0 = doy
  if( doy ge 354 )then doy0 = doy - 365
  spring_w = 1. - float(79-doy0)/90.
  winter_w = 1. - spring_w
  return

end

pro get_specific, model, iSeason, at, jt, ec, je

  nmlt  = n_elements(model.b1p(iseason,jt,at,*,0))
  nmlat = n_elements(model.b1p(iseason,jt,at,0,*))

  prob_curr = fltarr(nmlt,nmlat)
  if( at le 2 )then begin
     for i=0,nmlt-1 do begin
        for j=0,nmlat-1 do begin
           b1t = model.b1p(iseason,jt,at,i,j)
           b2t = model.b2p(iseason,jt,at,i,j)
           p = reform(model.Prob_all(iseason,jt,at,*,i,j))
           prob_curr(i,j) = prob_estimate(b1t,b2t,p,ec,at,i,j)
        endfor
     endfor
  endif

  if( at ge 3 )then begin
     for i=0,nmlt-1 do begin
        for j=0,nmlat-1 do begin
           prob_curr(i,j) = 1.
        endfor
     endfor
  endif

  je = fltarr(nmlt,nmlat)
  for i=0,nmlt-1 do begin
     for j=0,nmlat-1 do begin
        if( prob_curr(i,j) gt 0. )then begin
           je(i,j) = $
              (Ec*model.b2a(iseason,jt,at,i,j) + $
                  model.b1a(iseason,jt,at,i,j))*prob_curr(i,j)

           if( je(i,j) lt 0. )then je(i,j) = 0.
     

; These are screwed up because jt is NOT jtype....

;           if( (at le 2) and (jt eq 1) )then begin
;              if( je(i,j) gt 10. )then je(i,j) = 0.5
;              if( je(i,j) gt 5. )then je(i,j) = 5.
;           endif
 
           if( (at le 2) and (jt eq 1) )then begin
              if( je(i,j) lt 0. )then je(i,j) = 0.
              if( je(i,j) gt 2.0e10 )then je(i,j) = 0.
              if( je(i,j) gt 2.0e9 )then je(i,j) = 1.0e9
           endif

;           if( (at eq 3) and (jt eq 2) )then begin
;              if( je(i,j) lt 0. )then je(i,j) = 0.
;              if( je(i,j) gt 4. )then je(i,j) = 0.25
;              if( je(i,j) gt 2. )then je(i,j) = 2.
;           endif
;
;           if( (at eq 3) and (jt eq 4) )then begin
;              if( je(i,j) lt 0. )then je(i,j) = 0.
;              if( je(i,j) gt 5.0e8 )then je(i,j) = 0.
;              if( je(i,j) gt 1.0e8 )then je(i,j) = 1.0e8
;           endif
        endif
     endfor
  endfor

end

pro combine, jein, jeout

  nmlt  = n_elements(jein(*,0))
  nmlat = n_elements(jein(0,*))/2

  jeout = fltarr(nmlt, nmlat)

  for i = 0, nmlt-1 do for j = 0, nmlat-1 do begin
     j1 = jein(i,j)
     j2 = jein(i,j+nmlat)
     je = (j1+j2)/2.0
     if (j1*j2 eq 0) then je = j1 + j2
     jeout(i,j) = je
  endfor

end


pro run_ovation, model, ec, time, north, south

  c_r_to_a, itime, time
  doy = jday(itime(0), itime(1), itime(2))
  season_weights,doy,w0n,w1n,w2n,w3n
  doy = (doy + 182) mod 365
  season_weights,doy,w0s,w1s,w2s,w3s

  w0n = 0.25
  w1n = 0.25
  w2n = 0.25
  w3n = 0.25

  w0s = 0.25
  w1s = 0.25
  w2s = 0.25
  w3s = 0.25

;call with atype=0 for diffuse,  atype=1=mono,  2=wave,  3=ions
;jtype=1 for electron energy flux; 2=ion energy flux; 3=e- number flux
;jtype=4=ion number flux;  5=e- average energy  6=ion average energy
;jtype=7 = correlation coefficient, 8=corr coeff ions

  ; BUT - the code is broken, with only 2 jtypes working:
  ; jt = 0 => jtype = 1
  ; jt = 1 => jtype = 3

  ; what we really want right now is:
  ; diffuse (atype = 0) e- energy flux (jtype = 1; jt = 0)
  ; diffuse (atype = 0) e- number flux (jtype = 3; jt = 1)

  nmlt  = n_elements(model.b1p(0,0,0,*,0))
  nmlat = n_elements(model.b1p(0,0,0,0,*))

  je_all = fltarr(4, nmlt+1, nmlat)

  atype = 0
  jt = 0
  for iseason=0,3 do begin
     get_specific, model, iSeason, atype, jt, ec, je
     je_all(iSeason, 0:nmlt-1, *) = je
     je_all(iSeason, nmlt, *) = je(0,*)
  endfor

  je = w0n*je_all(0,*,*)+w1n*je_all(1,*,*)+w2n*je_all(2,*,*)+w3n*je_all(3,*,*)

  combine, reform(je), eflux

  ; Let's add some mono

  atype = 1
  jt = 0
  for iseason=0,3 do begin
     get_specific, model, iSeason, atype, jt, ec, je
     je_all(iSeason, 0:nmlt-1, *) = je
     je_all(iSeason, nmlt, *) = je(0,*)
  endfor

  je = w0n*je_all(0,*,*)+w1n*je_all(1,*,*)+w2n*je_all(2,*,*)+w3n*je_all(3,*,*)

  combine, reform(je), eflux_mono

  ; nFlux

  atype = 0
  jt = 1
  for iseason=0,3 do begin
     get_specific, model, iSeason, atype, jt, ec, je
     je_all(iSeason, 0:nmlt-1, *) = je
     je_all(iSeason, nmlt, *) = je(0,*)
  endfor
  je = w0n*je_all(0,*,*)+w1n*je_all(1,*,*)+w2n*je_all(2,*,*)+w3n*je_all(3,*,*)
  combine, reform(je), nflux

  ; Let's add some mono

  atype = 1
  jt = 1
  for iseason=0,3 do begin
     get_specific, model, iSeason, atype, jt, ec, je
     je_all(iSeason, 0:nmlt-1, *) = je
     je_all(iSeason, nmlt, *) = je(0,*)
  endfor
  je = w0n*je_all(0,*,*)+w1n*je_all(1,*,*)+w2n*je_all(2,*,*)+w3n*je_all(3,*,*)
  combine, reform(je), nflux_mono

  l = where(eflux_mono gt eflux and nflux_mono gt nflux, c)
  print, c
  if (c gt 0) then begin
     eflux(l) = eflux_mono(l)
     nflux(l) = nflux_mono(l)
  endif

  ; eflux is in ergs/cm2/s
  ; nflux is in part/cm2/s
  ; ergs to eV is 6.24150974e11

  l_low = where(nflux lt 0.01*median(nflux), c_low)
  if (c_low gt 0) then begin
     eflux(l_low) = 0.001*median(eflux)
     nflux(l_low) = 0.001*median(nflux)
  endif

  avee = eflux
  avee = avee/nflux * 6.24150974e11/1000.0
  if (c_low gt 0) then avee(l_low) = 0.1

  l_high = where(avee gt 25.0*median(avee), c_high)
  if (c_high gt 0) then avee(l_high) = 25.0*median(avee)
     
  north = {EFLUX:eflux, NFLUX:nflux, AVEE:avee, $
           LATS:model.lats,MLTS:[model.mlts,24.0]}

  je_all = fltarr(4, nmlt+1, nmlat)

  atype = 0
  jt = 0
  for iseason=0,3 do begin
     get_specific, model, iSeason, atype, jt, ec, je
     je_all(iSeason, 0:nmlt-1, *) = je
     je_all(iSeason, nmlt, *) = je(0,*)
  endfor

  je = w0s*je_all(0,*,*)+w1s*je_all(1,*,*)+w2s*je_all(2,*,*)+w3s*je_all(3,*,*)

  combine, reform(je), eflux

  atype = 0
  jt = 1
  for iseason=0,3 do begin
     get_specific, model, iSeason, atype, jt, ec, je
     je_all(iSeason, 0:nmlt-1, *) = je
     je_all(iSeason, nmlt, *) = je(0,*)
  endfor
  je = w0s*je_all(0,*,*)+w1s*je_all(1,*,*)+w2s*je_all(2,*,*)+w3s*je_all(3,*,*)
  combine, reform(je), nflux

  ; eflux is in ergs/cm2/s
  ; nflux is in part/cm2/s
  ; ergs to eV is 6.24150974e11

  l_low = where(nflux lt 0.01*median(nflux), c_low)
  if (c_low gt 0) then begin
     eflux(l_low) = 0.001*median(eflux)
     nflux(l_low) = 0.001*median(nflux)
  endif

  avee = eflux
  avee = avee/nflux * 6.24150974e11/1000.0
  if (c_low gt 0) then avee(l_low) = 0.1
     
  l_high = where(avee gt 25.0*median(avee), c_high)
  if (c_high gt 0) then avee(l_high) = 25.0*median(avee)
     
  south = {EFLUX:eflux, NFLUX:nflux, AVEE:avee, $
           LATS:model.lats,MLTS:[model.mlts,24.0]}

end
