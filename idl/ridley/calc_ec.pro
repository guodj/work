;boro.pro function to evaluate Borovsky's function
;Patrick Newell, March 2008

pro boro,sintc,n,v,B,Eb

  if(n le 0. )then begin
     Eb = 0.
     return
  endif

;mp = 1.67e-24                   ;grams
  mp = 1.
  rhom = 0.
  rho = n*mp
  VA = (2.18e6)*B/sqrt(n)       ; in cm/s
  VA = VA/(1.0e5)               ; in km/s
  MA = v/VA
  Mms = MA
  C = ((0.25)^6. + (1./(1. + 1.38*alog(MA)))^6.)^(-0.166667)
  betaS = 0.032*(MA^1.92)
  fact1 = 1.6*sintc
  fact2 = rho*v*v
  fact3 = (1 + 0.5/(MMS*MMS))
  fact4 = 1./sqrt(1. + betaS)
  fact5 = (C*rho + (1./sqrt(1.+betaS))*rhom)^(-0.5)
  fact6 = (sqrt(1.+betaS)+1.)^(-0.5)
  R = fact1*fact2*fact3*fact4*fact5*fact6
  Eb = R

  return
end

pro sol_coup,bx,by,bz,v,n,Ec
  ncoup = 33

  Ec = fltarr(ncoup)
  B = sqrt(bx*bx + by*by + bz*bz)
  if( (n*B*v) eq 0. )then begin
     return
  endif
  pi = 3.14159265
  BT = sqrt(by*by + bz*bz)
  bztemp = bz
  if( bz eq 0 ) then begin
     bztemp = 0.001
  endif
  tc = atan(by,bztemp)
  if( BT*cos(tc)*bz lt 0 ) then begin
     tc = tc + 3.14159265
  endif
  sintc = abs(sin(tc/2.))
  Bs = bz
  if( Bs gt 0. )then Bs = 0.
  EKL = v*BT*sintc*sintc
  Ewav = EKL*sintc*sintc
;  tau = (dipole_tilt(day,utsec) - 90.)*3.14159265/180.
;  sintcp = abs(sin(tc/2. + tau))
;  exptc = (exp(abs(tc)/pi) - 1.)/1.71828

  Ec[0] = bz
  Ec[1] = EKL
  eps = v*B*B*sintc*sintc*sintc*sintc
;Ec(2) = eps
  Ec[2] = v*sqrt(BT)*sintc*sintc ;Lyatsky
  Ec[3] = v*Bs
  p = n*v*v*(1.67e-6)
  Ec[4] = p                                   ;dynamic pressure in nPa
  Ec[5] =sintc*sintc*EKL*(v^0.3333)*(p^0.16667) ;Vasylinas
  Ec[6] =sintc*sintc*BT*(p^0.33333)             ;Another Vasyliunas formula
  Ec(7) = Ewav
;Ec[8] = n
  Ec[8] = sqrt(n)*v*v           ;8/2007 change for twovar paper plot
  Ec[9] = v
  Ec[10] = Bs
  Ec[11] = v*BT
  Ec[12] = v*BT*BT*sintc*sintc*sintc*sintc
  Ec[13] = v*B*sintc*sintc*sintc*sintc
  Ec[14] = Ewav*Ewav
  Ec[15] = sqrt(Ewav)
  Ec[16] = Ewav^0.6667
  Ec[17] = sqrt(Ekl)
  Ec[18] = EKL*(p^0.166667)*(v^0.3333)
  Ec[19] = (v^1.33333)*(sintc^2.66667)*(BT^0.66667)
;Ec[19] = (v^1.33333)*exptc*(BT^0.66667)
  Ec[20]  = Ec(19)/(v^0.33333)
  Ec[21] = (v^0.166667)*Ec(19)
  Ec[22] = (sintc^0.333333)*Ec(19)
  Ec[23] = 0.
  if( sintc ne 0. )then Ec[23] = Ec[19]/(sintc^0.66667)
  Ec[24] = (BT^0.3333)*Ec[19]
  Ec[25] = 0.
  if( BT gt 0.)then Ec[25] = Ec[19]/(BT^0.166667)
  Ec[26] = Ewav*sqrt(p)
  Ec[27] = v*(BT^0.67)*sintc*sintc ;Weimer
  Ec[28] = sqrt(n)*v*v*BT*sintc*sintc*sintc*sintc*sintc*sintc
  Ec[29] = Ec[19]*sqrt(p)
  Ec[30] = Ec[19]*(n^0.16667)
;Ec[31] = sintc*(n^(-0.22))*(v^0.56)*(B^1.44)       ;Borovosky
  boro,sintc,n,v,B,Eb
  Ec[31] = Eb
  Ec[32] = Ec[31]*(sintc^1.667) ;variant on Borovosky
;plot, [-1,0,1],[-1,0,-1]

  return
end

pro calc_ec, time, mag, vel, den, Ec

  nPts = n_elements(time)

  bx   = fltarr(nPts)
  by   = fltarr(nPts)
  bz   = fltarr(nPts)
  v    = fltarr(nPts)
  ni   = fltarr(nPts)
  Ec   = fltarr(nPts)
  EcUnSmoothed = fltarr(nPts)

  ; This is from Pat...
  nh = 1.0                       ;hours previous to integrate over
  wh = 0.65                    ;reduce weighting by factor of wh each hour back

  for iT = 0, nPts-1 do begin
     sol_coup,mag(0,iT),mag(1,iT),mag(2,iT),abs(vel(0,iT)),den(iT),Ec_single
     EcUnSmoothed(iT) = Ec_single(19)
  endfor

  for iT = 0, nPts-1 do begin

     l = where(time le time(iT) and time gt time(iT)-nh*3600.0, c)
     w = wh^(abs(time(iT)-time(l))/3600.0)
     tw = mean(w)
     bx(iT) = mean(mag(0,l)*w)/tw
     by(iT) = mean(mag(1,l)*w)/tw
     bz(iT) = mean(mag(2,l)*w)/tw
     v(iT)  = abs(mean(vel(0,l)*w)/tw)
     ni(iT) = mean(den(l)*w)/tw
     Ec(iT) = mean(EcUnSmoothed(l)*w)/tw

  endfor

end
