
;-----------------------------------------------------------------------------
;
; function geocgm
;   Geocgm calls the fortran subroutine which transforms between the
;   geographic coordinate system and the corrected geomagnetic
;   coordinate system.
;

function geocgm, year, lat, lon, alt, cgm_to_geo = cgm_to_geo

  if keyword_set(cgm_to_geo) eq 0 then cor = 1.0 else cor = -1.0

  yr = float(year)
  clar = 0.0		; cgm latitude
  clor = 0.0		; cgm longitude
  slar = 0.0		; geo latitude
  slor = 0.0		; geo longitude
  pla = 0.0		; cgm pole latitude
  plo = 0.0		; cgm longitude
  btr = 0.0		; IGRF Magnetic field H (nT)
  bfr = 0.0		; IGRF Magnetic field D (deg)
  brr = 0.0		; IGRF Magnetic Field Z (nT)
  dla = 0.0		; dipole latitude
  dlo = 0.0		; dipole longitude
  slac = 0.0		; conjugate geo latitude
  slor = 0.0		; conjugate geo longitude
  slaf = 0.0		; footprint geo latitude
  slof = 0.0		; footprint geo longitude
  rbm = 0.0		; Apex of magnetic field line (RE max)
  apx = 0.0		; Apex of magnetic field line h (km)
  ut = 0.0		; MLT Midnight in UT hours
  alf = 0.0		; Meridian angle (east +)

  if n_elements(alt) eq 0 then hi = 0.0 else hi = float(alt)
  if cor eq 1.0 then begin
    slar = float(lat)
    slor = float(lon)
  endif else begin
    clar = float(lat)
    clor = float(lon)
  endelse

  if (!version.os eq 'vms') or (!version.os eq 'VMS') then		$
    z = call_external('geocgm_exe','geocgm',				$
	cor, yr, clar, clor, hi, pla, plo, btr, bfr,			$
	brr, dla, dlo, slar, slor, slaf, slof, rbm,			$
	apx, ut, alf)

  info = {year:fix(yr),glat:slar, glon:slor, alt:hi, 			$
	plat:pla, plon:plo, bhor:btr, bdec:bfr, brad:brr, 		$
	dlat:dla, dlon:dlo, mlat:clar, mlon:clor, 			$
	fglat:slaf, fglon:slof, cglat:slac, cglon:slor,			$
	apexr:rbm, apex:apx, ut:ut, ma:alf}

;  info.year	year
;  info.mlat	cgm latitude
;  info.mlon	cgm longitude
;  info.glat	geo latitude
;  info.glon	geo longitude
;  info.plat 	cgm pole latitude
;  info.plon 	cgm longitude
;  info.bhor	IGRF Magnetic field H (nT)
;  info.bdec 	IGRF Magnetic field D (deg)
;  info.brad 	IGRF Magnetic Field Z (nT)
;  info.dlat 	dipole latitude
;  info.dlon 	dipole longitude
;  info.cglat 	conjugate geo latitude
;  info.cglon 	conjugate geo longitude
;  info.fglat 	footprint geo latitude
;  info.fglon 	footprint geo longitude
;  info.apexr 	Apex of magnetic field line (RE max)
;  info.apex 	Apex of magnetic field line h (km)
;  info.ut	MLT Midnight in UT hours
;  info.ma	Meridian angle (east +)

  return, info

end

