;
;  pro radar_dmsp
;
;  displays radar and dmsp observations of ionospheric convection
;

pro radar_dmsp, single, t, placement, npts

; have to identify where radar data is and where dmsp data is
; radar data is located in posive positions  - placement(i,0) > 0
; satellite data is located in negative positions  - placement(i,0) < 0

; how many satellites and radars do we have :

  nradar = max(placement(*,0)) - 1
  nsat   = min(placement(*,0))
  if nsat ge 0 then nstat = -1 else nstat = abs(nstat)-1

; first plot satellite data :


