;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;	MLT
;
; PURPOSE:
;
;	convert UT time to MLT
;
;	calling sequence:
;	  mt = mlt(year, ut_seconds, mag_long)
;		inputs:  year, time in seconds from Jan. 1 at 00;00:00 UT
;		         magnetic longitude of the observation point
;	
;	        the time in seconds can be found with routine "cnvtime"
;		the magnetic longitude of the obs. point can be found
;		  by using routine "cnvcoord"
;
;-----------------------------------------------------------------------------
function mlt, year, t, mlong
	year = fix(year)
	t = long(t)
	mlong = float(mlong)
	lib_dir = getenv("IDL_EXTRAS")
	if (n_elements(t) EQ 1) then mt=0.0 else mt=fltarr(n_elements(t))
	if (n_elements(t) EQ 1) then $
          mt = call_external(lib_dir+'libpgm.so',			$
		'mlt_idl', year, t, mlong,/f_value) 			$
	else 								$
	  for i=0,n_elements(t)-1 do 					$
	     mt(i)=call_external(lib_dir+'libpgm.so',			$
	       'mlt_idl', year, t(i), mlong, /f_value)

	return, mt
	end

