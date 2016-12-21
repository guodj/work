;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;	CNVCOORD
;
; PURPOSE:
;	Convert coordinates from geographic to PACE magnetic or
;	from PACE magnetic to geographic
;
;	calling sequence:
;	  pos = cnvcoord(inpos,[inlong],[height],[/GEO], [MODEL=year])
;	     the routine can be called either with a 3-element floating
;	     point array giving the input latitude, longitude and height
;	     or it can be called with 3 separate floating point values
;	     giving the same inputs.  The default conversion is from
;	     geographic to PACE geomagnetic coordinates.  If the keyword
;	     GEO is set (/GEO) then the conversion is from magnetic to
;	     geographic.  The optional keyword MODEL can be used
;            to select a particular magnetic coordinates model.
;            It should be set to the year of the model (for example,
;            MODEL=1990).  If the model is not specified, the most
;            recent model is used by default.  Once a model has
;            been specified it will continue to be used on subsequent
;            calls, even if the keyword is not included in the call.
;;
;	     The input array can also be given in the form 
;		inpos(3,d1,d2,. . .)
;	     The output array will be in the same form.
;----------------------------------------------------------------------------
function cnvcoord, in1, in2, in3, geo = geo, model = model
;
	if (keyword_set(geo)) then mgflag = 2 else mgflag = 1
        if (keyword_set(model)) then myear = model else myear = 1995
        myear = fix(myear)
        if (n_params() GE 3) then inp = float([in1,in2,in3]) $
		else inp = float(in1)
        if (n_elements(inp) MOD 3 NE 0) then begin
          print,'input position must be fltarr(3) [lat,long,height]'
          return,[0,0,0]
          end
	order=4
	err=0
	s0 = size(inp)
	tmp = reform(inp,3,n_elements(inp)/3)
	outpos=tmp
	s1 = size(tmp)
	lib_dir = getenv("IDL_EXTRAS")
	for j=0,s1(2)-1 do begin
	  invec = fltarr(3)
	  invec(0:2) = tmp(0:2,j)
	  outvec = invec
	  ret_val = call_external(lib_dir+'libpgm.so',			$
		'cnvcoord_idl',invec,order,outvec,mgflag,err,myear)
          if (ret_val NE 0) then begin 
            print,"cnvcoord error = ",ret_val,err
            return,0
            endif
	  outpos(*,j)=outvec
	endfor
	outpos=reform(outpos,s0(1:n_elements(s0)-3))
	return,outpos
	end

