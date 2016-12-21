;===========================================================================
pro getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,fmax,fmin,doask,$
                x,w,xreg,wreg,usereg,physics,eqpar,wnames,cut
;===========================================================================
   on_error,2

   for ifunc=0,nfunc-1 do begin
      if autoranges(ifunc) eq 'n' then begin
         if first then begin
            f_min=fmin(ifunc)
            f_max=fmax(ifunc)
            asknum,'Min value for '+funcs(ifunc),f_min,doask
            asknum,'Max value for '+funcs(ifunc),f_max,doask
            fmin(ifunc)=f_min
            fmax(ifunc)=f_max
         endif
      endif else begin
         if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            xreg,wreg,physics,eqpar,wnames,cut $
         else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            x,   w,   physics,eqpar,wnames,cut

         f_max=max(f)
         f_min=min(f)
         if first then begin
            fmax(ifunc)=f_max
            fmin(ifunc)=f_min
         endif else begin
            if f_max gt fmax(ifunc) then fmax(ifunc)=f_max
            if f_min lt fmin(ifunc) then fmin(ifunc)=f_min
         endelse
      endelse
   endfor
end

