;===========================================================================
pro readlimits,nfunc,funcs,autoranges,noautorange,fmax,fmin,doask
;===========================================================================
   on_error,2

   if n_elements(fmax) ne nfunc then fmax=dblarr(nfunc) else fmax=double(fmax)
   if n_elements(fmin) ne nfunc then fmin=dblarr(nfunc) else fmin=double(fmin)

   ; check if there is any function for which autorange is 'y'
   noautorange=1
   for ifunc=0,nfunc-1 do noautorange=noautorange and autoranges(ifunc) eq 'n'

   if(noautorange)then begin
      for ifunc=0,nfunc-1 do begin
         f_min=fmin(ifunc)
         f_max=fmax(ifunc)
         asknum,'Min value for '+funcs(ifunc),f_min,doask
         asknum,'Max value for '+funcs(ifunc),f_max,doask
         fmin(ifunc)=f_min
         fmax(ifunc)=f_max
      endfor
   endif

end
