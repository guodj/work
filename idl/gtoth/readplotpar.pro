;===========================================================================
pro readplotpar,ndim,cut,cut0,plotdim,nfunc,func,funcs,funcs1,funcs2,$
   plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask
;===========================================================================
   on_error,2

   ; Determine dimension of plots based on cut or ndim,
   ; calculate reduced cut0 array by eliminating degenerate dimensions
   if keyword_set(cut) then begin
      cut0=reform2(cut)
      siz=size(cut0)
      plotdim=siz(0)
   endif else begin
      plotdim=ndim
      cut0=0
   endelse

   askstr,'func(s) (e.g. rho p ux;uz bx+by -T) ',func,doask
   if plotdim eq 1 then begin
      print,'1D plotmode: plot'
      plotmode='plot'
   endif else begin
      if plotmode eq 'plot' then plotmode=''
      print,'2D scalar: ',$
            'shade/surface/contour/contlabel/contfill/contbar/tv/tvbar'
      print,'2D polar : polar/polarlabel/polarfill/polarbar'
      print,'2D vector: stream/stream2/vector/velovect/ovelovect'
      askstr,'plotmode(s)                ',plotmode,doask
   endelse
   askstr,'plottitle(s) (e.g. B [G];J)',plottitle,doask
   askstr,'autorange(s) (y/n)         ',autorange,doask

   nfunc=0
   str2arr,func,funcs,nfunc
   str2arr,plotmode,plotmodes,nfunc
   str2arr,plottitle,plottitles,nfunc,';'
   str2arr,autorange,autoranges,nfunc

   funcs1=strarr(nfunc)
   funcs2=strarr(nfunc)
   for ifunc=0,nfunc-1 do begin
      func12=str_sep(funcs(ifunc),';')
      funcs1(ifunc)=func12(0)
      if n_elements(func12) eq 2 then funcs2(ifunc)=func12(1)
   endfor

end
