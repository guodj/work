;===========================================================================
pro getfunc,f,f1,f2,func1,func2,x,w,physics,eqpar,wnames,cut
;===========================================================================
on_error,2

f1=funcdef(x,w,func1,physics,eqpar,wnames)

if keyword_set(cut) then f1=f1(cut)

if func2 eq '' then f=f1 else begin

   f2=funcdef(x,w,func2,physics,eqpar,wnames)

   if keyword_set(cut) then f2=f2(cut)

   ; Calculate f=sqrt(f1^2+f2^2)
   f=sqrt(f1^2+f2^2)
endelse

end

