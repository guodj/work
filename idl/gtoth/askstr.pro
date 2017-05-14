
;==========================================
pro askstr,prompt,var,doask
;==========================================
   on_error,2

   if var eq '' then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+var+' ? ',tmp
         if tmp ne '' then var=tmp
      endif else print,prompt,'=',var
   endelse
end

