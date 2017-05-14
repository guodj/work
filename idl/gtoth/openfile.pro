;==========================================
pro openfile,unit,filename,filetype
;==========================================
   on_error,2

   close,unit
   case filetype of
       'ascii' :openr,unit,filename
       'binary':openr,unit,filename,/f77_unf
       'real4' :openr,unit,filename,/f77_unf
       else    :print,'Openfile: unknown filetype:',filetype
   endcase
end

