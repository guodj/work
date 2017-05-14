;==========================================
pro gettype,filenames,filetypes,npictinfiles
;==========================================
   on_error,2

   filetypes=filenames
   npictinfiles=intarr(n_elements(filenames))
   for ifile=0,n_elements(filenames)-1 do begin
      ; Obtain filetype based on the length info in the first 4 bytes
      close,10
      openr,10,filenames(ifile)
      len=long(1)
      readu,10,len
      if len ne 79 then ftype='ascii' else begin
         ; The length of the 2nd line decides between real4 and binary
         ; since it contains the time, which is real*8 or real*4
         head=bytarr(79+4)
         readu,10,head,len
         case len of
            20: ftype='real4'
            24: ftype='binary'
            else: begin
               print,'Error in GetType: strange unformatted file:',$
                     filenames(ifile)
               retall
            end
         endcase
      endelse
      close,10

      ; Obtain file size and number of snapshots
      openfile,1,filenames(ifile),ftype
      status=fstat(1)
      fsize=status.size

      pointer=0
      pictsize=1
      npict=0
      while pointer lt fsize do begin
          ; Obtain size of a single snapshot
          point_lun,1,pointer
          gethead,1,ftype,pictsize=pictsize
          npict=npict+1
          pointer=pointer+pictsize
      endwhile
      close,1

      npictinfiles(ifile)=npict
      filetypes(ifile)   =ftype
   endfor
end

