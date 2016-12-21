;^CFG COPYRIGHT VAC_UM
;=============================================================================
;
; Read the log data from 1, 2, or 3 files into wlog, wlog0, wlog1, wlog2 arrays
;
; Original version written by M. Nauta, 
; later extended by G. Toth, for VAC and then 
; improved and simplified for BATSRUS
; 
;=============================================================================

pro getlogpro, logfilename, $
               nlogfile, nVars, wlognames, wlog, wlog0, wlog1, wlog2

nlogfile=0

str2arr,logfilename,logfilenames,nlogfile
if nlogfile gt 3 then begin
   print,'Error in GetLog: cannot handle more than 3 files.'
   retall
endif

for ifile=0,nlogfile-1 do begin
   data=fstat(1)
   if data.open ne 0 then close,1
   openr,1,logfilenames(ifile)

   headline=''
   readf,1,headline
   print,'headline       =',strtrim(headline,2)
   wlognames=''
   readf,1,wlognames

   wlognames=str_sep(strtrim(strcompress(wlognames),2),' ')
   nVars=n_elements(wlognames)

   if nlogfile eq 1 then begin
      print,'Reading array wlog:'
      for i=0,nVars-1 do $
         print,FORMAT='("  wlog(*,",I2,")= ",A)',i,wlognames(i)
   endif else begin
      print,'Reading array wlog',ifile,FORMAT='(a,i1)'
      for i=0,nVars-1 do $
         print,FORMAT='("  wlog",I1,"(*,",I2,")= ",A)',ifile,i,wlognames(i)
   endelse

   buf=long(10000)
   dbuf=long(10000)
   wlog=dblarr(nVars,buf)
   wlog_=dblarr(nVars)
   nt=long(0)
   while not eof(1) do begin
      on_ioerror,close_file
      readf,1,wlog_
      wlog(*,nt)=wlog_
      nt=nt+1
      if nt ge buf then begin
         buf=buf+dbuf
         wlog=[[wlog],[dblarr(nVars,buf)]]
      endif
   endwhile
close_file:close,1
   print,'Number of recorded timesteps: nt=',nt
   wlog=transpose(wlog(*,0:nt-1))

   if nlogfile gt 1 then begin
     case ifile of
     0: wlog0=wlog
     1: wlog1=wlog
     2: wlog2=wlog
     endcase
   endif
endfor

end
