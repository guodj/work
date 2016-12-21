;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

initial_guess = findfile('-t b* v*')
iFile = 0
while (strpos(initial_guess(iFile),"data") gt 0 or $
       strpos(initial_guess(iFile),"ps") gt 0 or $
       strpos(initial_guess(iFile),"sum") gt 0) do iFile = iFile + 1
initial_guess = initial_guess(iFile)
test = findfile(initial_guess)
if strlen(test(0)) eq 0 then initial_guess=initial_guess+'*'

filelist = ask('AMIE binary file name',initial_guess)

amiefiles = findfile(filelist)
nfiles = n_elements(amiefiles)

for nf = 0, nfiles-1 do begin

   amie_file = amiefiles(nf)

   print, "reading ",amie_file

   read_amie_binary, amie_file, data, lats, mlts, time, fields, $
                     imf, ae, dst, hp, cpcp, version

   if (imf(0,3) gt 0.0) then imf(*,3) = -imf(*,3)

   openw,1,'imf_'+amie_file

   printf,1,''
   printf,1,'Taken from AMIE file:'
   printf,1,'  '+amie_file
   printf,1,''
   printf,1,'#TIMEDELAY'
   printf,1,'0.0'
   printf,1,''
   printf,1,'#START'

   for i =0,n_elements(time)-1 do begin
      c_r_to_a, itime, time(i)
      printf,1,format='(i5,6i3,8f9.2)', itime, 0, $
             imf(i,*),0.0,0.0,5.0,50000.0
   endfor

   close,1

endfor

end
