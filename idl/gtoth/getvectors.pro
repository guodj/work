;===========================================================================
pro getvectors,physics,nvector,vectors
;===========================================================================
   physic=strtrim(physics)
   phys=strmid(physic,0,strlen(physic)-2)
   ndir=0
   reads,strmid(physic,strlen(physic)-1,1),ndir
   case phys of
   'rho':nvector=0
   'flx':nvector=0
   'hd' :begin
         nvector=1
         vectors=1
         end
   'hdadiab':begin
         nvector=1
         vectors=1
         end
   'mhdiso':begin
         nvector=2
         vectors=[1,ndir+1]
         end
   'mhd':begin
         nvector=2
         vectors=[1,ndir+2]
         end
   else:begin
      if nvector eq 0 then begin
         print,'Unrecognised physics: ',physics
         print,'Vector variables to transform for WREG'
         asknum,'nvector',nvector,doask
         if nvector gt 0 then begin
            vectors=intarr(nvector)
            read,PROMPT='Indices of first components in w? ',vectors
         endif
      endif
      end
   endcase
end

