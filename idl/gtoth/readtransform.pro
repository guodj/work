;===========================================================================
pro readtransform,ndim,nx,gencoord,transform,nxreg,xreglimits,wregpad,$
    physics,nvector,vectors,grid,doask
;===========================================================================
   on_error,2

   if (gencoord or transform eq 'unpolar') and ndim eq 2 then begin
      if transform eq '' then begin
        transform='none'
        askstr,"transform (r=regular/p=polar/u=unpolar/n=none)",transform,1
      endif else $
        askstr,"transform (r=regular/p=polar/u=unpolar/n=none)",transform,doask

      ; Complete name
      case transform of
          'r': transform='regular'
          'p': transform='polar'
          'u': transform='unpolar'
          'n': transform='none'
         else:
      endcase
      ; Get transformation parameters and calculate grid
      case 1 of
        transform eq 'regular':begin
           print,'Generalized coordinates, dimensions for regular grid'
           if n_elements(nxreg) ne 2 then nxreg=[0,0]
           if n_elements(xreglimits) ne 4 then xreglimits=dblarr(4) $
           else xreglimits=double(xreglimits)
           nxreg0=nxreg(0)
           nxreg1=nxreg(1)
           asknum,'nxreg(0) (use negative sign to limit x)',nxreg0,doask
           if nxreg0 lt 0 then begin
               nxreg0=abs(nxreg0)
               xmin=0 & xmax=0
               asknum,'xreglimits(0) (xmin)',xmin,doask
               asknum,'xreglimits(2) (xmax)',xmax,doask
               xreglimits(0)=xmin
               xreglimits(2)=xmax
           endif
           asknum,'nxreg(1) (use negative sign to limit y)',nxreg1,doask
           if nxreg1 lt 0 then begin
               nxreg1=abs(nxreg1)
               ymin=0 & ymax=0
               asknum,'xreglimits(1) (ymin)',ymin,doask
               asknum,'xreglimits(3) (ymax)',ymax,doask
               xreglimits(1)=ymin
               xreglimits(3)=ymax
           endif
           grid=lindgen(nxreg0,nxreg1)
           nxreg=[nxreg0,nxreg1]
           wregpad=0
        end
        transform eq 'polar' or transform eq 'unpolar':begin
           getvectors,physics,nvector,vectors
           grid=lindgen(nx(0),nx(1))
        end
        transform eq 'none':grid=lindgen(nx(0),nx(1))
        else: print,'Unknown value for transform:',transform
      endcase
   endif else if gencoord and ndim eq 3 then begin
      if transform eq '' then begin
         transform="none" & askstr,"transform (s=sphere/n=none)",transform,1
      endif else $
         askstr,"transform (s=sphere/n=none)",transform,doask
      case transform of
         's': transform='sphere'
         'n': transform='none'
        else:
      endcase
      if transform eq 'sphere' then getvectors,physics,nvector,vectors
      grid=lindgen(nx(0),nx(1),nx(2))
   endif else case ndim of
      1: grid=lindgen(nx(0))
      2: grid=lindgen(nx(0),nx(1))
      3: grid=lindgen(nx(0),nx(1),nx(2))
   endcase

   ;===== GRID HELPS TO CREATE A CUT, E.G.: cut=grid(*,4)

   help,grid
end

