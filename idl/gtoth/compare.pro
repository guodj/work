;==========================================
pro compare,w0,w1,wnames

; Compare all variables in w0 and w1 by calculating
; relative difference in the 1st norm.
;==========================================

on_error,2

sizew0=size(w0)
sizew1=size(w1)

if sizew0(0) ne sizew1(0) then begin
   print,'w0 and w1 have different dimensions:',sizew0(0),' and ',sizew1(0)
   retall
endif

ndim=sizew0(0)-1

if ndim eq 0 then begin
   ndim=1
   nw=1
endif else $
   nw=sizew0(ndim+1)

if max(abs(sizew0(1:ndim)-sizew1(1:ndim))) gt 0 then begin
   print,'w0 and w1 have different sizes:',sizew0(1:ndim),' /= ',sizew1(1:ndim)
   retall
endif

if keyword_set(wnames) then $
   print,'var max(|A-B|)/max(|A|+|B|) sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)' $
else $
   print,'ind max(|A-B|)/max(|A|+|B|) sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)'

for iw=0,nw-1 do begin
   case ndim of
   1: begin
      wsum=max(abs(w0(*,iw))+abs(w1(*,iw)))
      wdif=max(abs(w0(*,iw)-w1(*,iw)))
      wsum1=total(abs(w0(*,iw))+abs(w1(*,iw)))
      wdif1=total(abs(w0(*,iw)-w1(*,iw)))
      end
   2: begin
      wsum=max(abs(w0(*,*,iw))+abs(w1(*,*,iw)))
      wdif=max(abs(w0(*,*,iw)-w1(*,*,iw)))
      wsum1=total(abs(w0(*,*,iw))+abs(w1(*,*,iw)))
      wdif1=total(abs(w0(*,*,iw)-w1(*,*,iw)))
      end
   3: begin
      wsum=max(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))
      wdif=max(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
      wsum1=total(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))
      wdif1=total(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
      end
   endcase

   if keyword_set(wnames) then begin
      if wsum eq 0. then print,wnames(iw),' wsum=0' $
      else               print,wnames(iw),wdif/wsum,wdif1/wsum1,wsum
   endif else begin
      if wsum eq 0. then print,iw,' wsum=0' $
      else               print,iw,wdif/wsum,wdif1/wsum1,wsum
   endelse

endfor
end

