
;==========================================
pro getpict_asc,unit,npict,ndim,nw,nx,x,w
;==========================================
  on_error,2
  ;----------------------------------------
  ; Read coordinates and values row by row
  ;----------------------------------------
  wrow=dblarr(nw)
  xrow=dblarr(ndim)
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    x=dblarr(nx(0),ndim)
    w=dblarr(nx(0),nw)
    for x0=0,nx(0)-1 do begin
      readf,unit,xrow,wrow
      x(x0,0:ndim-1)=xrow(0:ndim-1)
      w(x0,0:nw-1)  =wrow(0:nw-1)
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    x=dblarr(nx(0),nx(1),ndim)
    w=dblarr(nx(0),nx(1),nw)
    for x1=0,nx(1)-1 do begin
      for x0=0,nx(0)-1 do begin
        readf,unit,xrow,wrow
        x(x0,x1,0:ndim-1)=xrow(0:ndim-1)
        w(x0,x1,0:nw-1)  =wrow(0:nw-1)
      endfor
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    x=dblarr(nx(0),nx(1),nx(2),ndim)
    w=dblarr(nx(0),nx(1),nx(2),nw)
    for x2=0,nx(2)-1 do begin
      for x1=0,nx(1)-1 do begin
        for x0=0,nx(0)-1 do begin
          readf,unit,xrow,wrow
          x(x0,x1,x2,0:ndim-1)=xrow(0:ndim-1)
          w(x0,x1,x2,0:nw-1)=wrow(0:nw-1)
        endfor
      endfor
    endfor
  end
  endcase
end

;==========================================
pro getpict_bin,unit,npict,ndim,nw,nx,x,w
;==========================================
  on_error,2
  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=dblarr(n1,ndim)
    w=dblarr(n1,nw)
    wi=dblarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=dblarr(n1,n2,ndim)
    w=dblarr(n1,n2,nw)
    wi=dblarr(n1,n2)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=dblarr(n1,n2,n3,ndim)
    w=dblarr(n1,n2,n3,nw)
    wi=dblarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;==========================================
pro getpict_real,unit,npict,ndim,nw,nx,x,w
;==========================================
  on_error,2
  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=fltarr(n1,ndim)
    w=fltarr(n1,nw)
    wi=fltarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=fltarr(n1,n2,ndim)
    w=fltarr(n1,n2,nw)
    readu,unit,x
    wi=fltarr(n1,n2)
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=fltarr(n1,n2,n3,ndim)
    w=fltarr(n1,n2,n3,nw)
    wi=fltarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end


;==========================================
pro getpict,unit,filetype,npict,x,w,$
    headline,physics,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,error
;==========================================

   on_error,2

   error=0

   if(eof(unit))then begin
      error=1
      return
   endif

   ; Get current pointer position
   point_lun,-unit,pointer

   ; Skip npict-1 snapshots
   ipict=0
   pictsize=1
   while ipict lt npict-1 and not eof(unit) do begin
      ipict=ipict+1
      gethead,unit,filetype,pictsize=pictsize
      pointer=pointer+pictsize
      point_lun,unit,pointer
   endwhile

   ; Backup 1 snapshot if end of file
   if eof(unit) then begin
       error=1
       point_lun,unit,pointer-pictsize
   endif

   ; Read header information
   gethead,unit,filetype,headline,physics,$
       it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables

   ; Read data
   case filetype of
   'ascii':  getpict_asc ,unit, npict, ndim, nw, nx, x, w
   'binary': getpict_bin ,unit, npict, ndim, nw, nx, x, w
   'real4':  getpict_real,unit, npict, ndim, nw, nx, x, w
    else:    begin
                print,'Getpict: unknown filetype:',filetype
                error=1
                close,unit
             end
   endcase

end

