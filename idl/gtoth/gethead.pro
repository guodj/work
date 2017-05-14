;=============================================================
pro gethead,unit,filetype,headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nx,eqpar,variables,pictsize=pictsize
;=============================================================
   on_error,2

;Type definitions
   headline='                                                                               '
   it=long(1)
   ndim=long(1)
   neqpar=long(1)
   nw=long(1)
   varname='                                                                               '
;Remember pointer position at beginning of header
   point_lun,-unit,pointer0
;Read header
   case filetype of
      'ascii': begin
                  time=double(1)
                  readf,unit,headline
                  readf,unit,it,time,ndim,neqpar,nw
                  gencoord=(ndim lt 0)
                  ndim=abs(ndim)
                  nx=lonarr(ndim)
                  readf,unit,nx
                  eqpar=dblarr(neqpar)
                  readf,unit,eqpar
                  readf,unit,varname
               end
      'binary':begin
                  time=double(1)
                  readu,unit,headline
                  readu,unit,it,time,ndim,neqpar,nw
                  gencoord=(ndim lt 0)
                  ndim=abs(ndim)
                  nx=lonarr(ndim)
                  readu,unit,nx
                  eqpar=dblarr(neqpar)
                  readu,unit,eqpar
                  readu,unit,varname
               end
      'real4': begin
                  time=float(1)
                  readu,unit,headline
                  readu,unit,it,time,ndim,neqpar,nw
                  gencoord=(ndim lt 0)
                  ndim=abs(ndim)
                  nx=lonarr(ndim)
                  readu,unit,nx
                  eqpar=fltarr(neqpar)
                  readu,unit,eqpar
                  readu,unit,varname
               end
      else: begin
                  print,'Gethead: unknown filetype',filetype
                  retall
            end
   endcase

   if keyword_set(pictsize) then begin
      ; Calculate the picture size
      ; Header length
      point_lun,-unit,pointer1
      headlen=pointer1-pointer0
      ; Number of cells
      nxs=1
      for idim=1,ndim do nxs=nxs*nx(idim-1)
      ; Snapshot size = header + data + recordmarks
      case filetype of
         'ascii' :pictsize = headlen + (18*(ndim+nw)+1)*nxs
         'binary':pictsize = headlen + 8*(1+nw)+8*(ndim+nw)*nxs
         'real4' :pictsize = headlen + 8*(1+nw)+4*(ndim+nw)*nxs
      endcase
   endif else begin
      ; Get variables and physics
      variables=str_sep(strtrim(strcompress(varname),2),' ')
      tmp=str_sep(strtrim(headline,2),'_')
      if n_elements(tmp) eq 2 then begin
         headline=tmp(0)
         physics=tmp(1)
      endif
   endelse
end

