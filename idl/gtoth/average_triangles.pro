;==============================================================================

pro average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
        triangles,triangconj

   wconj=dblarr(nxreg(0),nxreg(1))

   for iw=1,nw-1 do begin
      ; Calculate wreg with original triangulation
      wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
            [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

      ; Calculate wconj with conjugated triangulation
      wconj       =trigrid(xx,yy,w(*,*,iw),triangconj, $
            [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

      wreg(*,*,iw) = 0.5*(wreg(*,*,iw) + wconj)
   endfor
   print,'Using averaged conjugated triangulation'

end
