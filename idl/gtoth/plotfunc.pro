;===========================================================================
pro plotfunc,x,w,xreg,wreg,usereg,ndim,physics,eqpar,rBody,$
  variables,wnames,axistype,plotmodes,plottitles,$
  ax,az,contourlevel,linestyle,$
  velvector,velspeed,velseed,velpos,velx,vely,veltri,$
  cut,cut0,plotdim,$
  nfunc,multix,multiy,fixaspect,plotix,plotiy,funcs,funcs1,funcs2,fmin,fmax,f
;===========================================================================
   on_error,2

   ; Get grid dimensions and set irr=1 if it is an irregular grid

   if keyword_set(cut) then siz = size(cut)  $
   else if usereg then      siz = size(xreg) $
   else                     siz = size(x)
   nx=siz(1)
   if plotdim eq 1 then begin
      ny=1
      irr=0
   endif else begin
      ny=siz(2)
      irr=ny eq 1
   endelse

   if irr and axistype eq 'cells' then begin
      print,'Irregular grid, axistype must be set to coord'
      axistype='coord'
   endif

   if axistype eq 'coord' then begin
      if usereg then $
         getaxes,ndim,xreg,xx,yy,zz,cut,cut0,rSlice,plotdim,variables $
      else $
         getaxes,ndim,x   ,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
   endif

   ; Calculate plot spacing from number of plots per page (ppp) and charsize
   if !p.charsize eq 0.0 then !p.charsize=1.0
   ppp   = multix*multiy
   space = max([float(!d.y_ch_size)/float(!d.y_size),$
                float(!d.x_ch_size)/float(!d.x_size)])*3.0*!p.charsize
   set_space, ppp, space, sizes, nx = multix, ny = multiy

   ; Store x and y titles
   xtitle = !x.title
   ytitle = !y.title

   for ifunc=0,nfunc-1 do begin

      !p.title=plottitles(ifunc)
      if !p.title eq 'default' then !p.title=funcs(ifunc)

      plotmod=plotmodes(ifunc)

      i=strpos(plotmod,'grid')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showgrid=1
      endif else showgrid=0

      i=strpos(plotmod,'mesh')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showgrid=1
          if irr then showmesh=0 else showmesh=1
      endif else showmesh=0

      i=strpos(plotmod,'body')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showbody=1
      endif else showbody=0

      i=strpos(plotmod,'fill')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          fill=1
      endif else fill=0

      i=strpos(plotmod,'bar')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
          showbar=1
          fill=1
      endif else showbar=0

      i=strpos(plotmod,'irr')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
          irr=1
      endif

      i=strpos(plotmod,'label')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+5)
          label=1
      endif else label=0

      ; contour --> cont
      i=strpos(plotmod,'contour')
      if i ge 0 then plotmod=strmid(plotmod,0,i+4)+strmid(plotmod,i+7)

      ; Calculate the next p.multi(0) explicitly
      if !p.multi(0) gt 0 then multi0=!p.multi(0)-1 $
      else multi0=!p.multi(1)*!p.multi(2)-1

      ; Calculate subplot position indices
      if !p.multi(4) then begin
         ; columnwise
         plotix=multix-1-multi0/multiy
         plotiy=multi0 mod multiy
      endif else begin
         ; rowwise
         plotix=multix-1-(multi0 mod multix)
         plotiy=multi0/multix
      endelse

      if plotmod ne 'shade' and plotmod ne 'surface' then begin

        ; obtain position for flat plotmodes
        set_position, sizes, plotix, multiy-1-plotiy, pos, /rect

        ; shrink in X direction for a colorbar in any previous plot
        if strpos(plotmodes(ifunc mod ppp),'bar') ge 0 then $
          pos(2) = pos(2) - (pos(2) - pos(0))*0.15

        ; shrink in X direction for the Y axis of plot
        if plotmod eq 'plot' and multix gt 1 then $
          pos(0) = pos(0) + (pos(2) - pos(0))*0.15

        if keyword_set(fixaspect) and plotmod ne 'plot' then begin

	  if plotmod eq 'polar' then $
            aspectx=1 $
          else begin
            if !x.range(1) ne !x.range(0) then    $
               width=abs(!x.range(1)-!x.range(0)) $
            else if axistype eq 'coord' then      $
               width=  max(xx) - min(xx)          $
            else                                  $
               width=  nx-1.0

            if !y.range(1) ne !y.range(0) then    $
               height=abs(!y.range(1)-!y.range(0))$
            else if axistype eq 'coord' then      $
               height= max(yy) - min(yy)          $
            else                                  $
               height= ny-1.0

            aspectx = width/height
          endelse

          aspectpos = (pos(2)-pos(0))/(pos(3)-pos(1)) $
                    *float(!d.x_size)/float(!d.y_size)

          aspectratio = aspectpos/aspectx

          ;print,'aspectx,pos,ratio=',aspectx,aspectpos,aspectratio

          if aspectratio gt 1 then begin
             posmid=(pos(2)+pos(0))/2.
             posdif=(pos(2)-pos(0))/2.
             pos(0)=posmid - posdif/aspectratio
             pos(2)=posmid + posdif/aspectratio
          endif else begin
             posmid=(pos(3)+pos(1))/2.
             posdif=(pos(3)-pos(1))/2.
             pos(1)=posmid - posdif*aspectratio
             pos(3)=posmid + posdif*aspectratio
          endelse
        endif

        ; Omit X axis if unneeded
        if (plotiy gt 0) then begin
          !x.tickname = strarr(60)+' '
          !x.title = ' '
        endif

        ; Omit Y axis if unneeded
        if (plotix gt 0 and plotmod ne 'plot') then begin
          !y.tickname = strarr(60)+' '
          !y.title = ' '
        endif

        !p.position = pos

      endif

      if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             xreg,wreg,physics,eqpar,wnames,cut0 $
      else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             x,  w,   physics,eqpar,wnames,cut0

      f_min=fmin(ifunc)
      f_max=fmax(ifunc)
      if f_max eq f_min then begin
         f_max=f_max+1
         f_min=f_min-1
      endif

      if plotmod eq 'plot' then $
         if nfunc gt ppp                then lstyle=ifunc/ppp $
         else if keyword_set(linestyle) then lstyle=linestyle $
         else                                lstyle=!p.linestyle

      ; Skip minimum ad maximum levels
      if plotmod eq 'cont' or plotmod eq 'polar' then $
         levels=(findgen(contourlevel+2)-1)/(contourlevel-1) $
                *(f_max-f_min)+f_min

      if plotmod eq 'tv' then begin
         ; Calculate plotting position and size

         tvplotx=pos(0)*!d.x_size
         tvploty=pos(1)*!d.y_size
         tvsizex=(pos(2)-pos(0))*!d.x_size
         tvsizey=(pos(3)-pos(1))*!d.y_size
         ; recalculate f for tv mode
         if !d.name eq 'PS' then tvf=congrid(f,200,200) $
         else                    tvf=congrid(f,tvsizex,tvsizey)

         tvf=bytscl(tvf,MIN=f_min,MAX=f_max,TOP=!D.TABLE_SIZE-3)+1
      endif

      case axistype of
      'cells': case plotmod of
         'cont': contour,f>f_min,LEVELS=levels,$
                 FILL=fill,FOLLOW=label,XSTYLE=1,YSTYLE=1,/NOERASE
         'plot'     :plot,f,YRANGE=[f_min,f_max],XSTYLE=18,ystyle=18, $
                                                 LINE=lstyle,/NOERASE
         'shade'    :begin
                        shade_surf,f>f_min,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                        if showgrid then $
                           surface,f>f_min,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                     end
         'surface'  :surface,f>f_min,ZRANGE=[f_min,f_max],$
                        XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :begin
                        tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
                        contour,f,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
                     end
         'vel'      :vector,f1,f2,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'vector'   :vector,f1,f2,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'stream'   :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=6,LENGTH=0.06,HEAD=0.1,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'stream2' :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=6,LENGTH=0.012,HEAD=0.5,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'velovect' :velovect,f1,f2,/NOERASE
         'ovelovect':velovect,f1,f2,/NOERASE,$
            XRANGE=[0,n_elements(f1(*,0))-1],YRANGE=[0,n_elements(f1(0,*))-1]
         endcase
      'coord': case plotmod of
         'cont'     :if irr then begin
                       if not keyword_set(tri) then triangulate,xx,yy,tri
                       contour,f>f_min,xx,yy,$
                          FOLLOW=label, FILL=fill, TRIANGULATION=tri, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
                    endif else $
                       contour,f>f_min,xx,yy,$
                          FOLLOW=label, FILL=fill, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
	 'polar'    :polar_contour,f>f_min,yy*!pi/180,xx,$
                          FOLLOW=label, FILL=fill, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
         'plot'     :plot,xx,f,YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,$
                          LINE=lstyle,/NOERASE
         'shade'    :if irr then begin
                        shade_surf_irr,f>f_min,xx,yy,AX=ax,AZ=az
                        shade_surf,f>f_min,xx,yy,AX=ax,AZ=az,/NODATA,/NOERASE
                     endif else begin
                        shade_surf,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                        if showgrid then $
                           surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                     endelse
         'surface'  :surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                        XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :begin
                       tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
                       contour,f,xx,yy,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
                     end
         'vel'      :vector,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                        TRIANGLES=veltri,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'vector'   :vector,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                        TRIANGLES=veltri,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'stream'   :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=6,LENGTH=0.06,HEAD=0.1,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'stream2' :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=6,LENGTH=0.012,HEAD=0.5,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'velovect' :velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE
         'ovelovect':velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE,$
                        XRANGE=[min(xx),max(xx)],YRANGE=[min(yy),max(yy)]
         endcase
      else:print,'Unknown axistype:',axistype
      endcase

      if showbody and axistype eq 'coord' then $
      if rBody gt abs(rSlice) then begin
         theta = findgen(37)*!pi*2.0/36.0
         rBodySlice=sqrt(rBody^2-rSlice^2)
         polyfill, rBodySlice*cos(theta), rBodySlice*sin(theta), color = 0
         ; redraw box in case the body is at the edge
         plot,xx,yy,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
      endif

      if showbar then $
         plotct, [pos(2)+0.005, pos(1), pos(2)+0.025, pos(3)], [f_min,f_max]

      if showgrid and plotdim eq 2 and plotmod ne 'surface'    $
                                   and plotmod ne 'shade' then begin
         if(plotmod eq 'polar')then                                       $
           plotgrid,xx,yy*!pi/180,lines=showmesh,xstyle=1,ystyle=1,/polar $
         else if keyword_set(cut) then                                    $
             plotgrid,xx,yy,lines=showmesh,xstyle=1,ystyle=1,polar=polar  $
         else                                                             $
              plotgrid,x,lines=showmesh,xstyle=1,ystyle=1
      endif

      !p.multi(0) = multi0
      !p.position = 0
      !x.title    = xtitle
      !x.tickname = strarr(60)
      !y.title    = ytitle
      !y.tickname = strarr(60)
   endfor

   !p.position = 0

end
