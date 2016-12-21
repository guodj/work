
pro load_ovation, model

  winter_path = 'Aurora_Coef_Ovation/winter/'
  spring_path = 'Aurora_Coef_Ovation/spring/'
  summer_path = 'Aurora_Coef_Ovation/summer/'
  fall_path =   'Aurora_Coef_Ovation/fall/'
  season_path = [winter_path,spring_path,summer_path,fall_path]
  sname_a = ['winter','spring','summer','fall']
  atype_string = ['diff','mono','wave','ions']
  aname = ['diff_','mono_','wave_','ions_']
  n_or_s = 3

;;;grid size of model (fixed)
  nmlt = 96
  nmlat = 160
  ndF = 12

  naTypes = 4
  njTypes = 2
  Prob_all = fltarr(4,njTypes,naTypes,ndF,nmlt,nmlat)
  b1p      = fltarr(4,njTypes,naTypes,nmlt,nmlat)
  b2p      = fltarr(4,njTypes,naTypes,nmlt,nmlat)
  b1a      = fltarr(4,njTypes,naTypes,nmlt,nmlat)
  b2a      = fltarr(4,njTypes,naTypes,nmlt,nmlat)

  lats = findgen(nmlat/2)*0.5 + 50.0
  mlts = findgen(nmlt)/nmlt * 24.0

  for jt = 0,njTypes-1 do begin 

     if (jt eq 0) then jtype = 1
     if (jt eq 1) then jtype = 3

     for atype = 0,naTypes-1 do begin 

        print, 'atype : ',atype

        for iseason=0,3 do begin
 
           sname = sname_a(iseason)
; afile = season_path(iseason) + sname + '_' +  atype_string(atype)
           afile = season_path(iseason) + atype_string(atype)
           if( (jtype eq 3) or (jtype eq 4) )then afile = afile + '_n'
           afile = afile + '.txt'
; pfile = season_path(iseason) + sname + '_' + 'prob_b_'
           pfile = season_path(iseason) + 'prob_b_'
           pfile = pfile + atype_string(atype) + '.txt'

;need to calculate probability functions from linear regressions
;read those in

           if( atype le 2 )then begin
              b1 = 0. 
              b2 = 0.
              yend = 1900
              dend = 1
              y0 = 1900
              d0 = 1
              files_done = 0
              sf0 = 0
              openr,20,pfile
              readf,20,y0,d0,yend,dend,files_done,sf0

              for i =0, nmlt-1 do begin
                 for j=0, nmlat-1 do begin
                    readf,20,b1,b2
                    b1p(iseason,jt,atype,i,j) = b1
                    b2p(iseason,jt,atype,i,j) = b2
                 endfor
              endfor
              p = 0.0
              for i=0,nmlt-1 do begin
                 for j=0,nmlat-1 do begin
                    for k=0,ndF-1 do begin
                       readf,20,p
                       prob_all(iseason,jt,atype,k,i,j) = p
;                       readf,20,Prob_all(iseason,jt,atype,0:ndF-1,i,j)
                    endfor
                 endfor
              endfor
              close,20
           endif

;now read in regression coefficients for auroral flux
           openr,20,afile 
           readf,20,y0,d0,yend,dend,files_done,sf0
           i = 0
           j = 0
           b1 = 0.
           b2 = 0.
           rFa = fltarr(nmlt,nmlat)
           while (not eof(20) )do begin
              readf,20,i,j,b1,b2,rF
              rFa(i,j) = rF
              b1a(iseason,jt,atype,i,j) = b1
              b2a(iseason,jt,atype,i,j) = b2
           endwhile
           close,20

;Prob = Prob_all(iseason,0:2,0:nmlt-1,0:nmlat-1)
;           Prob = fltarr(3,ndF,nmlt,nmlat)
;           for i1 = 0, 2 do begin
;              for i2=0,ndF-1 do begin
;                 for i3=0,nmlt-1 do begin
;                    for i4=0,nmlat-1 do begin
;                       Prob(i1,i2,i3,i4) = Prob_all(iseason,jt,i1,i2,i3,i4)
;                    endfor
;                 endfor
;              endfor
;           endfor
           
        endfor
        
     endfor
   
  endfor


  model = {PROB_ALL:Prob_all, B1P:b1p, B2P:b2p, B1A:b1a, B2A:b2a, $
          LATS: lats, MLTS: mlts}

end
