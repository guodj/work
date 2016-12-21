
syear = '2004'

for im = 1,12 do begin

  smonth = chopr('0'+tostr(im),2)

  for id = 1,31 do begin
  
    sday = chopr('0'+tostr(id),2)
  
    dir = '/r/ridley/Data/cdaweb.gsfc.nasa.gov/pub/istp/ace/'
    
    mfi = dir+'mfi_h0/'+syear+'/ac_h0_mfi_'+syear+smonth+sday+'*.cdf'
    swe = dir+'swe_h0/'+syear+'/ac_h0_swe_'+syear+smonth+sday+'*.cdf'
    
    magfile = ''
    mfile = findfile(mfi)
    print, mfile
    if (strlen(mfile(0)) gt 0) then magfile = mfile(0)
    
    swfile = ''
    sfile = findfile(swe)
    print, sfile
    if (strlen(sfile(0)) gt 0) then swfile = sfile(0)
    
    if (strlen(magfile) gt 0 and strlen(swfile) gt 0) then begin
    
      read_cdf_imf, magfile, swfile, magfield, sw, magtime, swtime, isgse, delay
    
      ntimes = n_elements(magtime)
      ntsw   = n_elements(swtime)
    
      ; we are going to clean up the sw data here
  
      swclean = fltarr(5,ntsw)
      swtclean = dblarr(ntsw)
  
      npts=0
      for i=0,ntsw-1 do begin
  	if (max(abs(sw(*,i))) lt 1.0e10) then begin
  	  swclean(*,npts) = sw(0:4,i)
  	  swtclean(npts) = swtime(i)
  	  npts = npts + 1
  	endif
      endfor

      if (npts gt 0) then begin
        ntsw = npts
        sw = swclean(*,0:npts-1)
        swtime = swtclean(0:npts-1)
      endif else begin
        sw(0,*) = -400.0
        sw(1,*) =    0.0
        sw(2,*) =    0.0
        sw(3,*) =    5.0
        sw(4,*) = 100000.0
      endelse
  
      ; clean magfield data
  
      loc = where(abs(magfield(0,*)) lt 200.0, ntimes)
      magfield = magfield(*,loc)
      magtime  = magtime(loc)
  
      ; now linear interpolate the sw data into the magtime
  
      swnew = fltarr(5,ntimes)
      swtnew = dblarr(ntimes)
      
      for i=0,ntimes-1 do begin
      
  	d = abs(swtime - magtime(i))
  	loc = where(d eq min(d))
      
  	if (loc(0) eq 0) then begin
  	  swnew(*,i) = sw(0:4,0)
  	  swtnew(i) = swtime(0)
  	endif else begin
  	  if (loc(0) eq ntsw-1) then begin
  	    swnew(*,i) = sw(0:4,ntsw-1)
  	    swtnew(i) = swtime(ntsw-1)
  	  endif else begin
      
  	    j = loc(0)
  	    if (magtime(i)-swtime(j) lt 0) then j = j - 1
      
  	    x = (magtime(i) - swtime(j)) / (swtime(j+1) - swtime(j))
      
  	    swtnew(i) = (1.0-x) * swtime(j) + x * swtime(j+1)
  	    swnew(*,i) = (1.0-x) * sw(0:4,j) + x * sw(0:4,j+1)
      
  	  endelse
  	endelse
      
  	cdf_epoch, magtime(i), year, month, day, hour, minute, second, milli, /break
  	itime = [year,month,day,hour,minute,second]
  	c_a_to_r, itime, rtime
  	magtime(i) = rtime
  
      endfor
    
      b = magfield
      v = swnew(0:2,*)
      n = reform(swnew(3,*))
      t = reform(swnew(4,*))
      time = magtime
  
      bt = sqrt(b(1,*)^2+b(2,*)^2)
      theta = acos(b(2,*)/(bt+0.01))
  
      ekl = -v(0,*) * bt * (sin(theta/2.0)^2) * 1000.0 * 1.0e-9 * 1000.0
  
      mu0 = 4.0*!pi*1.0e-7
      k   = 1.3807e-23
      mp  = 1.6726e-27
  
      bmag = sqrt(b(0,*)^2+b(1,*)^2+b(2,*)^2)
  
      pdyn = n * 100.0^3 * mp * (v(0,*)*1000.0)^2
      pb   = 2.0*(bmag*1.0e-9)^2/mu0
      va = bmag*1.0e-9 / sqrt(mu0*n*100.0^3*mp)
      ma = -v(0,*)/(va/1000.0)
  
      save, b,v,n,t,time,isgse,ekl,bmag,pdyn,pb,va,ma,delay,$
  	    file='ace'+syear+smonth+sday+'.save'
  
    endif
  
  endfor

endfor

end
