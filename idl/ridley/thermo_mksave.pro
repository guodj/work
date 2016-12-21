alt = 500.0
v = sqrt(3.99e14/(6372000.0+alt*1000))
period = (6372000.0+alt*1000)*2*!pi/v
rps = 2*!pi/period

date = '130317'
filelist = findfile('armada??_'+date+'.bin')
nFiles = n_elements(filelist)

nSats = nFiles

iO_ = 4
iN2_ = 6
MaxValN2 = 1.0e21

for iSat = 0,nSats-1 do begin

   file = filelist(iSat)
   print,'reading file : ',file
   gitm_read_bin, file, data1, time1, nVars, Vars, version
   save, file=file+'.save',file, data1, time1, Vars

   Vars = [Vars,'O/N2']

   nT = n_elements(time1)
   on2 = fltarr(nT)

   a = data1(0,2,0,0,*)
   na = n_elements(a)
   da = a(2:na-2)-a(1:na-3)
   nas = n_elements(da)

   for iT=0,nT-1 do begin
      n2 = data1(iT,iN2_,0,0,*)
      o1 = data1(iT,iO_,0,0,*)

      mn2 = (n2(2:na-2)+n2(1:na-3))/2.0
      mo1 = (o1(2:na-2)+o1(1:na-3))/2.0
      
      iA = nas-1
      Done = 0
      n2int = 0.0
      o1int = 0.0
      while (not Done) do begin
         if (n2int + mn2(iA)*da(iA) lt MaxValN2) then begin
            n2int = n2int + mn2(iA)*da(iA)
            o1int = o1int + mo1(iA)*da(iA)
            iA = iA-1
         endif else begin
            dan = (MaxValN2 - n2int) / mn2(iA)
            n2int = n2int + mn2(iA)*dan
            o1int = o1int + mo1(iA)*dan
            on2(iT) = o1int/n2int
            Done = 1
         endelse
      endwhile
   endfor

   if (iSat eq 0) then begin

      alts = reform(data1(0,2,0,0,*))/1000.0
      l = where(alts gt alt)
      ia1 = l(0)-1
      r = (alt-alts(ia1))/(alts(ia1+1)-alts(ia1))

      dt = 2.0
      iShift = rps*dt
      stime = time1(0)

      tmax = max(time1-stime)
      nT = long(tmax/dt)
      t = dindgen(nT)*dt + stime

      vals = [1,4,5,6,15,16,24,35,36]
      nVals = n_elements(vals)

      data = fltarr(nSats, nVals+1, nT)

   endif

   print, 'interpolating...'
   for i=0,nVals-1 do begin

      iV = vals(i)
      dummy = reform(data1(*,iV,0,0,*))
      dummy1d = r*dummy(*,ia1+1) + (1.0-r)*dummy(*,ia1)
      data(iSat,i,*) = interpol(dummy1d,time1,t)
      
   endfor
   data(iSat,i,*) = interpol(on2,time1,t)

endfor

vals = [vals,nVars]
nVals = n_elements(vals)

save,file = 'all_data_1alt.save',data,t,vals,vars

end
