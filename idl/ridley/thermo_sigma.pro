
filelist = findfile('3DALL_t980106_??00*.bin')
nFiles = n_elements(filelist)

op_  = 24
o2p_ = 25
n2p_ = 26
np_  = 27
nop_ = 28
e_   = 33

setdevice, 'pedersen.ps', 'p', 5

ppp = 12
space = 0.05
pos_space, ppp, space, sizes

for iFile = 0, nFiles-1 do begin

   file = filelist(iFile)
   print, 'Reading file : ',file

   if (iFile mod ppp eq 0) then plotdumb

   gitm_read_bin, file, data, GitmTime, nVars, Vars, version

   s = size(data)
   nLons = s(2)
   nLats = s(3)
   nAlts = s(4)

   numden = fltarr(nlons, nlats, nalts)
   B0     = fltarr(nlons, nlats, nalts)

   lat_save = reform(data(1,0,*,0))

   for i=3,8 do numden = numden + reform(data(i,*,*,*))
   mmm = 0.0
   mass = [16.0,28.0,32.0,30.0,14.0,14.0]
   for i=3,8 do mmm = mmm + mass(i-3) * reform(data(i,*,*,*)) / numden

   altitude = reform(data( 2,*,*,*))
   tn       = reform(data(15,*,*,*))
   electron = reform(data(e_,*,*,*))
   te       = reform(data(e_+1,*,*,*))

   ec = 1.602e-19
   e2 = ec ^ 2
   mi = mmm * 1.6726e-27        ; pretend that mass ions = mass neutrals
   me = 9.1094e-31

   Vi = 2.6e-15 * (numden + electron)*(mmm^(-0.5))
   Ve = 5.4e-16 * (numden)*(TE^0.5)

   MeVe = me * ve
   MiVi = mi * vi

   B0_1d = 30000.0e-9 * (1.0 + 3.0*sin(lat_save*!dtor)^2)^0.5

   for i=0,nlons-1 do for k=0,nalts-1 do B0(i,*,k) = B0_1d

   GyroFrequency_Ion = ec*B0/Mi
   GyroFrequency_Electron = ec*B0/me

   VeOe = Ve^2 + GyroFrequency_Electron^2
   ViOi = Vi^2 + GyroFrequency_Ion^2


   Sigma_Pedersen = ((1.0/MeVe) * (Ve*Ve/VeOe) + $
                     (1.0/MiVi) * (Vi*Vi/ViOi)) * electron * e2

   Sigma_Hall = ((1.0/MeVe) * (Ve*GyroFrequency_Electron/VeOe) - $
                 (1.0/MiVi) * (Vi*GyroFrequency_Ion/ViOi)) * electron * e2

   Cond_Pedersen = fltarr(nlons,nlats)
   Cond_Hall     = fltarr(nlons,nlats)
   
   for k=0,nalts-2 do begin
      da = reform(altitude(*,*,k+1) - altitude(*,*,k))
      cond_pedersen(*,*) = cond_pedersen(*,*) + da * Sigma_Pedersen(*,*,k)
      cond_hall(*,*)     = cond_hall(*,*)     + da * Sigma_Hall(*,*,k)
   endfor

   lons = reform(data(0,*,0,0))/!dtor
   lats = reform(data(1,0,*,0))/!dtor

   value = cond_pedersen

;   determine_min_max, value, mini, maxi
   mini = 0.0
   maxi = 40.0
   nLevels = 31
   maxrange = 40.0

   print, maxi

   get_position, ppp, space, sizes, iFile mod ppp, pos
   
   makect,'mid'
   
   c_r_to_a, itime, GitmTime

   ut = itime(3) + itime(4)/60.0 + itime(5)/3600.0
   utrot = ut * 15.0

   if (iFile eq nFiles-1) then no06=1 else no06=0

   title = tostr(itime(3),2)+tostr(itime(4),2)+' UT'

   contour_circle, value, lons+utrot, lats, $
                   mini = mini, maxi = maxi, $
                   nLevels = nLevels, $
                    no00 = no00, no06 = no06, no12 = no12, no18 = no18, $
                   pos = pos, /nolines, $
                   maxrange = maxrange, title = title, $
                   /nomaxmin
;                    colorbar = colorbar, $
;                    latExtra = latExtra, lonExtra=lonExtra, $
;                    nocolor = nocolor, $
;                    nomaxmin = nomaxmin, thick = thick, average = average

   
   if ((iFile+1) mod ppp eq 0) then begin
      ctpos = pos
      ctpos(0) = ctpos(2)+0.01
      ctpos(2) = ctpos(0)+0.02
      plotct,254,ctpos,[mini,maxi],'Pedersen Conductance (mho)',/right
   endif

endfor



closedevice

end
