Pro Supermag_data_cvrt,dir,date,stnm

  ;convert SuperMag data file to AMIE format
  ;Written by Shasha Zou April 12, 2012

  cd,dir
  filename=date+'.0000.supermag.txt'
  openr,unit,filename,/get_lun
  openw,unit1,'mag'+date+'.final',/get_lun
  temp=''
  readf,unit,temp
  print, temp
  while (strpos(temp,'================================') eq -1) Do begin
     readf,unit,temp
     print,temp
  endwhile

;used to read daily data file from SuperMag with 1 min time resolution.
  For i=0,1439 Do begin
     readf,unit,temp
     print, 'should be time : ',temp
     result=strsplit(temp,string(9b),/extract)
     year=fix(result[0]) & month=fix(result[1]) & day=fix(result[2])
     hour=fix(result[3]) & mins=fix(result[4]) & sec=fix(result[5])
     stnm=fix(result[6])
     doy=julday(month,day,year)-julday(1,1,year)+1
     IAGA=make_array(stnm,/string)
     H=make_array(stnm,/long) & D=make_array(stnm,/long) & Z=make_array(stnm,/long)

     ; This is for the IMF line that SuperMAG people added::::
;     readf,unit,temp

     For j=0,stnm-1 Do begin
        readf,unit,temp
        result1=strsplit(temp,string(9b),/extract) ;data seperated by tap
        IAGA[j]=result1[0] & H[j]=double(result1[1]) & D[j]=double(result1[2]) & Z[j]=double(result1[3]) 
        if(H[j] eq 999999) then H[j]=-99999
        if(D[j] eq 999999) then D[j]=-99999
        if(Z[j] eq 999999) then Z[j]=-99999
        printf,unit1,year-2000,month,day,hour,mins,IAGA[j],'XYZ',H[j],D[j],Z[j],0,0,0,format='(3I2,1x,2I2,1x,a3,1x,a3,6I6)' ;format='(I2,1x,I3,1x,2I2,1x,a3,1x,a3,6I6)'
     Endfor
  Endfor

  close,unit
  free_lun,unit

  close,unit1
  free_lun,unit1

End

Pro Supermag_station_master,date,stnm
;generate masteryyyymmdd.dat file for AMIE run
;Written by Shasha Zou April 12, 2012

;read magnetometer information from Supermag_all_station.dat file
;magnetometer locations in geographic and AACGM coordinates 
;  dir = '/raid3/columbanus/Amie/Data/srcData/perm/'
  dir = '/raid1/Amie/src.v2.2/amie_files/perm/'
  file = 'Supermag_all_station.dat'
  print, 'Reading file : ',dir+file
  nlines=file_lines(dir+file)
  openr,unit,dir+file,/get_lun
  IAGA=make_array(nlines,/string) & Geog_lon=make_array(nlines,/float) & Geog_lat=make_array(nlines,/float)
  Geom_lon=make_array(nlines,/float) & Geom_lat=make_array(nlines,/float) & IAGA2=make_array(stnm,/string)
  st=''
  line = ''
  For i=0,nlines-1 Do begin

;IAGA     GLON     GLAT     MLON     MLAT   STATION-NAME
;SPA       0.00   -90.00    19.01   -74.08  South Pole Station

     readf,unit, line
     tmp = strsplit(line,' ',/extract)
;     readf,unit,st,glat,glon,mlat,mlon,da,oa,wt,format='(4X,A3,13X,2F7.2,F8.2,F7.2,2X,3F7.2)'
     IAGA[i]=tmp(0) & Geog_lat[i]=tmp(2) & Geog_lon[i]=tmp(1)
     Geom_lat[i]=tmp(4) & Geom_lon[i]=tmp(3)

  Endfor
  close,unit
  free_lun,unit

  openr,unit,'mag'+date+'.final',/get_lun
  For i=0,stnm-1 Do begin
     readf,unit,yr,mo,da,hr,mins,st,H,D,Z,her,der,zer,format='(3I2,1x,2I2,1x,a3,4x,6I6)'
     IAGA2[i]=st
  Endfor
  close,unit
  free_lun,unit

  openw,unit,'master'+date+'.txt',/get_lun
  printf,unit,'Header Line #1. The New AMIE ignores this line.'
  printf,unit,'Header Line #1. The New AMIE ignores this line.'
  For i=0,stnm-1 Do begin
     sub=where(IAGA2[i] eq IAGA)
     printf,unit,IAGA[sub],Geog_lat[sub],Geog_lon[sub],Geom_lat[sub],Geom_lon[sub],0,0,1,format='(4X,A3,13X,2F7.2,F8.2,F7.2,2X,3F7.2)'
  Endfor
  close,unit
  free_lun,unit

End

;___________main____________
;Written by Shasha Zou April 12, 2012

filelist = findfile('*0000.supermag.txt')

display, filelist
file = ask('file to process (full filename!)',filelist(0))

date=strmid(file,0,8)

dir='./'
Supermag_data_cvrt,dir,date,stnm
print,'Number of magnetometers for use: ',stnm
Supermag_station_master,date,stnm

end
