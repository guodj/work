!***************************************************************************************************
!
!GAMDM
!GLOBAL AVERAGE MASS DENSITY MODEL VERSION 1.0
!
!JOHN EMMERT 10/26/09
!
!EVALUATES THE GLOBAL AVERAGE MASS DENSITY MODEL DESCRIBED IN
!  Emmert, J. T., and J. M. Picone, (2010), Climatology of globally averaged 
!    thermospheric mass density, J. Geophys Res., doi:10.1029/2010JA015298.
!
!CALLING SEQUENCE:  CALL GAMDM(f107, f107a, day, kp, dens)
!
!INPUT ARGUMENTS
!     f107:      The daily F10.7 radio flux for the previous day, at the
!                Earth-Sun distance and in units of 10^-22 W/m^2/Hz.
!                When comparing with density derived from Two-Line Element sets
!                as described in Picone et al. [JGR, 2005, doi:10.1029/2004JA010585],
!                use a lag of three days instead of one day. 
!     f107a:     The 81-day average of f107, centered on the input day.
!     day:       Day of year (1-366; Day 1 = January 1). No distinction is made
!                between normal and leap years.
!     kp:        The daily average Kp index, 12 hours prior to the input time.
!
!                F10.7 and Kp indices are available from 
!                ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/
!
!OUTPUT ARGUMENT
!     dens:      A 3-element array containing the global average log-density at 
!                250, 400, and 550 km. The units are ln(kg/m^3). Model results at  
!                other altitudes are not currently available.
!
!***************************************************************************************************

subroutine gamdm(f107, f107a, day, kp, dens)

  implicit none
  
  real(4), intent(in)  :: f107, f107a, day, kp
  real(4), intent(out) :: dens(1:3)
  
  real(4)              :: terms(1:25), M(1:7), df107, df107a, dayrad
  real(4), parameter   :: daytorad = 2.0 * 3.14159265 / 366.0
  real(4), parameter   :: coeff(1:75) = &
    !250 km
    (/-2.48825E+01,-2.41365E+01,-2.38984E+01,-2.35296E+01,-2.32255E+01,-2.28893E+01,-2.31351E+01, &
      -3.31283E-03,-3.23169E-03, 7.13764E-05, &
       6.40899E-02, 2.58398E-02,-9.47404E-02,-5.11076E-02, &
      -8.78496E-03,-7.21274E-03, 9.48979E-03, 1.23627E-03, &
       8.03911E-05, 2.80653E-04,-8.08317E-05,-1.08498E-05, &
       3.43837E-02,-7.51443E-04,-1.92187E-04, &
    !400 km
      -2.93258E+01,-2.79776E+01,-2.74964E+01,-2.67615E+01,-2.60838E+01,-2.53981E+01,-2.57489E+01, &
      -5.85699E-03,-6.19354E-03, 1.38039E-04, &
       1.11026E-01, 4.61075E-02,-1.55642E-01,-8.14988E-02, &
      -1.69781E-02,-9.67317E-03, 1.51681E-02,-2.70232E-03, &
       2.20494E-04, 4.83850E-04,-3.42466E-04, 7.12279E-06, &
       5.42633E-02,-1.16037E-03,-2.56061E-04, &
     !550 km
      -3.24696E+01,-3.09152E+01,-3.03305E+01,-2.93446E+01,-2.83397E+01,-2.73605E+01,-2.75448E+01, &
      -7.02575E-03,-8.60290E-03, 1.78200E-04, &
       1.48822E-01, 6.31580E-02,-1.91211E-01,-1.08137E-01, &
      -2.13547E-02,-9.35550E-03, 1.79448E-02,-2.93148E-03, &
       4.16032E-04, 6.15307E-04,-7.19811E-04, 7.09470E-05, &
       7.24543E-02,-1.07551E-03,-3.18109E-04/)

  call f107spl3(f107, M)
  terms(1:7) = M
  
  df107 = f107 - f107a
  df107a = f107a - 150.0
  if (df107 .lt. 0) then
    terms(8) = df107
    terms(9:10) = 0.0
  else
    terms(8) = 0.0
    terms(9) = df107
    terms(10) = df107 * df107a
  end if
  
  dayrad = day * daytorad
  terms(11) = cos(dayrad)
  terms(12) = sin(dayrad)
  terms(13) = cos(2.0*dayrad)
  terms(14) = sin(2.0*dayrad)
  terms(15) = cos(3.0*dayrad)
  terms(16) = sin(3.0*dayrad)
  terms(17) = cos(4.0*dayrad)
  terms(18) = sin(4.0*dayrad)
  terms(19:22) = terms(11:14)*df107a
  
  terms(23) = kp - 1.6
  if (df107a .lt. 0) then
    terms(24) = terms(23) * df107a
    terms(25) = 0.0
  else
    terms(24) = 0.0
    terms(25) = terms(23) * df107a
  end if

  dens(1) = dot_product(coeff(1:25), terms) 
  dens(2) = dot_product(coeff(26:50), terms) 
  dens(3) = dot_product(coeff(51:75), terms) 
  
  return
  
end subroutine gamdm

!***************************************************************************************************

subroutine f107spl3(f107, M)

  implicit none
  
  real(4), intent(in)  :: f107
  real(4), intent(out) :: M(1:7)

  integer(4)           :: i, j
  real(4)              :: N(1:12)
  real(4), parameter   :: node(1:13) = (/30,40,50,60,70,80,100,150,220,320,420,520,620/)
  
  do i = 1, 12
    N(i) = 0.0
    if ((f107 .ge. node(i)) .and. (f107 .lt. node(i+1))) N(i) = 1.0
  end do
  do j = 2, 4
    do i = 1, 13-j
      N(i) =   N(i)   * (f107 - node(i))   / (node(i+j-1) - node(i)) &
             + N(i+1) * (node(i+j) - f107) / (node(i+j)   - node(i+1))
    end do
  end do

  M(1) = N(1) + 0.5*N(2)
  M(2) = 0.5*N(2) + N(3)
  M(3:5) = N(4:6)
  M(6) = N(7) + 0.526316*N(8)
  M(7) = 0.473684*N(8) + N(9)
  
  return
  
end subroutine f107spl3

!***************************************************************************************************

!TEST DRIVER

program gamdm_test

  implicit none
  real(4)      :: f107, f107a, day, kp, dens(1:3), df107
  integer(4)   :: if107, if107a, iday, ikp, idf107

!F10.7 PROFILE
  day = 305.0
  kp = 3.7
  print *, 'DAILY F10.7 PROFILE'
  print '(a4,f0.0, a5,f3.1)', 'DAY=',day, ', KP=',kp
  print '(2a8,3a12)', 'F10.7', 'F10.7a','250 km', '400 km', '550 km'
  do if107 = 70, 220, 10
    f107 = float(if107)
    f107a = f107
    call GAMDM(f107,f107a,day,kp,dens)
    print '(2f8.0,3f12.4)', f107, f107a, dens
  end do
  print *
  print *

!DELTA-F10.7 PROFILES
  day = 90.0
  kp = 1.4
  print *, 'DELTA-F10.7 PROFILES'
  print '(a4,f0.0, a5,f3.1)', 'DAY=',day, ', KP=',kp
  print '(3a8,3a12)', 'F10.7', 'F10.7a','DF10.7','250 km', '400 km', '550 km'
  do if107a = 80, 200, 60
    do idf107 = -60, 60, 30
      f107a = float(if107a)
      df107 = float(idf107)
      f107 = df107 + f107a
      if (f107 .lt. 70) cycle
      call GAMDM(f107,f107a,day,kp,dens)
      print '(3f8.0,3f12.4)', f107, f107a, df107, dens
    end do
  end do
  print *
  print *

!DAY-OF-YEAR PROFILE
  f107 = 170.0
  f107a = 170.0
  kp = 5.4
  print *, 'DAY-OF-YEAR PROFILE'
  print '(a6,f0.0, a9,f0.0, a5,f3.1)', 'F10.7=',f107, ', F10.7a=',f107a, ', KP=',kp
  print '(a4,3a12)', 'DAY','250 km', '400 km', '550 km'
  do iday = 0, 360, 20
    day = float(iday)
    call GAMDM(f107,f107a,day,kp,dens)
    print '(f4.0,3f12.4)',day, dens
  end do
  print *
  print *

!KP PROFILE
  f107 = 100.0
  f107a = 100.0
  day = 200.0
  print *, 'KP PROFILE'
  print '(a6,f0.0, a9,f0.0, a6,f0.0)', 'F10.7=',f107, ', F10.7a=',f107a, ', DAY=',day
  print '(a4,3a12)', 'KP','250 km', '400 km', '550 km'
  do ikp = 0, 7, 1
    kp = float(ikp)
    call GAMDM(f107,f107a,day,kp,dens)
    print '(f4.1,3f12.4)', kp, dens
  end do
  print *
  print *

end program gamdm_test

!***************************************************************************************************

!TEST DRIVER OUTPUT

! DAILY F10.7 PROFILE
!DAY=305., KP=3.7
!   F10.7  F10.7a      250 km      400 km      550 km
!     70.     70.    -23.8518    -27.5384    -30.4179
!     80.     80.    -23.6626    -27.1606    -29.9534
!     90.     90.    -23.5277    -26.8807    -29.5835
!    100.    100.    -23.4167    -26.6479    -29.2585
!    110.    110.    -23.3264    -26.4530    -28.9768
!    120.    120.    -23.2536    -26.2897    -28.7331
!    130.    130.    -23.1945    -26.1521    -28.5213
!    140.    140.    -23.1457    -26.0347    -28.3353
!    150.    150.    -23.1037    -25.9315    -28.1691
!    160.    160.    -23.0539    -25.8191    -28.0017
!    170.    170.    -23.0078    -25.7159    -27.8484
!    180.    180.    -22.9660    -25.6227    -27.7096
!    190.    190.    -22.9289    -25.5400    -27.5859
!    200.    200.    -22.8972    -25.4685    -27.4776
!    210.    210.    -22.8714    -25.4090    -27.3852
!    220.    220.    -22.8520    -25.3619    -27.3093
! 
! 
! DELTA-F10.7 PROFILES
!DAY=90., KP=1.4
!   F10.7  F10.7a  DF10.7      250 km      400 km      550 km
!     80.     80.      0.    -23.8685    -27.4850    -30.3268
!    110.     80.     30.    -23.7273    -27.1757    -29.9139
!    140.     80.     60.    -23.7415    -27.1557    -29.8362
!     80.    140.    -60.    -23.6390    -27.0697    -29.8114
!    110.    140.    -30.    -23.3503    -26.4605    -28.9770
!    140.    140.      0.    -23.2170    -26.1404    -28.4778
!    170.    140.     30.    -23.1690    -26.0095    -28.2656
!    200.    140.     60.    -23.1601    -25.9689    -28.1856
!    140.    200.    -60.    -22.9931    -25.7342    -27.9700
!    170.    200.    -30.    -22.9261    -25.5518    -27.6571
!    200.    200.      0.    -22.8982    -25.4597    -27.4762
!    230.    200.     30.    -22.8134    -25.2771    -27.2184
!    260.    200.     60.    -22.7829    -25.1986    -27.0942
! 
! 
! DAY-OF-YEAR PROFILE
!F10.7=170., F10.7a=170., KP=5.4
! DAY      250 km      400 km      550 km
!  0.    -23.0892    -25.8494    -28.0062
! 20.    -23.1015    -25.8670    -28.0260
! 40.    -23.0682    -25.8032    -27.9458
! 60.    -22.9964    -25.6785    -27.7930
! 80.    -22.9296    -25.5744    -27.6663
!100.    -22.9192    -25.5704    -27.6621
!120.    -22.9764    -25.6715    -27.7880
!140.    -23.0652    -25.8150    -27.9703
!160.    -23.1438    -25.9423    -28.1371
!180.    -23.1972    -26.0355    -28.2631
!200.    -23.2239    -26.0858    -28.3317
!220.    -23.2115    -26.0628    -28.3030
!240.    -23.1480    -25.9490    -28.1558
!260.    -23.0529    -25.7873    -27.9456
!280.    -22.9753    -25.6618    -27.7800
!300.    -22.9525    -25.6270    -27.7302
!320.    -22.9813    -25.6724    -27.7849
!340.    -23.0322    -25.7542    -27.8878
!360.    -23.0789    -25.8321    -27.9851
! 
! 
! KP PROFILE
!F10.7=100., F10.7a=100., DAY=200.
!  KP      250 km      400 km      550 km
! 0.0    -23.9454    -27.4826    -30.2435
! 1.0    -23.8735    -27.3703    -30.1172
! 2.0    -23.8015    -27.2580    -29.9910
! 3.0    -23.7295    -27.1457    -29.8648
! 4.0    -23.6576    -27.0334    -29.7386
! 5.0    -23.5856    -26.9212    -29.6123
! 6.0    -23.5137    -26.8089    -29.4861
! 7.0    -23.4417    -26.6966    -29.3599

!***************************************************************************************************
