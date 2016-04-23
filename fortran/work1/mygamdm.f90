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

program mygamdm

  implicit none
  character *100 fname1,fname2,fname3
  real(4)      :: year,doy,f107, f107a,kp, ikp=2.2, dens(1:3),denskp(1:3)
  integer(4)   :: u1=1,u2=2,u3=3,i,j,iost
    fname1='/home/gdj/study/graduation_project/data/processed_data/gamdm_index.txt'
    fname2='/home/gdj/study/graduation_project/data/processed_data/gamdm_density.txt '
    fname3='/home/gdj/study/graduation_project/data/processed_data/gamdm_density_kp2.txt'
    open(unit=u1,file=fname1,iostat=iost)
    open(unit=u2,file=fname2)
    open(unit=u3,file=fname3)
    write(2,*) 'year  doy  rho250 rho400 rho550'
    write(3,*) 'year  doy  rho250 rho400 rho550,kp=2.2'
    do while (iost==0)
        read (1,*)year,doy,f107, f107a,kp
        call GAMDM(f107,f107a,doy,kp,dens)
        call GAMDM(f107,f107a,doy,ikp,denskp)
        dens=exp(dens)
        denskp=exp(denskp)
        write(2,100)int(year),int(doy),dens
        write(3,100)int(year),int(doy),denskp
    end do
100 format(i7,i7,3E12.4)
    close(u1)
    close(u2)
    close(u3)
end
