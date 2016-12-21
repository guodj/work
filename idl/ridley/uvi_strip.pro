;$Id: norm.pro,v 1.6 1994/11/29 18:44:16 beth Exp $
;
; Copyright (c) 1994, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.
;+
; NAME:
;       NORM
;
; PURPOSE:
;       1) This function computes the Euclidean norm of a vector.
;          OR
;       2) This function computes the Infinity norm of an array.
;
; CATEGORY:
;       Complex Linear Algebra.
;
; CALLING SEQUENCE:
;       Result = NORM(A)
;
; INPUTS:
;       A:      An N-element real or complex vector.
;               An M by N real or complex array.
;
; KEYWORD PARAMETERS:
;       DOUBLE: If set to a non-zero value, computations are done in
;               double precision arithmetic.
;
; EXAMPLE:
;       1) Define an N-element complex vector (a).
;            a = [complex(1, 0), complex(2,-2), complex(-3,1)]
;          Compute the Euclidean norm of (a).
;            result = norm(a)
;
;       2) Define an M by N complex array (a). 
;            a = [[complex(1, 0), complex(2,-2), complex(-3,1)], $
;                 [complex(1,-2), complex(2, 2), complex(1, 0)]]
;          Compute the Infinity norm of the complex array (a).
;            result = norm(a)
;
; PROCEDURE:
;       NORM.PRO computes the Euclidean norm of an N-element vector.
;       NORM.PRO computes the Infinity norm of an M by N array
;
; REFERENCE:
;       ADVANCED ENGINEERING MATHEMATICS (seventh edition)
;       Erwin Kreyszig
;       ISBN 0-471-55380-8
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, April 1992
;       Modified:    GGS, RSI, February 1994
;                    Computes the Euclidean norm of an N-element vector.
;                    Accepts complex inputs. Added DOUBLE keyword.
;       Modified:    GGS, RSI, September 1994
;                    Added support for double-precision complex inputs.
;-

function NORM, a, double = double
  
  on_error, 2  ;Return to caller if error occurs.

  type = size(a) 

  if type(0) eq 1 then begin ;If vector, compute the Euclidean norm.
    if keyword_set(DOUBLE) ne 0 then begin 
      if type(2) ne 6 and type(2) ne 9 then a = double(a) $
      else a = dcomplex(a) 
    endif
    return, sqrt(total(abs(a)^2)) 
  endif else if type(0) eq 2 then begin ;If matrix, compute the Infinity norm.
    if keyword_set(DOUBLE) ne 0  then begin
      if type(3) ne 6 and type(3) ne 9 then a = double(a) $
      else a = dcomplex(a)
    endif
    return, max(total(abs(a),1))   ;Create the matrix of absolute values,
                                   ;add the elements of each row,
                                   ;and return the maximum.
  endif else message, $
    'Input must be an n-element vector or an M by N array.'

end
;
; IC_GCI_TRANSF
;
; Routines to transform from GCI coordinates to GEO coordinates.
; The following ISTP ICSS routines are appended in this file.
; 
; IC_CONV_MATRIX.FOR
; IC_GCI_TO_GEO.FOR
; IC_GET_NUT_ANGLES.FOR
; IC_GRNWCH_SIDEREAL.FOR
;
; In addition, the NAG routine F01CKF is emulated by a subroutine
; included in this file.  Consequently, this code should not be
; linked with the NAG libraries if available.  The IC_ routines
; have not been modified in any fashion.
;
; A good reference for the routines below is 'An Explanatory
; Supplement to the Astronomical Almanac,' P. Kenneth Seidelmann, ed.
; University Science Books, 1992.  ISBN 0-935702-68-7
;
; G. GERMANY   8/9/95
;
;------------------------------------------------------------------


;
; IC_GCI_TO_GEO - return a transformation matrix
;
; PURPOSE:  Calculate the transformation matrix from GCI 
;           coordinates to GEO coordinates at a given date and time.
;
; UNIT TYPE:  SUBROUTINE
;
; INVOCATION METHOD:  CALL IC_GCI_TO_GEO (orb_pos_time, 
;                                         transform_matrix)
;
; argUMENT LIST:
;
; NAME	                  TYPE   USE  DESCRIPTION
; ----                   ----   ---  -----------
; orb_pos_time(2)        I*4    I    time OF ORB. VECTOR, year-day-MILLI OF day
; transform_matrix(3,3)  R*8    O    TRANSFORMATION MATRIX
;
; FILE/RECORD REFERENCES:  NONE
;
; EXTERNAL VARIABLES:  NONE
;
; EXTERNAL REFERENCES:
;	IC_GRNWCH_SIDEREAL - Returns the Greenwich sidereal time in radians
;	F01CKF - Multiplies two matrices
;	IC_CONV_MATRIX - Returns the conversion matrix to rotate
;                       from mean of date to true of date
;
; ABNORMAL TERMINATION CONDITIONS, ERROR MESSAGES:  NONE
;
; ASSUMPTIONS, CONstrAINTS, REstrICTIONS:  NONE
;
; DEVELOPMENT HISTORY
;
; AUTHOR	CHANGE ID	RELEASE	  DATE	    DESCRIPTION OF CHANGE
; ------	---------	-------   ----	    ---------------------
; J. LUBELCZYK                 B1R1      11/21/90  INITIAL PDL
; J. Lubelczyk			B1R1      12/10/90  CODING
; J. Lubelczyk                           09/18/91  Updated to return true
;                                                  of date trans matrix
; J. LUBELCZYK ICCR #83, CCR #'S 130, 137 11/91    B3 update
;
; NOTES:
;

PRO ic_gci_to_geo, orb_pos_time, transform_matrix


        mean_matrix = DBLARR(3,3)
	transform_matrix = DBLARR(3,3) ;transformation matrix
	cmatrix = DBLARR(3,3)	      ;the conversion matrix to rotate from 
                                      ;mean of date to true of date

	ic_grnwch_sidereal, orb_pos_time, grnwch_sidereal_time

;
;   calculate the sin and cos of the greenwich mean sidereal time
;
	sin_gst = sin(grnwch_sidereal_time)
	cos_gst = cos(grnwch_sidereal_time)

;
;   Fill the mean of date transformation matrix using the sin and cos
;    of the greenwich mean sidereal time
;
	mean_matrix(0,0) = cos_gst
	mean_matrix(0,1) = -sin_gst
	mean_matrix(0,2) = 0
	mean_matrix(1,0) = sin_gst
	mean_matrix(1,1) = cos_gst
	mean_matrix(1,2) = 0
	mean_matrix(2,0) = 0
	mean_matrix(2,1) = 0
	mean_matrix(2,2) = 1

	ic_conv_matrix,orb_pos_time, cmatrix

        transform_matrix = TRANSPOSE( $
                           TRANSPOSE(mean_matrix) # TRANSPOSE(cmatrix) )

END


;      STATEMENT FUNCTION DEFINITION FOR dxjul -- JULIAN EPHEMERIS
;      DATE AT BEGINNING OF year.
;

FUNCTION dxjul,i

	RETURN,DOUBLE((-32075+1461*(i+4800-13/12)/4  $
                     +367*(-1+13/12*12)/12-3         $
                     *((i+4900-13/12)/100)/4)-0.5	)

END


; IC_CONV_MATRIX - Returns the conversion matrix that is necessary to rotate
;		    from mean of date to true of date
;
; PURPOSE:  THIS SUBROUTINE CALCULATES, THROUGH APPROPRIATE ANALYTIC
;           EXPRESSIONS, VALUES FOR THE PRECESSION AND NUTATION
;           ANGLES AND THE MATRIX REQUIRED TO ROTATE FROM MEAN OF
;           JULIAN 2000 TO TRUE OF DATE.
; NAME	               TYPE   USE  DESCRIPTION
; ----                ----   ---  -----------
; orb_pos_time(2)     I*4    I    time OF ORB. VECTOR, year-day-MILLI OF day
; CMATRIX(3,3)        R*8    O    Matrix to rotate from J2000 to true of date
;
; EXTERNAL REFERENCES:
;	F01CKF - NAG routine that multiplies two matrices
;      IC_GET_NUT_ANGLES          Routine to compute the nutation angles 
;                            

PRO ic_conv_matrix,orb_pos_time, cmatrix

	cmatrix = DBLARR(3,3)    ;conversion matrix
	nutmat = DBLARR(3,3)	    ;Nutation matrix
	premat = DBLARR(3,3)	    ;precession matrix

	Jdj2000 = DOUBLE (2451545.0)
        R = DOUBLE (1296000.0)


;   Convert the given millisecond of day [orb_pos_time(1)] to second of day.
;   Convert the packed form into year and day-of-year

      secs = (DOUBLE(orb_pos_time(1)))/DOUBLE(1000.0)
      year = orb_pos_time(0)/1000
      day  = orb_pos_time(0) MOD 1000

;
;   Calculate the julian date and the time in Julian centuries from J2000
;
      fday = secs/DOUBLE(86400.00)
      jul_day = dxjul(year) + DOUBLE(day)+fday
      time = (jul_day - Jdj2000)/DOUBLE(36525.0)

      T2 = time*time
      T3 = time*T2

;    CALCULATE CONVERSION FACTORS: DEGREES TO RADANS (dtr), SECONDS TO 
;    RADIANS (str)                                                     
;                                                                      
      PI= DOUBLE(4.0 * ATAN(1.0))
      dtr=PI/DOUBLE(180.0)
      str=dtr/DOUBLE(3600.0)


;    CALCULATE PRECESSION ANGLES

      zeta   = DOUBLE(                     	$
               0.11180860865024398D-01*time	$
             + 0.14635555405334670D-05*T2	$
             + 0.87256766326094274D-07*T3 )
      theta  = DOUBLE(                     	$
               0.97171734551696701D-02*time	$
             - 0.20684575704538352D-05*T2	$
             - 0.20281210721855218D-06*T3 )
      zee    = DOUBLE(                     	$
               0.11180860865024398D-01*time	$
             + 0.53071584043698687D-05*T2	$
             + 0.88250634372368822D-07*T3 )

      sinzet = sin(zeta)                                               
      coszet = cos(zeta)                                               
      sinzee = sin(zee)   
      coszee = cos(zee)   
      sinthe = sin(theta) 
      costhe = cos(theta) 
;                                                                      
;    COMPUTE THE TRANSFORMATION MATRIX BETWEEN MEAN EQUATOR AND        
;    EQUINOX OF 1950 AND MEAN EQUATOR AND EQUINOX OF DATE. THIS        
;    MATRIX IS CALLED premat.                                          
;                                                                      
      premat(0,0) = -sinzet*sinzee  + coszet*coszee*costhe              
      premat(0,1) =  coszee*sinzet  + sinzee*costhe*coszet              
      premat(0,2) =  sinthe*coszet                                      
      premat(1,0) = -sinzee*coszet  - coszee*costhe*sinzet              
      premat(1,1) =  coszee*coszet  - sinzee*costhe*sinzet              
      premat(1,2) = -sinthe*sinzet                                      
      premat(2,0) = -coszee*sinthe                                      
      premat(2,1) = -sinzee*sinthe                                      
      premat(2,2) =  costhe                                             

;    CALCULATE MEAN OBLIQUITY OF DATE (epso). WHERE TIME IS MEASURED IN
;    JULIAN CENTURIES FROM 2000.0.

      epso=DOUBLE( (1.813E-3*T3-5.9E-4*T2	$
                    -4.6815E+1*time+8.4381448E+4)*str )


;    CALL IC_GET_NUT_ANGLES TO COMPUTE NUTATION IN OBLIQUITY AND LONGITUDE

      ic_get_nut_angles,time,deleps,delpsi,eps

      cosep=cos(eps)                               
      cosepO=cos(epso)                             
      cospsi=cos(delpsi)                           
      sinep=sin(eps)                               
      sinepO=sin(epso)                             
      sinpsi=sin(delpsi)                           

      nutmat(0,0)=cospsi                            
      nutmat(1,0)=-sinpsi*cosepO                    
      nutmat(2,0)=-sinpsi*sinepO                    
      nutmat(0,1)=sinpsi*cosep                      
      nutmat(1,1)=cospsi*cosep*cosepO+sinep*sinepO  
      nutmat(2,1)=cospsi*cosep*sinepO-sinep*cosepO  
      nutmat(0,2)=sinpsi*sinep                      
      nutmat(1,2)=cospsi*sinep*cosepO-cosep*sinepO  
      nutmat(2,2)=cospsi*sinep*sinepO+cosep*cosepO  

;    CALCULATE ELEMENTS OF nutmat * premat.  THIS MATRIX IS THE 
;    ANALYTICALLY CALCULATED TRANSFORMATION MATRIX, WHICH WILL  
;    TRANSFORM THE MEAN EARTH EQUATOR AND EQUINOX OF J2000 INTO
;    THE TRUE EARTH EQUATOR AND EXQUINOX OF DATE.

;     cmatrix = nutmat # premat
      cmatrix = TRANSPOSE(TRANSPOSE(premat) # TRANSPOSE(nutmat))

END




; IC_GET_NUT_ANGLES - Returns angles that are necessary to adjust the
;		       Greenwich Hour angle to true of date
;
; PURPOSE : THIS SUBROUTINE CALCULATES, THROUGH APPROPRIATE ANALYTIC
;           EXPRESSIONS, VALUES FOR THE NUTATION
;           ANGLES TO ROTATE FROM MEAN OF
;           JULIAN 2000 TO TRUE OF DATE.
;
; NAME	      TYPE    USE	DESCRIPTION
; ----       ----    ---       -----------
; time       R*8     I         time IN JULIAN CENTURIES OF 36525.0
;                              MEAN SOLAR dayS FROM J2000. (NOTE: THIS
;                              CAN BE POSITIVE OR NEGATIVE.)
; deleps     R*8     O	        DELTA epsILON, Nutation in obliquity
; delpsi     R*8     O         DELTA PSI, Nutation in longitude
; eps        R*8     O         epsILON



PRO ic_get_nut_angles,time,deleps,delpsi,eps


	TOL = DOUBLE(0.5)/DOUBLE(36525.0)


        isinco = FLTARR(2,106)   ;Array used in nutation calculations
	icosco = FLTARR(2,106)   ;Array used in nutation calculations
	fund = DBLARR(1,5)	    ;The fundamental arguments
	T = DBLARR(1,2)	    ;Time
	ifunar = intarr(5,106)   ;array used in nutation calculations
        oldtim = DOUBLE(0.0)



;                                                                     
;    INITIALIZE VALUES OF IFUNAR, isinco, AND icosco ARRAYS FOR USE IN
;    NUTATION CALCULATIONS.                                           
;                                                                     
;    THE 1980 IAU THEORY OF NUTATION,CONTAINED IN JPL
;    DE200 PLANETARY EPHEMERIS.


      ifunar(0,0:79) = [0,0,2,-2,2, 0,2,-2,2,0, 2,2,2,0,2, 0,0,2,0,0, $
                2,0,2,0,0, -2,-2,0,0,2, 2,0,2,2,0, 2,0,0,0,2, $
                2,2,0,2,2,  2,2,0,0,2,  0,2,2,2,0, 2,0,2,2,0, $
                0,2,0,-2,0, 0,2,2,2,0,  2,2,2,2,0, 0,0,2,0,0]
      ifunar(0,80:105) = [2,2,0,2,2,  2,4,0,2,2,  0,4,2,2,2, 0,-2,2,0,-2, $
                2,0,-2,0,2, 0]

      ifunar(1,0:79) = [1,2,1,0,2, 0,1,1,2,0,  2,2,1,0,0, 0,1,2,1,1, $
                1,1,1,0,0, 1,0,2,1,0,  2,0,1,2,0, 2,0,1,1,2, $
                1,2,0,2,2, 0,1,1,1,1,  0,2,2,2,0, 2,1,1,1,1, $
                0,1,0,0,0, 0,0,2,2,1,  2,2,2,1,1, 2,0,2,2,0]
      ifunar(1,80:105) = [2,2,0,2,1, 2,2,0,1,2,  1,2,2,0,1, 1,1,2,0,0, $
                1,1,0,0,2, 0]

      ifunar(2,0:79) = [0,0,0,0,0, -1,-2,0,0,1,  1,-1,0,0,0, 2,1,2,-1,0, $
               -1,0,1,0,1,  0,1,1,0,1,   0,0,0,0,0,  0,0,0,0,0, $
                0,0,0,0,0,  0,0,0,0,0,	 1,1,-1,0,0, 0,0,0,0,0, $
               -1,0,1,0,0,  1,0,-1,-1,0, 0,-1,1,0,0, 0,0,0,0,0]
      ifunar(2,80:105) = [0,0,0,1,0,  0,0,-1,0,0,  0,0,0,0,1, -1,0,0,1,0,	$
               -1,1,0,0,0,  1]

      ifunar(3,0:79) = [0,0,-2,2,-2, 1,0,2,0,0,    0,0,0,2,0,   0,0,0,0,-2, $
                0,2,0,1,2,   0,0,0,-1,0,   0,1,0,1,1,  -1,0,1,-1,-1, $
                1,0,2,1,2,   0,-1,-1,1,-1, 1,0,0,1,1,   2,0,0,1,0, $
                1,2,0,1,0,   1,1,1,-1,-2,  3,0,1,-1,2,  1,3,0,-1,1]
      ifunar(3,80:105) = [-2,-1,2,1,1, -2,-1,1,2,2,   1,0,3,1,0,  -1,0,0,0,1, $
                0,1,1,2,0,   0]

      ifunar(4,0:79) = [0,0,0,0,0,    -1,-2,0,-2,0, -2,-2,-2,-2,-2, 0,0,-2,0,2, $
               -2,-2,-2,-1,-2, 2,2,0,1,-2,   0,0,0,0,-2,    0,2,0,0,2, $
                0,2,0,-2,0,    0,0,2,-2,2,  -2,0,0,2,2,    -2,2,2,-2,-2, $
                0,0,-2,0,1,    0,0,0,2,0,    0,2,0,-2,0,    0,0,1,0,-4]
      ifunar(4,80:105) = [2,4,-4,-2,2,   4,0,-2,-2,2,  2,-2,-2,-2,0,  2,0,-1,2,-2, $
                0,-2,2,2,4,    1]


      isinco(0,0:79) = [-171996.,2062.,46.,11.,-3.,   -3.,-2.,1.,-13187.,1426.,	$
                -517.,217.,129.,48.,-22.,     17.,-15.,-16.,-12.,-6., $
                -5.,4.,4.,-4.,1.,              1.,-1.,1.,1.,-1., $
                -2274.,712.,-386.,-301.,-158., 123.,63.,63.,-58.,-59., $
                -51.,-38.,29.,29.,-31.,        26.,21.,16.,-13.,-10., $
                -7.,7.,-7.,-8.,6.,             6.,-6.,-7.,6.,-5., $
                 5.,-5.,-4.,4.,-4.,           -3.,3.,-3.,-3.,-2., $
                -3.,-3.,2.,-2.,2.,            -2.,2.,2.,1.,-1.]
       isinco(0,80:105) = [1.,-2.,-1.,1.,-1.,           -1.,1.,1.,1.,-1., $
                -1.,1.,1.,-1.,1.,              1.,-1.,-1.,-1.,-1., $
                -1.,-1.,-1.,1.,-1.,            1.]

       isinco(1,0:38) = [-174.2,.2,0.,0.,0.,            0.,0.,0.,-1.6,-3.4, $
                 1.2,-.5,.1,0.,0.,            -.1,0.,.1,0.,0., $
                 0.,0.,0.,0.,0.,               0.,0.,0.,0.,0., $
               -.2,.1,-.4,0.,0.,               0.,0.,.1,-.1]

      icosco(0,0:79) = [92025.,-895.,-24.,0.,1.,   0.,1.,0.,5736.,54., $
                224.,-95.,-70.,1.,0.,      0.,9.,7.,6.,3., $
                3.,-2.,-2.,0.,0.,          0.,0.,0.,0.,0., $
                977.,-7.,200.,129.,-1.,   -53.,-2.,-33.,32.,26., $
                27.,16.,-1.,-12.,13.,     -1.,-10.,-8.,7.,5., $
                0.,-3.,3.,3.,0.,          -3.,3.,3.,-3.,3., $
                0.,3.,0.,0.,0.,            0.,0.,1.,1.,1., $
                1.,1.,-1.,1.,-1.,          1.,0.,-1.,-1.,0.]
      icosco(0,80:105) = [-1.,1.,0.,-1.,1.,           1.,0.,0.,-1.,0., $
                0.,0.,0.,0.,0.,            0.,0.,0.,0.,0., $
                0.,0.,0.,0.,0.,            0.]

      icosco(1,0:33) = [8.9,.5,0.,0.,0.,           0.,0.,0.,-3.1,-.1, $
                -.6,.3,0.,0.,0.,           0.,0.,0.,0.,0., $
                 0.,0.,0.,0.,0.,           0.,0.,0.,0.,0.,	$
               -.5,0.,0.,-.1]

      R = DOUBLE(1296000.0)

      IF ( ABS ( time - oldtim ) LE TOL ) THEN BEGIN
         deleps = olddep 
         delpsi = olddps 
         eps    = oldeps
         RETURN
      ENDIF

      T2 = time*time
      T3 = time*T2

;    CONVERT IFUNAR, isinco, AND icosco ARRAYS TO REAL*8 ARRAYS FUNarg,
;    sincof, AND coscof, RESPECTIVELY.                                 
;                                                                      
      funarg = DOUBLE(ifunar)
      sincof = DOUBLE(isinco)
      coscof = DOUBLE(icosco)


;    CALCULATE CONVERSION FACTORS: DEGREES TO RADANS (dtr), SECONDS TO 
;    RADIANS (str)                                                     

      PI= DOUBLE(4.0 * ATAN(1.0))
      dtr=PI/DOUBLE(180.0)
      str=dtr/DOUBLE(3600.0)

;    BEGIN COMPUTATION OF NUTATION IN OBLIQUITY AND LONGITUDE          

;    CALCULATE FUNDAMENTAL argUMENTS FOR USE IN NUTATION CALCULATIONS
;    time IS REFERENCED TO J2000.0.
;    fund(1,1)= F
;    fund(2,1)= OMEGA
;    fund(3,1)= L PRIME
;    fund(4,1)= L
;    fund(5,1)= D
	
      fund(0,0)=str*(335778.877E0+(1342.0E0*R+295263.137E0)*time	$
      -13.257E0*T2+1.1E-2*T3)
      fund(0,1)=str*(450160.280E0-(5.E0*R+482890.539E0)*time+	$
      7.455E0*T2+8.0E-3*T3)
      fund(0,2)=str*(1287099.804E0+(99.0E0*R+1292581.224E0)*time-	$
      5.77E-1*T2-1.2E-2*T3)
      fund(0,3)=str*(485866.733E0+(1325.0E0*R+715922.633E0)*time+	$
      31.310E0*T2+6.4E-2*T3)
      fund(0,4)=str*(1072261.307E0+(1236.0E0*R+1105601.328E0)*time-	$
      6.891E0*T2+1.9E-2*T3)


;    CALCULATE MEAN OBLIQUITY OF DATE (epso). WHERE time IS MEASURED IN
;    JULIAN CENTURIES FROM 2000.0.

      epso=(1.813E-3*T3-5.9E-4*T2-4.6815E+1*time+8.4381448E+4)*str


;                                                                      
;    CALCULATE NUTATION IN LONGITUDE (delpsi) AND NUTATION IN OBLIQUITY
;    (deleps).  THIS IS A THREE STEP PROCESS:                          
;    (1) CALCULATE argUMENTS OF sinE (FOR delpsi) AND COsinE (FOR deleps)
;        THESE ARE OF THE FORM                                        
;                                                                     
;        arg = SUMMATION ( A(I) * fund(I,1) ), I = 1,5                 
;                                                                      
;        WHERE THE A(I)'S ARE ELEMENTS OF FUNarg.                      
;                                                          
;      arg = funarg # fund
;      arg = fund # funarg
      arg = TRANSPOSE(TRANSPOSE(funarg) # TRANSPOSE(fund))

;                                                                      
;    (2) CALCULATE COEFFICIENTS OF sinE AND COsinE, WHICH ARE THE PRODUCTS
;        OF sincof * T AND coscof * T.  THESE COEFFICIENTS ARE IN UNITS
;        OF 0.0001 SECONDS OF ARC.                                     
;                                                                      
      T(0,0)=DOUBLE(1.0)
      T(0,1)=time       
                         
;      cofcos = coscof # T
;      cofsin = sincof # T                       
;      cofcos = T # coscof
;      cofsin = T # sincof
      cofcos = TRANSPOSE(TRANSPOSE(coscof) # TRANSPOSE(T))
      cofsin = TRANSPOSE(TRANSPOSE(sincof) # TRANSPOSE(T))

      cofcos=cofcos*DOUBLE(1.E-4)
      cofsin=cofsin*DOUBLE(1.E-4)

;                                                                      
;    (3) CALCULATE THE sinES AND COsinES OF THE argUMENTS AND MULTIPLY 
;        BY THEIR COEFFICIENTS, THEN ADD.  COMPUTE delpsi AND deleps.  
;                                                                      
      sumpsi=DOUBLE(0.0)
      sumeps=DOUBLE(0.0)

      sinp=sin(arg)
      cose=cos(arg)

      FOR E=0,105 DO BEGIN
         prodps=cofsin(0,E)*sinp(0,E)                                    
         prodep=cofcos(0,E)*cose(0,E)                                    
         sumpsi=sumpsi+prodps                                              
         sumeps=sumeps+prodep  
      ENDFOR

      deleps=sumeps*str                                                 
      delpsi=sumpsi*str                  

;                                        
;    CALCULATE TRUE OBLIQUITY OF DATE (eps).       
;                                                  
      eps=epso+deleps
      olddep = deleps
      olddps = delpsi
      oldeps = eps
      oldtim = time

END



; IC_GRNWCH_SIDEREAL - return the greenwich true sidereal time in radians
;
; PURPOSE:  Calculate the true of date greenwich sidereal time in radians.
;
; NAME	                  TYPE   USE  DESCRIPTION
; ----                   ----   ---  -----------
; orb_pos_time(2)        I*4    I    time OF ORB. VECTOR, year-day-MILLI OF day
; gst		          R*8    O    GREENWICH MEAN SIDEREAL time
; EXTERNAL REFERENCES:
;	IC_GET_NUT_ANGLES - Returns angles necessary to adjust the Greenwich
;                          hour angle to true of date
; NOTES:
; 1)  THE ORIGINAL ALGORITHM USED WAS COPIED FROM A SHORT PROGRAM BY 
;     G. D. MEAD, INCLUDED IN 'GEOPHYSICAL COORDINATE 
;                      TRANSFORMATIONS' BY CHRISTOPHER T. RUSSELL
; 2)  THIS VERSION INCORPORATES SEVERAL CHANGES TO CALCULATE THE GREENWICH
;     MEAN SIDEREAL time CORRECTLY ON THE J2000 COORDINATE SYSTEM.  THE
;     PREVIOUS VERSION WAS ONLY CORRECT IN THE B1950 COORDINATE SYSTEM.
; 3)  NOW RETURNS THE TRUE OF DATE GREENWICH SIDEREAL time ON THE J2000 SYS

PRO ic_grnwch_sidereal,orb_pos_time, gst

	half  = DOUBLE(0.50)
        C0    = DOUBLE(1.7533685592332653)     ;Polynomial Coef.
        C1    = DOUBLE(628.33197068884084)     ;Polynomial Coef.
        C2    = DOUBLE(0.67707139449033354E-05)  ;Polynomial Coef.
        C3    = DOUBLE(6.3003880989848915)     ;Polynomial Coef.
        C4    = DOUBLE(-0.45087672343186841E-09) ;Polynomial Coef.
        TWOPI = DOUBLE(6.283185307179586)      ;Two PI


        year = orb_pos_time(0)/1000
        day  = orb_pos_time(0) MOD 1000

;
;   Convert the given millisecond of day [orb_pos_time(1)] to second of day.
;
	secs = (DOUBLE(orb_pos_time(1)))/DOUBLE(1000.0)

;
;    Begin calculating the greenwich mean sidereal time **
;
	fday = secs/86400.00
	dj = DOUBLE(365*(year-1900)+(year-1901)/4+day-half)


;
;	      THE NEXT STATEMENT CAUSES THE REFERENCE EPOCH TO BE SHIFTED	
;	   TO THE J2000 REFERENCE EPOCH.
;

	T = (dj-DOUBLE(36525.0))/DOUBLE(36525.0)
        gst = DOUBLE(C0 + T*(C1 + T*(C2 + C4*T)) + C3*fday)
        gst = DOUBLE(gst MOD TWOPI)
        IF (gst LT DOUBLE(0.0)) THEN gst = gst + TWOPI

;
;   Convert gst to true of date by adjusting for nutation
;
	ic_get_nut_angles, T, deleps, delpsi, eps
	gst = gst + delpsi*cos(eps+deleps)

        END


       FUNCTION atan2d,x1,x2
          RETURN,DOUBLE(ATAN(x1,x2) * 57.2958)    
       END

      FUNCTION sind,x
         RETURN,DOUBLE(sin(x/57.2958))
      END

      FUNCTION cosd,x
         RETURN,DOUBLE(cos(x/57.2958))
      END

      FUNCTION dtand,x
          RETURN,DOUBLE(tan(x/57.2958))
      END

      FUNCTION datand,x
          RETURN,DOUBLE(atan(x)*57.2958)
      END




;-------------------------------------------------------------------------
;  ROUTINE:	fgeode
;  PURposE:	EVALUATE FUNCTIONAL SHAPE OF GEODE
;  FILE:	POINTING.FOR
; 
;  NAME	                TYPE    USE     DESCRIPTION
;  ----                 ----    ---     -----------
;  A			R*8	 I	EQUATORIAL RADIUS
;  B			R*8	 I	POLAR RADIUS
;  V1(3),V2(3)		R*8	 I	CARTESION COORDINATES
;  fgeode		R*8	 O	FUNCTIONAL SHAPE
;  			
;  AUTHOR	DATE	    DESCRIPTION OF CHANGE
;  ------	----	    ---------------------
;  G. GERMANY   3/31/95     ORIGINAL CODE
;
;-------------------------------------------------------------------------
      FUNCTION fgeode,a,b,v1,v2

      RETURN,v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2) * a*a/(b*b)
      END

;
;------------------------------------------------------------
; MODULE NAME: dfmag
; PURposE: CALCULATE MAGNITUDE OF 3-D VECTOR
;
; INPUTS:
; NAME           TYPE            DEFINITION
; ----           -----           ----------
; X1             REAL*8          X COMPONENT
; X2             REAL*8          Y COMPONENT
; X3             REAL*8          Z COMPONENT
;
; OUTPUTS:
; NAME           TYPE            DEFINITION
; ----           -----           ----------
; dfmag           REAL*8          MAGNITUDE
;
;------------------------------------------------------------

      FUNCTION dfmag,x1,x2,x3
         RETURN,SQRT(x1*x1+x2*x2+x3*x3)
      END


;-------------------------------------------------------------------------
;  ROUTINE:	uviptg
;  PURposE:	CALCULATE GEODETIC LAT,LON FOR A GIVEN PIXEL
;  FILE:	uviptg.FOR
; 
;  NAME	       DESCRIPTION
;  ----        -----------
;  cnv_flag    =1 on initial call then zero
;  system      =1 SEC, =2 PRIM
;  row,col     PIXEL SPECIFICATION
;  L0	       GCI LOOK DIRECTION
;  versStr     S/W Version String
;  att	       GCI attITUDE
;  orb	       GCI posITION
;  emis_hgt    LOS altitude
;  gdlat       GEODETIC LATITUDE
;  gdlon       GEODETIC LONGITUDE
;  rotm        conversion rotation matrix
;  ra          gci right ascension
;  dec         gci declination
;  			
;  FUNCTION/ROUTINE CALLS:	dcross_product
;				dunit
;				get_rotm
;				dmmult
;				get_scalar
;				ic_gci_to_geo
;				drtoll
;                               vector_to_ra_dec
;
;  AUTHOR	DATE	    DESCRIPTION OF CHANGE
;  ------	----	    ---------------------
;  G. GERMANY   7/31/95     ORIGINAL CODE
;  G. GERMANY   2/7/96      Fixed bug with pixel angles
;  G. GERMANY   2/13/96     Added cnv_flag & emis_hgt to call
;  G. GERMANY   2/28/96     Added ra & dec to call
;  G. GERMANY   3/29/96     Fixed bug with zrot,yrot orientation.
;                           Added versStr to call.
;  G. GERMANY   7/19/96     Fixed bug with cnv_flag call
;  G. GERMANY   7/25/96     Added trap for undefined rotm
;
;  NOTES:
;
;  1.     uviptg returns geodetic latitude and longitude
;         for a single pixel.
;
;  2.     The earth is assumed to be an ellipsoid of revolution.
;         See the routine fgeode for details.  Geodetic coordinates
;         are defined from the normal to the geode surface.  
;         Geographic (geocentric) coordinates assume the earth is
;         a sphere and are defined from the center of the sphere.
;
;  3.     The look direction for the specified pixel (Lpix) is
;         calculated from the look direction of the center of the 
;         UVI field of view (L0) by two successive rotations in
;         row and column directions.  Each pixel is assumed to have
;         a fixed angular width.  The angular distance from the center
;         of the pixel to the center of the fov is calculated and then
;         L0 is rotated into Lpix.  
;
;  4.     The secondary and primary detectors have different orientations
;         and require different rotations between L0 and Lpix.
;
;  5.     Geocentric lat/lon values are calculated for the intersection
;         of the look direction for the specified pixel (Lpix) and 
;         the surface of the earth.  The geocentric values are then
;         transformed into geodetic values.
;
;  6.     The intersection of Lpix and the earth is calculated first
;         in GCI coordinates and then converted to geographic coordinates.
;         The conversion is by means of ic_gci_to_geo.  This routine
;         and its supporting routines, was taken from the CDHF and is
;         part of the ICSS_TRANSF_orb call.
;
;  7.     The viewed emissions are assumed to originate emis_hgt km
;         above the surface of the earth.  See get_scalar for details.
;
;  8.     Set cnv_flag=1 for the first call of an image.  All subsequent
;         calls for the same image set cnv_flag=0.
;
;-------------------------------------------------------------------------

      PRO uviptg,cnv_flag,system,row,col,time,L0,att,orb,emis_hgt $
             ,versStr,gdlat,gdlon,rotm,ra,dec

;7/96
      versStr = 'V1.1c 7/96'

      Lpix = DBLARR(3)
      pos  = DBLARR(3)
      p_geo = DBLARR(3)
      dtemp = DBLARR(3)
      y_img_axis = DBLARR(3)
      m          = DBLARR(3,3)

;     time(1)=yyyyddd, time(2)=msec of day

      fov        = DOUBLE(8.0)
      fill_value = -1D31

      primary    = 1
      secondary  = 2
      nrows      = 228
      ncols      = 200

;...  single pixel angular resolution
      pc = DOUBLE(fov/ncols)
      pr = DOUBLE(fov/nrows)
      
 
;...  initialize output values to standard fill
      gdlat = fill_value
      gdlon = fill_value

; 2/28/96
      ra  = fill_value
      dec = fill_value

;...  check for invalid row/col values
;     return if encountered
      IF((row LT 1) OR (row GT nrows)  OR  $
         (col LT 1) OR (col GT ncols)) THEN RETURN


;...  calculate image y-axis in uvi coords.
;     This will be one of the axes of rotation
;     for the transform from L0 to Lpix below.
;     The other axis is the image z axis.
;     z = s/c spin axis
;     y = (s/c spin axis) vector cross (look direction)
;     x = look direction
      dtemp=DOUBLE(CROSSP(att,L0))
      dunit,dtemp,y_img_axis

;...  find rotation angles for specified pixel
      IF (system EQ secondary) THEN BEGIN
          zrot =  (col-100-0.5)*pc
          yrot = -(row-114-0.5)*pr
      ENDIF ELSE BEGIN
          IF(system EQ primary) THEN BEGIN
             zrot =  (col-100-0.5)*pc
;2/7/96      yrot = -(row-114-0.5)*pr
             yrot =  (row-114-0.5)*pr
          ENDIF ELSE BEGIN
;            error trap          
             RETURN
          ENDELSE
      ENDELSE

;3/96
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;4/22/97 M. Brittnacher--rotation offset correction
;      zrot = - zrot
;      yrot = - yrot

      zrotp = - zrot
      yrotp = - yrot

;...  rotate to actual spacecraft coordinates
;     presumed rotation offset is 2.5 degrees
      thetad = 2.5
      thetar = (!PI/180.)*thetad
      sinth = sin(thetar)
      costh = cos(thetar)
      zrot = zrotp*costh - yrotp*sinth
      yrot = zrotp*sinth + yrotp*costh

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

;...  rotate the central look vector (L0) to produce
;     the look vector for the current pixel (Lpix).
;     The get_rotm gets the rotation matrix for
;     a single rotation.  
;     Two rotations are performed - one in the row direction
;     (about z axis) and one in the column direction (about
;     y axis).
      get_rotm,att,zrot,m
      Lpix = TRANSPOSE(TRANSPOSE(m) # L0)
      Lpix = REFORM(Lpix)

      get_rotm,y_img_axis,yrot,m
      Lpix = TRANSPOSE(TRANSPOSE(m) # Lpix)
      Lpix = REFORM(Lpix)


;...  Convert look direction to unit vector
      dunit,Lpix,Lpix_out
      Lpix = Lpix_out

;...  2/28/96
;     calculate right ascension and declination
      vector_to_ra_dec,Lpix(0),Lpix(1),Lpix(2),ra,dec

;...  Find scalar (s) such that s*L0 points to
;     the imaged emission source.  If the line of
;     sight does not intersect the earth s=0.0
;2/96 get_scalar,orb,Lpix,s,f
      get_scalar,orb,Lpix,emis_hgt,s,f

      IF(s GT 0) THEN BEGIN

;...     Scale look vector & get position vector of imaged
;        source location (in GCI coords)
         pos = orb + S*Lpix

;...     Convert from GCI to GEO coordinates.  ROTM is the
;        rotation matrix.
;7/96  
         IF(N_ELEMENTS(rotm) EQ 0) THEN cnv_flag = 1
         IF(cnv_flag EQ 1) THEN ic_gci_to_geo,time,rotm 
         p_geo = TRANSPOSE(TRANSPOSE(rotm) # pos)

;...     Get geocentric lat/lon.  this converts from
;        a 3 element vector to two angles: lat & longitude
         drtoll,p_geo(0),p_geo(1),p_geo(2),gclat,gclon

;...     Convert to geodetic lat/lon.  F is the flattening
;        factor of the Earth.  See get_scalar for details.
;        Ref: Spacecraft Attitude Determination and Control,
;        J.R. Wertz, ed., 1991, p.821.
         gdlon=gclon
         IF(gclat GE 90)THEN BEGIN
            gdlat=DOUBLE(90)
         ENDIF ELSE BEGIN
            gdlat=datand(dtand(gclat)/(1-f)*(1-f))
         ENDELSE

      ENDIF

      END




;------------------------------------------------------------
; MODULE NAME: 	drtoll
; PURposE: 	CONVERT 3-D posITION vector TO LATITUDE/LONGITUDE
; FILE:		POINTING.FOR
;
;
; INPUTS:
; NAME           TYPE            DEFINITION
; ----           -----           ----------
; X              REAL*8          X COMPONENT (KM)
; Y              REAL*8          Y COMPONENT
; Z              REAL*8          Z COMPONENT
;
; OUTPUTS:
; NAME           TYPE            DEFINITION
; ----           -----           ----------
; LAT            REAL*8          GEO LATITUDE (DEG)
; LONG           REAL*8          GEO LONGITUDE
;
; NOTES:
;
;  The math function atan2d can handle x=0 but will fail
;      for both x and y zero.  This is trapped below.
;
;  This routine assumes the input vector (x,y,z) is in
;      geographic coordinates.
;
;------------------------------------------------------------

      PRO drtoll,x,y,z,lat,long

      IF((x EQ y) AND (x EQ 0)) THEN BEGIN
        lat  = DOUBLE(90*z/ABS(z))
        long = DOUBLE(0)
      ENDIF ELSE BEGIN
        lat  = atan2d(z,SQRT(x*x+y*y))
        long = atan2d(y,x)
      ENDELSE

      IF(long LT 0) THEN long=long+360

       RETURN
       END


;-------------------------------------------------------------------------
;  ROUTINE:	dunit
;  PURposE:	CONVERTS A VECTOR TO A UNIT VECTOR
;  FILE:	POINTING.FOR
; 
;  NAME	                TYPE    USE     DESCRIPTION
;  ----                 ----    ---     -----------
;  VECTOR(3)		R*8	 I	VECTOR TO BE CONVERTED
;  UVECTOR(3)           R*8      O      UNIT VECTOR
;  			
;  FUNCTION/ROUTINE CALLS:	DMAG
;
;  AUTHOR	DATE	    DESCRIPTION OF CHANGE
;  ------	----	    ---------------------
;  G. GERMANY   4/17/95     ORIGINAL CODE
;
;  NOTE: A temporary variable TMP is used for the calculation
;        and then copied into the output variable UVECTOR at the end.
;        This is to prevent inadvertently changing VECTOR if VECTOR
;        and UVECTOR are the same variable.
;
;-------------------------------------------------------------------------
      PRO dunit,vector,uvector
;-------------------------------------------------------------------------

      uvector = DBLARR(3)
      tmp     = DBLARR(3)
 
      dmag=dfmag(vector(0),vector(1),vector(2))

      IF (dmag NE 0) THEN BEGIN
         tmp=vector/dmag
      ENDIF ELSE BEGIN
         tmp=0.0
      ENDELSE

      uvector = tmp

      RETURN
      END

;
;-------------------------------------------------------------------------
;  ROUTINE:	get_rotm
;  PURposE:	CALCULATE ROTATION MATRIX ABOUT AN ARBITRARY axis
;  MODULE:	POINTING.FOR
; 
;  NAME	                TYPE    USE     DESCRIPTION
;  ----                 ----    ---     -----------
;  axis(3)		R*8   	 I	axis OF ROTATION
;  angle		R*8   	 I	angle OF ROTATION
;  M(3,3)		R*8	 O	ROTATION MATRIX
;  			
;  AUTHOR	DATE	    DESCRIPTION OF CHANGE
;  ------	----	    ---------------------
;  G. GERMANY   4/3/95      ORIGINAL CODE
;
;  REF: Spacecraft Attitude Determination and Control,
;       J.R.Wertz, ed., 1991, eq (E-6).
;-------------------------------------------------------------------------
      PRO get_rotm,axis,angle,m
;-------------------------------------------------------------------------


      m = DBLARR(3,3)

      s=DOUBLE(sind(angle))
      c=DOUBLE(cosd(angle))
      c1=1-c

      e1=axis(0)
      e2=axis(1)
      e3=axis(2)

      m(0,0) = c1*e1*e1 + c
      m(1,0) = c1*e1*e2 + e3*s
      m(2,0) = c1*e1*e3 - e2*s

      m(0,1) = c1*e1*e2 - e3*s
      m(1,1) = c1*e2*e2 + c
      m(2,1) = c1*e2*e3 + e1*s

      m(0,2) = c1*e1*e3 + e2*s
      m(1,2) = c1*e2*e3 - e1*s
      m(2,2) = c1*e3*e3 + c

      RETURN
      END



;-------------------------------------------------------------------------
;  ROUTINE:	get_scalar
;  PURposE:	RETURNS SCALAR S USED IN POINTING
;  FILE:	POINTING.FOR
; 
;  NAME	       DESCRIPTION
;  ----        -----------
;  O(3)	       orbITAL posITION vector
;  L(3)	       LOOK DIRECTION vector
;  emis_hgt    assumed height(km) of observed emissions
;  S	       SCALAR
;  F	       FLattENING FACTOR
;  			
;  FUNCTION/ROUTINE CALLS:
;
;  AUTHOR	DATE	    DESCRIPTION OF CHANGE
;  ------	----	    ---------------------
;  G. GERMANY   3/31/95     ORIGINAL CODE
;  G. GERMANY   2/13/96     made emis_hgt input
;
;-------------------------------------------------------------------------
      PRO get_scalar,O,L,emis_hgt,s,f
;-------------------------------------------------------------------------

;...  Equatoral radius (km) and polar flattening of the earth
;     Ref: Table 15.4, 'Explanatory Supplement to the
;          Astronomical Almanac,' K. Seidelmann, ed. (1992).
      re_eq = DOUBLE(6378.136)
      inv_f = DOUBLE(298.257)

;...  initialize output
      s=0.0

;...  get polar radius
      re_po=re_eq*(1-1./inv_f)

;...  get radii to assumed emission height
      ree = re_eq + emis_hgt
      rep = re_po + emis_hgt

;...  get flattening factor based on new radii
      f = (ree - rep)/ree

;...  get elements of quadratic formula
      a = fgeode(ree,rep,L,L)
      b = fgeode(ree,rep,L,O) * 2
      c = fgeode(ree,rep,O,O) - ree*ree

;...  check solutions to quadratic formula
      determinant = b*b - 4*a*c
      IF (determinant LT 0) THEN BEGIN
         s=0.0
         RETURN
      ENDIF

;...  solve quadratic formula (choose smallest solution)
      s1 = ( -b +  SQRT(determinant) ) / ( 2*a )
      s2 = ( -b -  SQRT(determinant) ) / ( 2*a )
      IF (s1 LT s2) THEN BEGIN
          s = s1
      ENDIF ELSE BEGIN
          s = s2
      ENDELSE

      RETURN
      END



      PRO vector_to_ra_dec,x,y,z,ra,dec

      IF((x EQ y) AND (x EQ 0)) THEN BEGIN
        dec = 90.*z/ABS(z)
        ra  = 0
      ENDIF ELSE BEGIN
        dec = atan2d(z,SQRT(x*x+y*y))
        ra  = atan2d(y,x)
      ENDELSE

      IF(ra LT 0) THEN ra=ra+360

      END

;------------------------------------------------------------
      FUNCTION atan2d,x1,x2
         RETURN,DOUBLE(ATAN(x1,x2) * 57.2958)    
      END
      FUNCTION sind,x
         RETURN,DOUBLE(sin(x/57.2958))
      END
      FUNCTION cosd,x
         RETURN,DOUBLE(cos(x/57.2958))
      END

;
;------------------------------------------------------------
; MODULE NAME: dfmag
; PURposE: CALCULATE MAGNITUDE OF 3-D VECTOR
;
; INPUTS:
; NAME           TYPE            DEFINITION
; ----           -----           ----------
; X1             REAL*8          X COMPONENT
; X2             REAL*8          Y COMPONENT
; X3             REAL*8          Z COMPONENT
;
; OUTPUTS:
; NAME           TYPE            DEFINITION
; ----           -----           ----------
; dfmag           REAL*8          MAGNITUDE
;
;------------------------------------------------------------

      FUNCTION dfmag,x1,x2,x3
         RETURN,SQRT(x1*x1+x2*x2+x3*x3)
      END

;-------------------------------------------------------------------------
;  ROUTINE:	dunit
;  PURposE:	CONVERTS A VECTOR TO A UNIT VECTOR
;  FILE:	POINTING.FOR
; 
;  NAME	                TYPE    USE     DESCRIPTION
;  ----                 ----    ---     -----------
;  VECTOR(3)		R*8	 I	VECTOR TO BE CONVERTED
;  UVECTOR(3)           R*8      O      UNIT VECTOR
;  			
;  FUNCTION/ROUTINE CALLS:	DMAG
;
;  AUTHOR	DATE	    DESCRIPTION OF CHANGE
;  ------	----	    ---------------------
;  G. GERMANY   4/17/95     ORIGINAL CODE
;
;  NOTE: A temporary variable TMP is used for the calculation
;        and then copied into the output variable UVECTOR at the end.
;        This is to prevent inadvertently changing VECTOR if VECTOR
;        and UVECTOR are the same variable.
;
;-------------------------------------------------------------------------
      PRO dunit,vector,uvector
;-------------------------------------------------------------------------

      uvector = DBLARR(3)
      tmp     = DBLARR(3)
 
      dmag=dfmag(vector(0),vector(1),vector(2))

      IF (dmag NE 0) THEN BEGIN
         tmp=vector/dmag
      ENDIF ELSE BEGIN
         tmp=0.0
      ENDELSE

      uvector = tmp

      RETURN
      END

;------------------------------------------------------------
      PRO vector_to_ra_dec,x,y,z,ra,dec

      IF((x EQ y) AND (x EQ 0)) THEN BEGIN
        dec = 90.*z/ABS(z)
        ra  = 0
      ENDIF ELSE BEGIN
        dec = atan2d(z,SQRT(x*x+y*y))
        ra  = atan2d(y,x)
      ENDELSE

      IF(ra LT 0) THEN ra=ra+360

      END 


;-------------------------------------------------------------------------
;  ROUTINE:	get_rotm
;  PURposE:	CALCULATE ROTATION MATRIX ABOUT AN ARBITRARY axis
;  MODULE:	POINTING.FOR
; 
;  NAME	                TYPE    USE     DESCRIPTION
;  ----                 ----    ---     -----------
;  axis(3)		R*8   	 I	axis OF ROTATION
;  angle		R*8   	 I	angle OF ROTATION
;  M(3,3)		R*8	 O	ROTATION MATRIX
;  			
;  AUTHOR	DATE	    DESCRIPTION OF CHANGE
;  ------	----	    ---------------------
;  G. GERMANY   4/3/95      ORIGINAL CODE
;
;  REF: Spacecraft Attitude Determination and Control,
;       J.R.Wertz, ed., 1991, eq (E-6).
;-------------------------------------------------------------------------
      PRO get_rotm,axis,angle,m
;-------------------------------------------------------------------------


      m = DBLARR(3,3)

      s=DOUBLE(sind(angle))
      c=DOUBLE(cosd(angle))
      c1=1-c

      e1=axis(0)
      e2=axis(1)
      e3=axis(2)

      m(0,0) = c1*e1*e1 + c
      m(1,0) = c1*e1*e2 + e3*s
      m(2,0) = c1*e1*e3 - e2*s

      m(0,1) = c1*e1*e2 - e3*s
      m(1,1) = c1*e2*e2 + c
      m(2,1) = c1*e2*e3 + e1*s

      m(0,2) = c1*e1*e3 + e2*s
      m(1,2) = c1*e2*e3 - e1*s
      m(2,2) = c1*e3*e3 + c

      RETURN
      END



;-------------------------------------------------------------------------
;  routine:  	uvilook
;  purpose:  	calculate central look direction of uvi field of view
;               (in gci coordinates).
; 
;  name	                type    use     description
;  ----                 ----    ---     -----------
; time(2)               i*4      i      time of image
; o_gci(3)              r*8      i      orbital position vector 
; a_gci3(3)             r*8      i      s/c attitude vector (gci)
; dsp_angle		r*8	 i	despun nadir offset angle
; filter                i*4      i      filter index
; versStr               str      o      sw version string
; l0			r*8	 o	look direction of uvi fov
; ra			r*8	 o	right ascension of look direction
; dec			r*8	 o	declination of look direction
;  			                        
;  function/routine calls:	dcross_product
;				dunit
;				dmmult
;				get_rotm
;				vector_to_ra_dec
;
;  author	date	    description of change
;  ------	----	    ---------------------
;  g. germany   4/3/96      original code
;  g. germany   4/9/96      fixed bug with dsp orientation
;                           added versStr output
;  g. germany   4/18/96     undid previous change to dsp orientation
;  g. germany   6/3/96      added filter input & cross_track corrections
;  g. germany   6/7/96      added mean cross_track correction
;  g. germany   8/13/96     changed filter indices
;  g. germany   1/97        added secondary cross track corrections
;
;  notes:
;
;	time(1) = 1000*year + day
;	time(2) = msec
;
;  filter = 0  BKG
;         = 1  1304
;         = 2  1356
;         = 3  LBHS
;         = 4  LBHL
;         = 5  SOLR
;
;  system  =1 PRI, =2 SEC 
;
;-------------------------------------------------------------------------

PRO uvilook,time,o_gci,a_gci3,dsp_angle,filter,versStr,l0,ra,dec $
                ,system=system

      versStr = 'v1.1  1/97'

;1/97
      primary    = 1
      secondary  = 2
      IF(NOT KEYWORD_SET(system)) THEN BEGIN
         MESSAGE,/INFO,'Assuming primary detector in use.'
      ENDIF ELSE BEGIN
         CASE system OF
          primary :
          secondary :
          ELSE : MESSAGE,'ERROR: Unsupported system value: ' $
                      +STRTRIM(STRING(system),2)
         ENDCASE
      ENDELSE 

;     calculate look direction in spacecraft frame (lsc).
;     defined as a unit vector positioned at the spacecraft,
;     perpendicular to the attitude vector, coplanar with
;     the attitude and position vectors, and pointing in 
;     the opposite sense of the position vector.  the
;     position vector is assumed to point from the center
;     of the earth to the spacecraft.  the look vector will
;     be converted to a unit vector and positioned at the
;     spacecraft after it has been converted to the image
;     reference frame.
;     lsc = (r x a) x a
      dtemp = DOUBLE(CROSSP(o_gci,a_gci3))
      lsc   = DOUBLE(CROSSP(dtemp,a_gci3))

;     calculate look direction in despun platform frame (ldsp).
;     defined as look direction in spaceraft frame rotated
;     about the z-axis by the despun offset angle.
      get_rotm,a_gci3,dsp_angle,m
      ldsp = TRANSPOSE(TRANSPOSE(m) # lsc)
      ldsp = REFORM(ldsp)

;     transform from dsp reference frame to uvi reference frame.
;     this corresponds to the alignment cube placement relative
;     to the dsp.  this is the mean pointing error relative to
;     the dsp, but does not include the filter-dependent pointing
;     error.
      ;conversion to change from pixels to degrees 
      pix_conv = 8.0D/228.

;1/97
      ;mean cross track error in pixels
      CASE system OF
       primary   : mean_cross_track = 7.0 
       secondary : mean_cross_track = 41.6  
       ELSE : MESSAGE,'ERROR: Unknown system value: '+STRTRIM(STRING(system),2)
      ENDCASE

      dtemp=DOUBLE(CROSSP(a_gci3,ldsp))
      dunit,dtemp,y_img_axis
      get_rotm,y_img_axis,mean_cross_track*pix_conv,m

      luvi = TRANSPOSE(TRANSPOSE(m) # ldsp)
      luvi = REFORM(luvi)


;     transform from uvi reference frame to image plane reference frame.
;     the corresponds to the look direction for the center of the image
;     plane and is different for each detector and filter combination.
      ;cross track error in pixels
;1/97
  
      CASE system OF
       primary   : BEGIN
                   CASE filter OF
                    0 :   cross_track =  0.0 ;bkg
                    1 :   cross_track =  0.9 ;1304
                    2 :   cross_track =  2.1 ;1356
                    3 :   cross_track =  0.0 ;LBHS
                    4 :   cross_track = -5.3 ;LBHL
                    5 :   cross_track =  0.0 ;SOLR
                    ELSE: cross_track =  0.0
                   ENDCASE
                   END
       secondary : BEGIN
                   CASE filter OF
                    0 :   cross_track =  0.0 ;bkg
                    1 :   cross_track =  0.7 ;1304
                    2 :   cross_track =  2.7 ;1356
                    3 :   cross_track =  0.0 ;LBHS
                    4 :   cross_track = -5.3 ;LBHL
                    5 :   cross_track =  0.0 ;SOLR
                    ELSE: cross_track =  0.0
                   ENDCASE
                   END
       ELSE : MESSAGE,'ERROR: Unknown system value: '+STRTRIM(STRING(system),2)
      ENDCASE
    
      dtemp=DOUBLE(CROSSP(a_gci3,luvi))
      dunit,dtemp,y_img_axis
      get_rotm,y_img_axis,cross_track*pix_conv,m

      limg = TRANSPOSE(TRANSPOSE(m) # luvi)
      limg = REFORM(limg)


;     convert look direction to unit vector
      dunit,limg,l0

;     calculate right ascension & declination (gci)
      vector_to_ra_dec,l0(0),l0(1),l0(2),ra,dec 

END
;--------------------------------------------------------------

PRO get_dsp_angle,year,doy,hour,minute,second,dsp_angle,dsp_flag,path

RESTORE, path+'dsp_angle.sav'

siz = SIZE(yeaArr)
arrSiz = siz(1)

; verify that data is available for this year
IF (year LT yeaArr(0)) OR (year GT yeaArr(arrSiz-1)) THEN BEGIN
  PRINT,'DSP angles not available for this year'
  RETURN
ENDIF

; search for lower and upper indices for given day of Polar
dop = CEIL((year - 1996.)*365.25) + doy
lindx = 0
uindx = arrSiz-1
indx = FLOOR(0.5*(lindx + uindx))
savindx = -1
WHILE (dopArr(indx) NE dop) DO BEGIN
  IF (dopArr(indx) GT dop) THEN uindx = indx ELSE lindx = indx
  indx = FLOOR(0.5*(lindx+uindx))
  IF (indx EQ savindx) THEN BEGIN
    err0 = WIDGET_MESSAGE(/ERROR,'Pitch angles not available for this day')
    dsp_flag = 1
    RETURN
  ENDIF
  savindx = indx
ENDWHILE
savindx = indx
WHILE (doyArr(indx) EQ doy) DO indx = indx-1
lindx = indx+1
indx = savindx
IF (doy EQ doyArr(arrSiz-1)) THEN uindx = arrSiz-1 $
ELSE BEGIN
  WHILE (doyArr(indx) EQ doy) DO indx = indx+1
  uindx = indx-1
ENDELSE

; search for the dsp slew angle for the given time of day
sod = hour*3600L + minute*60L + second*1L
CASE 1 OF

  (sod LT sodArr(lindx)): BEGIN
    IF (dopArr(lindx-1) EQ dop-1) THEN $ 
      dsp_angle = dspArr(lindx-1) $
    ELSE BEGIN
      err0 = WIDGET_MESSAGE(/ERROR,'Pitch angles for previous day required')
      dsp_flag = 1
      RETURN
    ENDELSE
                          END

  (sod GT sodArr(uindx)): BEGIN
    IF (sod LE 86400L) THEN dsp_angle = dspArr(uindx) $
    ELSE BEGIN
      err0 = WIDGET_MESSAGE(/ERROR,'Time exceeds length of day')
      dsp_flag = 1
      RETURN
    ENDELSE
                          END

  ELSE: BEGIN
    indx = lindx
    WHILE (sod GT sodArr(indx)) DO indx = indx+1
    dsp_angle = dspArr(indx-1)
        END

ENDCASE

dsp_flag = 0  ; no errors

END

PRO find_record,flag,id,epoch,record

  record=-1
  cdf_struct = CDF_INQUIRE(id)
;cdf_struct.maxrec is the address (from 0) of the last record
;therefore, the total number of records is maxrec+1
  numrecs=cdf_struct.maxrec+1
  epoch_ary = DBLARR(numrecs)
  IF(flag EQ 'orb') THEN BEGIN
    CDF_VARGET,id,'Epoch',tmp,REC_COUNT=numrecs,COUNT=[1]
  ENDIF ELSE BEGIN
    CDF_VARGET,id,'EPOCH',tmp,REC_COUNT=numrecs,COUNT=[1,1,1]
  ENDELSE
  epoch_ary = REFORM(tmp)
  diff = ABS(epoch_ary - epoch)
  r_ary = WHERE(diff EQ MIN(diff))
  record=r_ary(0)

END

pro get_mlat, fin, fsize, glat, glon, mlat, mlon
 
  glatr = 0.25*float(round(glat/0.25))
  if glon lt 0.0 then glon = glon + 360.0
  glonr = 0.25*float(round(glon/0.25))
  dlat = long(4.0*(90.0-glatr))
  dlon = long(4.0*glonr)
  total = dlat*(360*4) + dlon
  total = long(16)*total
  if (total le fsize-16) then begin
    point_lun,fin,total
    readu,fin,glat,glon,mlat,mlon
    mlon = mlon-12.0
    if mlon lt 0.0 then mlon = mlon + 24.0
  endif else begin
    mlat = -999.0
    mlon = -999.0
  endelse
 
  return
 
end

pro strip

  path = '/d/ridley/d.data/d.polar/d.software/d.savefiles/'

  dirname = ''
  print, 'Enter directory of UVI cdf files (like /l2/polar_data/) :'
  read, dirname
  if strmid(dirname,strlen(dirname)-1,1) ne '/' then dirname = dirname+'/'

  date = ''
  print, 'Enter date in form yymmdd (as in 960724):
  read, date

  atfile = dirname+'po_at_def_19'+date+'_v01.cdf'
  l1file = dirname+'po_l1_uvi_19'+date+'_v01.cdf'
  orfile = dirname+'po_or_def_19'+date+'_v01.cdf'

  at_id = cdf_open(atfile)
  l1_id = cdf_open(l1file)
  or_id = cdf_open(orfile)

  info = cdf_inquire(l1_id)
  nrecs = info.maxrec

  ftime = intarr(2,6)

  CDF_VARGET1,l1_id,'EPOCH',epoch,REC_START=0
  CDF_EPOCH,epoch,year,month,day,hour,minute $ 
    ,second,millisecond,/BREAKDOWN_EPOCH     
  ftime(0,0) = year
  ftime(0,1) = month
  ftime(0,2) = day
  ftime(0,3) = hour
  ftime(0,4) = minute
  ftime(0,5) = second

  CDF_VARGET1,l1_id,'EPOCH',epoch,REC_START=nrecs-1
  CDF_EPOCH,epoch,year,month,day,hour,minute $ 
    ,second,millisecond,/BREAKDOWN_EPOCH     
  ftime(1,0) = year
  ftime(1,1) = month
  ftime(1,2) = day
  ftime(1,3) = hour
  ftime(1,4) = minute
  ftime(1,5) = second

  c_a_to_s, ftime(0,*), stime
  c_a_to_s, ftime(1,*), etime

  print, 'Start time : ',stime
  print, 'Enter start time to start stripping images (in same format) :'
  srtime = ''
  read, srtime
  c_s_to_a, itime, srtime
  cdf_epoch,sepoch,itime(0), itime(1), itime(2), itime(3), itime(4),	$
	itime(5), /compute

  print, 'End time   : ',etime
  print, 'Enter end time to end stripping images (in same format) :'
  ertime = ''
  read, ertime
  c_s_to_a, itime, ertime
  cdf_epoch,eepoch,itime(0), itime(1), itime(2), itime(3), itime(4),	$
	itime(5), /compute

  print, ''
  print, '1. 1304'
  print, '2. 1356'
  print, '3. LBHS'
  print, '4. LBHL'
  print, '5. Solar'
  print, 'Enter filter number which you would like : '
  filter = 0
  read, filter

  print, 'Enter number of times to skip in between each image : '
  print, '(0 for no skipping, 1 for every other image, ...)'
  nskip = 0
  read, nskip
  ntimes = nskip + 1

  print, 'You can save a bunch of time by stating that the center'
  print, 'pixel has to map to above a certain latitude (ex. 0 - for north)'
  print, 'If you want all images, enter -90.0'
  print, 'Enter minimum center latitute :'
  minlat = 0.0
  read, minlat

  print, 'You can also include CGM coordinates in the file. This will not'
  print, 'affect the plotting routines which have been created. It will'
  print, 'take about 30 seconds more per image.'
  print, 'Would you like cgm coordinates (y/n) ?'
  cgm = ''
  read, cgm
  if (strmid(mklower(cgm),0,1) eq 'y') then cgm = 1 else cgm = 0

  fac = 1

  srec = -1
  erec = -1

  i = 0
  done = 0
  nimages = 0
  while not done do begin
    CDF_VARGET1,l1_id,'EPOCH',epoch,REC_START=i
    if (epoch ge sepoch) and (srec eq -1) then srec = i
    if ((epoch gt eepoch) or (i ge nrecs-1)) and (erec eq -1) then begin
      erec = i-1
      done = 1
    endif
    if (srec ge 0) and (erec eq -1) then begin
      CDF_VARGET1,l1_id,'FILTER',ffilter,REC_START=i
      if ffilter eq filter then begin
        if nimages eq 0 then 						$
	  image_rec = i							$
	else image_rec = [image_rec,i]
	nimages = nimages + 1
      endif
    endif
    i = i + 1
  endwhile

  print, 'There are ',tostr(erec-srec), ' images between those two times.'
  print, tostr(nimages),' of those were the taken filter which you selected.'

; from studio.pro

  RESTORE, path+'cal_pri.sav'
  RESTORE, path+'cal_sec.sav'
  RESTORE, path+'f13_pri.sav'
  RESTORE, path+'f13_sec.sav'

; -----

  for n=0,nimages-1,ntimes do begin

    print, 'Working on image number : ',tostr(n+1),' out of ',tostr(nimages)

    CDF_VARGET1,l1_id,'EPOCH',epoch,REC_START=image_rec(n)
    CDF_EPOCH,epoch,year,month,day,hour,minute, $ 
      second,milli,/BREAKDOWN_EPOCH     
    doy=jday(year,month,day)

    yr = year
    if yr gt 1900 then yr = yr - 1900
    if yr gt 100 then yr = yr - 100
    time_tag = chopr('0'+tostr(yr),2) +					$
	       chopr('0'+tostr(month),2) +				$
	       chopr('0'+tostr(day),2) + '_' +				$
	       chopr('0'+tostr(hour),2) +				$
	       chopr('0'+tostr(minute),2) +				$
	       chopr('0'+tostr(second),2)

; taken from studio.pro :

    CDF_VARGET,l1_id,'INT_IMAGE',image,REC_START=image_rec(n)

    CDF_VARGET1,l1_id,'SYSTEM',tmp,REC_START=image_rec(n)
    system = tmp + 1

    CASE system OF
      1: BEGIN
           cal_img = cal_pri
	   flat_field = f13_pri
         END
      2: BEGIN
           cal_img = cal_sec
	   flat_field = f13_sec
         END
    ENDCASE

    ; initialize calibration units and string

    background = 29.0
    cal_unit = 0

    CASE cal_unit OF
    0: BEGIN
         units = 30.0
         unitstr = 'Photons cm!U-2!N s!U-1!N'
       END
    1: BEGIN
         units = 30.166
         unitstr = '       Rayleighs'
       END
    ENDCASE

    enhance = 100.0
    ; obtain integration period
    CDF_VARGET1,l1_id,'FRAMERATE',rate,REC_START=image_rec(n)
    ; correct for integration period
    ip_corr = 4./rate
    CDF_VARGET1,l1_id,'GAIN',gain,REC_START=image_rec(n)
    CDF_VARGET1,l1_id,'DOORPOS',doorpos,REC_START=image_rec(n)
    CASE doorpos OF
      1: door_trans = 1.          ; full extent
      2: door_trans = 1.          ; open
      3: BEGIN                  ; closed
           CASE filter OF
             1: door_trans = 0.42
             2: door_trans = 0.42
             3: door_trans = 0.50
             4: door_trans = 0.54
             5: door_trans = 0.60
           ENDCASE
         END
    ENDCASE

;   get platform lock flag
    lock = 0B
    CDF_VARGET1,l1_id,'PFLCK',lock,REC_START=image_rec(n)
    IF (lock) THEN $
      lokw = WIDGET_MESSAGE('WARNING: PLATFORM NOT LOCKED. MAPPING INCORRECT')

    img_factor = 0.002*enhance*units*ip_corr/door_trans
    new_img = img_factor*flat_field*cal_img(gain-1,filter-1)*(image-background)
    image = new_img

    get_dsp_angle,year,doy,hour,minute,second,dsp_angle,dsp_flag,path

    find_record,'orb',or_id,epoch,orb_record
    CDF_VARGET,or_id,'GCI_POS',o_gci,REC_START=orb_record,COUNT=[3]
    find_record,'orb',at_id,epoch,att_record
    CDF_VARGET1,at_id,'GCI_R_ASCENSION',ra,REC_START=att_record
    CDF_VARGET1,at_id,'GCI_DECLINATION',dec,REC_START=att_record
    a_gci3 = fltarr(3)
    a_gci3(0)=COS(dec)*COS(ra) ;x
    a_gci3(1)=COS(dec)*SIN(ra) ;y
    a_gci3(2)=SIN(dec)         ;z

    time = LONARR(2)
    time(0) = 1000*year + doy
    time(1) = (hour*3600 + minute*60 + second)*1000
    radE = 6357.
    rdist = NORM(o_gci)/radE
    satp = [rdist,0,0]

    ic_gci_to_geo,time,rotm 
    o_geo = TRANSPOSE(TRANSPOSE(rotm) # o_gci)
    a_geo = TRANSPOSE(TRANSPOSE(rotm) # a_gci3)

    ; Get geocentric lat/lon for the spacecraft position

    drtoll,o_geo(0),o_geo(1),o_geo(2),gclat,gclon
    p0lat = gclat
    p0lon = gclon

    ; Transform to satellite coordinate system:
    ; origin at center of Earth, satellite is on z axis
    deg = 180./acos(-1.0)
    rad = acos(-1.0)/180.
    ph = (gclon+90)*rad
    th = (90-gclat)*rad
    cosph = cos(ph)
    sinph = sin(ph)
    costh = cos(th)
    sinth = sin(th)
    geo2sat = [[cosph,sinph,0], $
              [-costh*sinph,costh*cosph,sinth], $
              [sinth*sinph,-sinth*cosph,costh]]

    o_sat = TRANSPOSE(geo2sat ## [o_geo(0),o_geo(1),o_geo(2)])
    z_sat = TRANSPOSE(geo2sat ## [0,0,1.])
    a_sat = TRANSPOSE(geo2sat ## [a_geo(0),a_geo(1),a_geo(2)])
    rotation = atan(-a_sat(0),a_sat(1))*deg

    ; Obtain boresight direction in RA/DEC
    uvilook,time,o_gci,a_gci3,dsp_angle,filter,versStr,l0,ra,dec,system=system

; ------

    cnv_flag = 1
    emis_hgt = 120.0

    xs = n_elements(image(*,0))/fac
    ys = n_elements(image(0,*))/fac
    lats = fltarr(xs,ys)
    lons = fltarr(xs,ys)
    imas = fltarr(xs,ys)
    mlat = fltarr(xs,ys)
    mlt  = fltarr(xs,ys)

    uviptg,cnv_flag,system,ys*fac/2,xs*fac/2,time,l0,a_gci3,o_gci,emis_hgt, $
           versStr,lat_x,lon_x,rotm,ra,dec

    print, 'Center pixel maps to ',lat_x,' deg. latitude'

    if lat_x gt minlat then begin

      fname = time_tag
      if (filter eq 1) then fname = fname+'_1304'
      if (filter eq 2) then fname = fname+'_1356'
      if (filter eq 3) then fname = fname+'_lbhs'
      if (filter eq 4) then fname = fname+'_lbhl'
      if (filter eq 5) then fname = fname+'_solar'
      if (cgm) then fname = fname + '.cgm' else fname = fname + '.bin'
      print, 'writing file : ',fname
      openw,1,fname

      printf,1, xs, ys
      printf,1, time_tag, ip_corr

      for i=1,xs do begin
        for j=1,ys do begin

          x1 = (i-1)*fac
          y1 = (j-1)*fac

          if image(x1,y1) gt 0 then begin

            uviptg,cnv_flag,system,y1,x1,time,l0,a_gci3,o_gci,emis_hgt, $
              versStr,lat_x,lon_x,rotm,ra,dec

            cnv_flag = 0

          endif else begin

            lat_x = 0.0
            lon_x = 0.0
            image(x1,y1) = -1.0

          endelse

          if abs(lat_x) gt 90.0 then lat_x = 0.0
          if abs(lon_x) gt 360.0 then lon_x = 0.0

	  lats(i-1,j-1) = lat_x
          lons(i-1,j-1) = lon_x
	  imas(i-1,j-1) = image(x1,y1)

        endfor

      ;print,'Row ',tostr(i),' of ',tostr(xs),' completed.'

      endfor

      writeu,1, lats
      writeu,1, lons
      writeu,1, imas

      if cgm then begin

	print, 'converting file to cgm coordinates'

        fin = 12
        close,fin
        openu,fin,path+'cgm_bin.dat'
        a = fstat(fin)
        fsize = a.size
        dx = 0.25

        for i=0,xs-1 do for j=0,ys-1 do 			$
          if imas(i,j) gt 0.0 then begin
            get_mlat, fin, fsize, lats(i,j),lons(i,j), ml, mm
	    mlat(i,j) = ml
	    mlt(i,j)  = mm
	  endif

	writeu,1,mlat
	writeu,1,mlt

      endif

      close,1
      print, 'Wrote file.'

    endif else print, 'Image center maps to too low latitude, skipping'

  endfor

  cdf_close, at_id
  cdf_close, l1_id
  cdf_close, or_id

  return

end

strip

end
