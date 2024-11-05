      SUBROUTINE ENERGY(PS,U,V,NX,NY,SL,WLAT,SINLAT,COSLAT,LEVS)
      IMPLICIT NONE

      INTEGER NX, NY, LEVS, TOP100
! ==>   REAL PS(NX,NY), U(NX,NY,LEVS), V(NX,NY,LEVS) ==> add psib array
      REAL PS(NX,NY), U(NX,NY,LEVS), V(NX,NY,LEVS)
      REAL PSIB(NX,NY)
      REAL*4 ALAND(NX,NY)
      INTEGER*4 LEN1, LEN2
      INTEGER LAND(NX,NY)
      REAL SL(LEVS), W(LEVS), WLAT(NY), SINLAT(NY), COSLAT(NY)
      REAL SINX(NX), COSX(NX), PI
      REAL UAVE(NX,NY), VAVE(NX,NY)
      REAL U100AVE(NX,NY), V100AVE(NX,NY)
      REAL UC(NX,NY), US(NX,NY), VC(NX,NY), VS(NX,NY)
      REAL U100C(NX,NY), U100S(NX,NY), V100C(NX,NY), V100S(NX,NY)
      REAL XP1(NY), XP2(NY), XP3(NY), XW1(NY), XW2(NY), XW3(NY)
      REAL XPIB1(NY), XPIB2(NY), XPIB3(NY)
      REAL X100W1(NY), X100W2(NY), X100W3(NY)
      REAL P1, P2, P3, W1, W2, W3, PIB1, PIB2, PIB3
      REAL W1001, W1002, W1003
      REAL XP1NH, XP2NH, XP3NH
      REAL XP1SH, XP2SH, XP3SH
      REAL XPIB1NH, XPIB2NH, XPIB3NH
      REAL XPIB1SH, XPIB2SH, XPIB3SH
      REAL XW1NH, XW2NH, XW3NH
      REAL XW1SH, XW2SH, XW3SH
      REAL X100W1NH, X100W2NH, X100W3NH
      REAL X100W1SH, X100W2SH, X100W3SH

!C==>  add  REAL ib terms


      REAL VW, NHW, SHW, T100W, TVW
      REAL SUM, SUM2, SUM3, SUMPT, MEAN_O_P
      REAL R, C, CMA, G, O
      REAL CW1, CP1, CW2, CP2, CW3, CP3
      INTEGER I, J, K, IDEBUG
      DATA IDEBUG/1/

      IF (IDEBUG.GE.1) THEN
         PRINT *, 'NX = ',NX,'  NY = ',NY,'  LEVS = ', LEVS
      ENDIF

      DO I = 1, LEVS
         IF (SL(I).GT.0.1) TOP100 = I
      ENDDO
      PRINT *, 'TOP100 = ', TOP100

!     print *, 'clat = ', coslat
!     print *, 'slat = ', sinlat
!     print *, 'wlat = ', wlat

!     DO K = 1, LEVS
!        print *, 'levs = ', k
!        print *,'u',' max min:',maxval(u(:,:,k)),minval(u(:,:,k))
!        print *,'v',' max min:',maxval(v(:,:,k)),minval(v(:,:,k))
!     ENDDO

!C
!C  READ LAND SEA MASK
!C
!c      IF (NX.EQ.512) THEN
!c         OPEN(21,FILE='land_mask_high.txt',STATUS='OLD')
!c      ELSE IF (NX.EQ.768) THEN
!c         OPEN(21,FILE='land_mask_veryhigh.txt',STATUS='OLD')
!c      ELSE IF (NX.EQ.1152) THEN
!c         OPEN(21,FILE='land_mask_ultrahigh.txt',STATUS='OLD')
!c      ELSE IF (NX.EQ.1760) THEN
!c         OPEN(21,FILE='land_mask_t574.txt',STATUS='OLD')
!c      ENDIF

      OPEN(21,FILE='landmask_1536.dat',FORM='UNFORMATTED',STATUS='OLD')
      READ(21) ALAND
!     PRINT *, 'ALAND'
!     PRINT *, ALAND(1,1:NY)

!C
!C     CONVERT SURFACE PRESSURE TO PASCAL if it is not in Pascal
!C
      print *, ' PS center point = ', PS(nx/2, ny/2)
      DO J = 1, NY
         DO I = 1, NX
!           PS(I,J) = PS(I,J) * 100.0
            LAND(I,J) = ALAND(I,J)
         ENDDO
      ENDDO

!C==>  add a subroutine here to calculate mean ocean pressure and  ib pressure
!C
!C     CALCULATE MEAN OCEAN PRESSURE
!C
      SUM2 = 0.0
      SUM3 = 0.0
      DO J = 1, NY
         SUM = 0.0
         SUMPT = 0.0
         DO I = 1, NX
            SUM = SUM + PS(I,J) * ( 1.0 - LAND(I,J) )
            SUMPT = SUMPT + ( 1.0 - LAND(I,J) )
         ENDDO
         SUM2 = SUM2 + SUM * COSLAT(J)
         SUM3 = SUM3 + SUMPT * COSLAT(J)
      ENDDO
      MEAN_O_P = SUM2 / SUM3

!C
!C     CALCULATE IB PRESSURE
!C
      DO J = 1, NY
         DO I = 1, NX
         PSIB(I,J)=PS(I,J)*LAND(I,J)+MEAN_O_P*(1.0-LAND(I,J))
         ENDDO
      ENDDO

!C
!C     CALCULATE SINX AND COSX FOR ALL NX POINTS
!C
      PI = ATAN(1.0) * 4.0
      DO I = 1, NX
         SINX(I)=SIN(2.0*PI*(I-1)/FLOAT(NX))
         COSX(I)=COS(2.0*PI*(I-1)/FLOAT(NX))
      ENDDO
      IF (IDEBUG.GE.2) THEN
         PRINT *,'    I      SINX      COSX'
         DO I = 1, NX
            PRINT 110, I, SINX(I), COSX(I), (SINX(I)**2+COSX(I)**2)
  110       FORMAT(I5,3F10.6)
         ENDDO
      ENDIF

!C
!C     CALCULATE VERTICAL WEIGHT
!C
      DO K = 2, LEVS - 1
         W(K) = 0.5*(SL(K-1) - SL(K+1))
      ENDDO
      W(LEVS) = 0.5*(SL(LEVS-1) - SL(LEVS)) + SL(LEVS)
      W(1) = (1.0 - SL(1)) + 0.5*(SL(1) - SL(2))

!     print *, ' w = ', w

      VW = 0.0
      DO K = 1, LEVS
        VW = VW + W(K)
      ENDDO
      IF (IDEBUG.GE.1) THEN
         PRINT *, 'TOTAL VERTICAL WEIGHT IS ', VW
!        PRINT 111, ((SL(I), W(I)), I = LEVS, 1, -1)
  111    FORMAT(F8.4,F10.5)
      ENDIF

      IF (IDEBUG.GE.1) THEN
         SUM = 0.0
         SUM2 = 0.0
         DO J = 1, NY
            SUM = SUM + WLAT(J)
            SUM3 = SINLAT(J)**2 + COSLAT(J)**2
!           PRINT 112, J, WLAT(J), SUM3
  112       FORMAT(I4,5E15.4)
         ENDDO
         PRINT *, 'TOTAL LATITUDE WEIGHT IS ', SUM
      ENDIF

!C
!C     CALCULATE VERTICAL AVERAGE FOR U AND V
!C
      DO J = 1, NY
         DO I = 1, NX
!C
!C     AVERAGE U AND V FOR THE WHOLE COLUMN
!C
            UAVE(I,J) = 0.0
            VAVE(I,J) = 0.0
            TVW = 0.0
            DO K = 1, LEVS
               UAVE(I,J) = UAVE(I,J) + U(I,J,K)*W(K)
!C              IF ((K.NE.52).AND.(K.NE.63)) THEN
                  VAVE(I,J) = VAVE(I,J) + V(I,J,K)*W(K)
!C              ENDIF
               TVW = TVW + W(K)
            ENDDO
            UAVE(I,J) = UAVE(I,J) / TVW
            VAVE(I,J) = VAVE(I,J) / TVW
!C
!C     AVERAGE U AND V UP TO 100MB (SIGMA LEVEL 32)
!C
            T100W = 0.0
            U100AVE(I,J) = 0.0
            V100AVE(I,J) = 0.0
            DO K = 1, TOP100
               U100AVE(I,J) = U100AVE(I,J) + U(I,J,K)*W(K)
               V100AVE(I,J) = V100AVE(I,J) + V(I,J,K)*W(K)
               T100W = T100W + W(K)
            ENDDO
            U100AVE(I,J) = U100AVE(I,J) / T100W
            V100AVE(I,J) = V100AVE(I,J) / T100W
         ENDDO
      ENDDO

!C
!C     WEIGHT VERTICAL AVERAGE BY SURFACE PRESSURE
!C
      DO J = 1, NY
         DO I = 1, NX
            UAVE(I,J) = UAVE(I,J) * PS(I,J)
            VAVE(I,J) = VAVE(I,J) * PS(I,J)
            U100AVE(I,J) = U100AVE(I,J) * PS(I,J)
            V100AVE(I,J) = V100AVE(I,J) * PS(I,J)
         ENDDO
      ENDDO

!C
!C     PRE CALCULATE THE TERM : U * COS(PHI), U * SIN(PHI)
!C                              V * COS(LAMBDA), V * SIN(LAMBDA)
!C
      DO J = 1, NY
         DO I = 1, NX
            UC(I,J) = UAVE(I,J) * COSLAT(J)
            US(I,J) = UAVE(I,J) * SINLAT(J)
            VC(I,J) = VAVE(I,J) * COSX(I)
            VS(I,J) = VAVE(I,J) * SINX(I)
            U100C(I,J) = U100AVE(I,J) * COSLAT(J)
            U100S(I,J) = U100AVE(I,J) * SINLAT(J)
            V100C(I,J) = V100AVE(I,J) * COSX(I)
            V100S(I,J) = V100AVE(I,J) * SINX(I)
         ENDDO
      ENDDO

!C
!C     CALCULATE ZONAL AVERAGE
!C
      DO J = 1, NY
         P1 = 0.0
         P2 = 0.0
         P3 = 0.0
         PIB1 = 0.0
         PIB2 = 0.0
         PIB3 = 0.0
         W1 = 0.0
         W2 = 0.0
         W3 = 0.0
         W1001 = 0.0
         W1002 = 0.0
         W1003 = 0.0
         DO I = 1, NX
            P1 = P1 + PS(I,J) * COSX(I)
            P2 = P2 + PS(I,J) * SINX(I)
            P3 = P3 + PS(I,J)
            W1 = W1 + US(I,J) * COSX(I) - VS(I,J)
            W2 = W2 + US(I,J) * SINX(I) + VC(I,J)
            W3 = W3 + UC(I,J)
!C==>  add psib calculations
            PIB1 = PIB1 + PSIB(I,J) * COSX(I)
            PIB2 = PIB2 + PSIB(I,J) * SINX(I)
            PIB3 = PIB3 + PSIB(I,J)
            W1001 = W1001 + U100S(I,J) * COSX(I) - V100S(I,J)
            W1002 = W1002 + U100S(I,J) * SINX(I) + V100C(I,J)
            W1003 = W1003 + U100C(I,J)
         ENDDO
         XP1(J) = P1 / FLOAT(NX)
         XP2(J) = P2 / FLOAT(NX)
         XP3(J) = P3 / FLOAT(NX)
         XPIB1(J) = PIB1 / FLOAT(NX)
         XPIB2(J) = PIB2 / FLOAT(NX)
         XPIB3(J) = PIB3 / FLOAT(NX)
         XW1(J) = W1 / FLOAT(NX)
         XW2(J) = W2 / FLOAT(NX)
         XW3(J) = W3 / FLOAT(NX)
         X100W1(J) = W1001 / FLOAT(NX)
         X100W2(J) = W1002 / FLOAT(NX)
         X100W3(J) = W1003 / FLOAT(NX)
         IF (IDEBUG.EQ.2) THEN
!C           PRINT *, '   J            XP1           XP2            XP3'
            PRINT 112, J, XP1(J), XP2(J), XP3(J)
         ELSE IF (IDEBUG.EQ.3) THEN
            PRINT 112, J, XW1(J), XW2(J), XW3(J)
         ENDIF
      ENDDO
!C
!C     CALCULATE HEMISPHERIC AVERAGE, NORTHERN HEMISHPERE FIRST
!C
      XP1NH = 0.0
      XP2NH = 0.0
      XP3NH = 0.0
      XPIB1NH = 0.0
      XPIB2NH = 0.0
      XPIB3NH = 0.0
      XW1NH = 0.0
      XW2NH = 0.0
      XW3NH = 0.0
      X100W1NH = 0.0
      X100W2NH = 0.0
      X100W3NH = 0.0
      NHW = 0.0
      DO J = 1, NY/2
         XP1NH = XP1NH + XP1(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XP2NH = XP2NH + XP2(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XP3NH = XP3NH + XP3(J) * COSLAT(J)**2 * WLAT(J)
         XPIB1NH = XPIB1NH + XPIB1(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XPIB2NH = XPIB2NH + XPIB2(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XPIB3NH = XPIB3NH + XPIB3(J) * COSLAT(J)**2 * WLAT(J)
         XW1NH = XW1NH + XW1(J) * WLAT(J)
         XW2NH = XW2NH + XW2(J) * WLAT(J)
         XW3NH = XW3NH + XW3(J) * WLAT(J)
         X100W1NH = X100W1NH + X100W1(J) * WLAT(J)
         X100W2NH = X100W2NH + X100W2(J) * WLAT(J)
         X100W3NH = X100W3NH + X100W3(J) * WLAT(J)
         NHW = NHW + WLAT(J)
      ENDDO
      XP1NH = XP1NH / NHW
      XP2NH = XP2NH / NHW
      XP3NH = XP3NH / NHW
      XPIB1NH = XPIB1NH / NHW
      XPIB2NH = XPIB2NH / NHW
      XPIB3NH = XPIB3NH / NHW
      XW1NH = XW1NH / NHW
      XW2NH = XW2NH / NHW
      XW3NH = XW3NH / NHW
      X100W1NH = X100W1NH / NHW
      X100W2NH = X100W2NH / NHW
      X100W3NH = X100W3NH / NHW
!C
!C     NOW CALCULATE SOUTHERN HEMISHPERE
!C
      XP1SH = 0.0
      XP2SH = 0.0
      XP3SH = 0.0
      XPIB1SH = 0.0
      XPIB2SH = 0.0
      XPIB3SH = 0.0
      XW1SH = 0.0
      XW2SH = 0.0
      XW3SH = 0.0
      X100W1SH = 0.0
      X100W2SH = 0.0
      X100W3SH = 0.0
      SHW = 0.0
      DO J = NY/2+1, NY
         XP1SH = XP1SH + XP1(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XP2SH = XP2SH + XP2(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XP3SH = XP3SH + XP3(J) * COSLAT(J)**2 * WLAT(J)
         XPIB1SH = XPIB1SH + XPIB1(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XPIB2SH = XPIB2SH + XPIB2(J) * SINLAT(J) * COSLAT(J) * WLAT(J)
         XPIB3SH = XPIB3SH + XPIB3(J) * COSLAT(J)**2 * WLAT(J)
         XW1SH = XW1SH + XW1(J) * WLAT(J)
         XW2SH = XW2SH + XW2(J) * WLAT(J)
         XW3SH = XW3SH + XW3(J) * WLAT(J)
         X100W1SH = X100W1SH + X100W1(J) * WLAT(J)
         X100W2SH = X100W2SH + X100W2(J) * WLAT(J)
         X100W3SH = X100W3SH + X100W3(J) * WLAT(J)
         SHW = SHW + WLAT(J)
      ENDDO
      XP1SH = XP1SH / SHW
      XP2SH = XP2SH / SHW
      XP3SH = XP3SH / SHW
      XPIB1SH = XPIB1SH / SHW
      XPIB2SH = XPIB2SH / SHW
      XPIB3SH = XPIB3SH / SHW
      XW1SH = XW1SH / SHW
      XW2SH = XW2SH / SHW
      XW3SH = XW3SH / SHW
      X100W1SH = X100W1SH / SHW
      X100W2SH = X100W2SH / SHW
      X100W3SH = X100W3SH / SHW
      IF (IDEBUG.GE.1) THEN
         PRINT *,'TOTAL NH LATITUDE WEIGHT = ', NHW
         PRINT *,'TOTAL SH LATITUDE WEIGHT = ', SHW
      ENDIF
!C
!C     CALCULATE THE CONSTANT TERMS
!C     R = RADIUS OF EARTH
!C     C = AXIAL MOMENT OF INERTIA (MOI) OF EARTH
!C     CMA = AXIAL MOI OF EARTH MINUS EQUATORIAL MOI OF EARTH
!C     G = GRAVITATIONAL CONSTANT
!C     O = MEAN ROTATION RATE OF EARTH
!C

      R = 6.37E+6
      C = 7.04E+37
      CMA = 0.00333*C
      G = 9.81
      O = 7.29E-5

!C     CW1=((-1.43*R**3)/(O*CMA*G))*1.0E7*2*PI
!C     CW2=CW1
!C     CW3=((R**3)/(C*O*G))*1.0E7*2*PI
!C     CP1=((-1.0*R**4)/(CMA*G))*1.0E7*2*PI
!C     CP2=CP1
!C     CP3=((0.7*R**4)/(C*G))*1.0E7*2*PI
      CW1=-0.1385226351E-03
      CW2=CW1
      CW3=0.3225736975E-06
      CP1=-0.4498334229E-01
      CP2=CP1
      CP3=0.1048561608E-03


!C
!C     MULTIPLY EACH TERM BY PROPER CONSTANT
!C
      XW1NH = XW1NH * CW1
      XW1SH = XW1SH * CW1
      XW2NH = XW2NH * CW2
      XW2SH = XW2SH * CW2
      XW3NH = XW3NH * CW3
      XW3SH = XW3SH * CW3
      X100W1NH = X100W1NH * CW1
      X100W1SH = X100W1SH * CW1
      X100W2NH = X100W2NH * CW2
      X100W2SH = X100W2SH * CW2
      X100W3NH = X100W3NH * CW3
      X100W3SH = X100W3SH * CW3
      XP1NH = XP1NH * CP1
      XP1SH = XP1SH * CP1
      XP2NH = XP2NH * CP2
      XP2SH = XP2SH * CP2
      XP3NH = XP3NH * CP3
      XP3SH = XP3SH * CP3
      XPIB1NH = XPIB1NH * CP1
      XPIB1SH = XPIB1SH * CP1
      XPIB2NH = XPIB2NH * CP2
      XPIB2SH = XPIB2SH * CP2
      XPIB3NH = XPIB3NH * CP3
      XPIB3SH = XPIB3SH * CP3

!C
!C     PRINT OUT RESULT
!C


!C==>  print out hemispheric values for analysis files
      WRITE(51,1000) X100W1NH, XW1NH, XP1NH,XPIB1NH,X100W1SH, XW1SH, XP1SH,XPIB1SH
      WRITE(51,1000) X100W2NH, XW2NH, XP2NH,XPIB2NH,X100W2SH, XW2SH, XP2SH,XPIB2SH
      WRITE(51,1000) X100W3NH, XW3NH, XP3NH,XPIB3NH,X100W3SH, XW3SH, XP3SH,XPIB3SH
 1000 FORMAT(8(1X,F9.5))

      RETURN
      END
