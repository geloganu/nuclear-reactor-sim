      program chicago_pile_1
            COMMON IX,ANUMBE
            real:: coord(3),vcoord(3)
            IX=99999

            open(1,file='lattice_output.out',status='UNKNOWN')
            !open(2,file='exit_energies.out',status='UNKNOWN')
            !open(3,file='unmod_energies.out',status='UNKNOWN')
            write(1,*) "backsc       transm       absorb        therma"

            !slab of thickness d
            rad=6.35!cm for UO2 bulk density=2g/cm3 and lump mass=2140g assuming perfect sphere
            spacing=21!cm for lattice spacing
            density=1.6
            ANUMBE=12.

            do 80 dist=1,25
            spacing=dist
            !initialize counters and position
            backsc=0.
            transm=0.
            absorb=0.
            therma=0.

            RMC=1.E4
            do 73 I=1,RMC

            do J=1,3
                  coord(J)=0.
            end do

            call ENERGY(ENE)
            !write(3,*)ENE
            call RANDM(RN)
            angle=ACOS(RN)
            call RANDM(RN)
            azim=2.*3.1415926*RN

            !random unit vector
            vcoord(1)=sin(angle)*cos(azim)
            vcoord(2)=sin(azim)*sin(angle)
            vcoord(3)=cos(angle)

34          call XSECT(ENE,XELAST,XABSOR)
            XTOTAL=XELAST+XABSOR
            call RANDM(RN)
            ell=-(log(RN)*ANUMBE)/(XTOTAL*density*6.02E23)

            coord=coord+vcoord*ell

            !impose periodic boundary conditions
            do J=1,3
                  if(coord(J).LT.-rad)coord(J)=2*spacing

                  if(coord(J).GT.rad)coord(J)=spacing

                  if(coord(J).GT.spacing+rad)coord(J)=0

                  if(coord(J).LT.spacing-rad)coord(J)=0
            end do
            !check backscatter
            !if(ABS(coord(1)).LT.rad.AND.ABS(coord(2)).LT.rad)then
            !      if(ABS(coord(3)).LT.rad)then
            !            backsc=backsc+1.
            !      goto 73
            !      end if
            !end if

            !loop through each coordinate
            do XJ=0,spacing,spacing
            do YJ=0,spacing,spacing
            do ZJ=0,spacing,spacing
                  if(ABS(coord(1)+XJ).LT.XJ+rad)then
                        if(ABS(coord(2)+YJ).LT.rad)then
                              if(ABS(coord(3)+ZJ).LT.rad)then
                                    transm=transm+1.
                                    if(ENE.LT.5e-4)therma=therma+1.
                                    goto 73
                              end if
                        end if
                  end if
                  
            end do
            end do
            end do
            
            !check if neutron is absorbed
            ratio=XABSOR/XTOTAL
            call RANDM(RN)
            if(RN.LT.ratio) then
                  absorb=absorb+1.
                  goto 73
            end if

            call ELOSS(ENE,angle,ENEAFT)
            ENE=ENEAFT

            call EULER(vcoord(1),vcoord(2),vcoord(3),angle,Sx,Sy,Sz)
            vcoord(1)=Sx
            vcoord(2)=Sy
            vcoord(3)=Sz
            goto 34

73          continue
            !print*,backsc/RMC,transm/RMC,absorb/RMC,therma/RMC
            write(1,*) backsc/RMC,transm/RMC,absorb/RMC,therma/RMC
80          end do

      contains
C ************************************************
      SUBROUTINE ELOSS(ENE,ANGL,ENEAFT)
      COMMON IX,ANUMBE
C *** CALCULATE ENERGY LOSS IN COLLISION, USE ISOTROPIC (S-WAVE) 
C *** APPROXIMATION. ROUND UP ANY ENERGY BELOW 0.1 EV TO 0.02 EV (THERMAL)
      IF(ENE.LT.1.E-7)THEN
      ENEAFT=2.E-8
      ELSE
C *** THIS IS THE MAXIMUM RECOIL ENERGY THAT CAN BE IMPARTED
      EREMAX=((4.*ANUMBE)/((1.+ANUMBE)**2.))*ENE
      CALL RANDM(RN1)
C *** THE PROBABILITY DISTRIBUTION FOR RECOIL ENERGIES IS FLAT
C *** FROM ZERO UP TO THE MAXIMUM 
      EREMAX=EREMAX*RN1
      ENEAFT=ENE-EREMAX
      ENDIF

      IF(ENE.LT.1.E-7)THEN
      CALL RANDM(RN1)
      ANGL=ACOS((2.*RN1)-1.)
      ELSE
      ETA=ACOS(SQRT((EREMAX/ENE)*(((1.+ANUMBE)**2.)/(4.*ANUMBE))))
      ANGL=ATAN((SIN(2.*ETA))/((1./ANUMBE)-COS(2.*ETA)))
      ENDIF
      
      RETURN
      END
C ************************************************
      SUBROUTINE ENERGY(ENE)
      COMMON IX
C *** USE VON NEUMANN'S METHOD ("HIT AND MISS") TO CALCULATE THE 
C *** ENERGY OF THE NEUTRON EMITED FROM A FISSION SOURCE.
   7  CALL RANDM(RN1)
C *** P=0 ABOVE ~10 MEV
      RN1=RN1*9.999
C *** APPLY MAXWELL SPECTRUM
      PROB=SQRT(RN1)*EXP(-RN1/1.4)
      CALL RANDM(RN2)
C *** THE MAXIMUM OF THE DISTRIBUTION AS WRITTEN IS LESS THAN 0.5 
      RN2=RN2*0.5
      IF(RN2.GT.PROB)GOTO 7    
      ENE=RN1   
      RETURN
      END
C ************************************************
      SUBROUTINE RANDM(RN1)
      COMMON IX
      IY=IX*6539
      IF(IY)5,6,6
    5 IY=IY+2147483647+1
    6 RN1=IY
      RN1=RN1*.4656613E-9
      IX=IY
      RETURN
      END
C ************************************************
      SUBROUTINE XSECT(ENE,XELAST,XABSOR)
C *** RETURNS THE APPROXIMATE CROSS SECTIONS (BARNS) FOR ABSORPTION AND 
C *** ELASTIC SCATTERING IN CARBON.
C *** HERE WE NEED THE ENERGY IN EV, NOT MEV
      X=ENE*1.E6

C *** CHECK THESE APPROXIMATIONS VS THE PLOTS IN THE HANDOUTS
      IF(X.LT.1.E4)THEN
      XELAST=5.
      ELSE
      XELAST=10.5-(1.346*LOG10(X))
      ENDIF
      XELAST=XELAST*1.E-24

      IF(X.LT.1.E3)THEN
      XABSOR=(6.442E-4)*(X**(-0.49421))
      ELSE IF(X.LT.1.E5) THEN
      XABSOR=1.5E-5
      ELSE IF(X.LT.5.E6) THEN
      XABSOR=2.E-5
      ELSE
      XABSOR=(4.E-06)*EXP(X*3.2189E-07) 
      ENDIF
      XABSOR=XABSOR*1.E-24
      RETURN
      END
C*************************************************
      SUBROUTINE EULER(EX,EY,EZ,ANGL,SX,SY,SZ)
!     THIS SUBROUTINE TAKES THE ORIGINAL LINEAR TRAJECTORY,
!     ROTATES IT TO LIE ALONG THE Z-AXIS, GENERATES A VECTOR
!     AT ZENITH ANGLE THETA = SCATTERING ANGLE 
!     AND AZIMUTHAL ANGLE FI = RANDOM * 2PI. THE
!     ORIGINAL AXIS IS NOW ROTATED BACK TAKING THE SCATTERING
!     VECTOR WITH IT. NOW WE HAVE THE SCATTERED
!     DIRECTION VECTOR (SX,SY,SZ).
!     WE USE EULER ANGLES TO PERFORM THE TRANSFORMATION.
      COMMON IX
!     NORMALIZE THE DIRECTION TO A UNIT VECTOR (IN CASE IT WASN'T)

      S=SQRT(EX**2+EY**2+EZ**2)
      EX=EX/S
      EY=EY/S
      EZ=EZ/S
      BET=ACOS(EZ)

!     THESE APPROXIMATIONS ARE ONLY NEEDED FOR COMPTON SCATTERING 
!     FOR GAMMAS (BUT THEY WILL NOT HURT HERE)

      IF(ABS(BET).LT.0.027)ALF=0.0
      IF(ABS(BET).LT.0.027)GO TO 44
      ARG=EY/SIN(BET)
      AARG=ABS(ARG)
      IF(AARG.LT.1.0)GOTO 344
      ARG=ARG/(1.0001*AARG)
344   ALF=ASIN(ARG)
 44   CONTINUE
      SCO1=COS(ALF)*SIN(BET)+EX
      SCO1=ABS(SCO1)
      SCO2=ABS(EX)
      IF(SCO1.LT.SCO2)BET=-BET
      IF(SCO1.LT.SCO2)ALF=-ALF
      GAM=0.0
!     WE NOW HAVE THE EULER ANGLES OF ROTATION BETWEEN
!     Z-AXIS TO DIRECTION OF INITIAL PARTICLE.
      THET = ANGL
      CALL RANDM(RN1)
      FI = 6.2831853 * RN1
!     WE NOW HAVE SCATTERED THE PARTICLE FROM THE
!     Z-AXIS AND MUST ROTATE IT TO THE ORIGINAL UNSCATTERED
!     PARTICLE DIRECTION. CACULATE THE ROTATION MATRIX.
      R11 = COS(ALF)*COS(BET)*COS(GAM)-SIN(ALF)*SIN(GAM)
      R12 = COS(BET)*SIN(ALF)*COS(GAM)+COS(ALF)*SIN(GAM)
      R13 =-SIN(BET)*COS(GAM)
      R21 =-SIN(GAM)*COS(BET)*COS(ALF)-SIN(ALF)*COS(GAM)
      R22 =-SIN(GAM)*COS(BET)*SIN(ALF)+COS(ALF)*COS(GAM)
      R23 = SIN(BET)*SIN(GAM)
      R31 = SIN(BET)*COS(ALF)
      R32 = SIN(ALF)*SIN(BET)
      R33 = COS(BET)
      SOX = SIN(THET)*COS(FI)
      SOY = SIN(THET)*SIN(FI)
      SOZ = COS(THET)
      SX =  R11*SOX+R21*SOY+R31*SOZ
      SY =  R12*SOX+R22*SOY+R32*SOZ
      SZ =  R13*SOX+R23*SOY+R33*SOZ
!     WE NOW HAVE THE UNIT PROPAGATION VECTOR OF THE
!     SCATTERED PARTICLE IN THE *ORIGINAL* FRAME.
      RETURN
      END
C ************************************************
      end program