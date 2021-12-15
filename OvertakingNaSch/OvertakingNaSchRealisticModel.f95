! In this module we declare all constants used throughout the program and declare an object with Knospes' model properties
! If one wish to vary some of these constants (the deceleration probability, for example) you may comment the respective line in this module and declare it in the GLOBAL_VARIABLES module
MODULE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  TYPE :: PARTICLES
    INTEGER :: POSITION
    INTEGER :: VELOCITY
    INTEGER :: SPECIES
    INTEGER :: BRAKE_LIGHT
  END TYPE PARTICLES
  INTEGER, PARAMETER :: L = 500000
  INTEGER, PARAMETER :: MEASUREMENT_TIME = 100000
  INTEGER, PARAMETER :: TRANSIENT_TIME = 10000
  INTEGER, PARAMETER :: VEHICLE_LENGTH = 5
  INTEGER, PARAMETER :: MAXIMUM_VELOCITY = 20
  INTEGER, PARAMETER :: SECURITY_GAP = 7
  REAL, PARAMETER :: CUT_OFF_TIME = 6.0
  REAL, PARAMETER :: P_f = 0.1
  REAL, PARAMETER :: P_0 = 0.5
  REAL, PARAMETER :: P_B = 0.94
END MODULE EXTERNAL_PARAMETERS
! Here we declare all the global variables   
MODULE GLOBAL_VARIABLES
  USE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  TYPE(PARTICLES), DIMENSION(:), ALLOCATABLE :: VEHICLES
  INTEGER, DIMENSION(:), ALLOCATABLE :: STREET
  INTEGER :: ISEED, N
END MODULE GLOBAL_VARIABLES
! Functions used in the main program and routines
MODULE FUNCTIONS
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  IMPLICIT NONE
  CONTAINS
! This function generates random numbers in the interval [0:1) and random integer seeds in the interval [1:134456)
  FUNCTION RandomGenerator(Seed)
    INTEGER, INTENT(INOUT) :: Seed
    REAL :: RandomGenerator
    ISEED = mod(8121*seed+28411, 134456)
    RandomGenerator = real(ISEED)/134456.
    RETURN
  END FUNCTION RandomGenerator
END MODULE FUNCTIONS

MODULE ROUTINES
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  IMPLICIT NONE
  CONTAINS
! This routine initializes the vehicles' vector with random position and maximum velocity. The velocity is corrected afterwards 
  SUBROUTINE  RandomPositionInitialization(SizeOfSystem)
    INTEGER, INTENT(IN) :: SizeOfSystem
    INTEGER :: i, acum, aux
    LOGICAL, DIMENSION(0:INT(L/VEHICLE_LENGTH)) :: AuxiliaryVector
    STREET = 0
    AuxiliaryVector = .FALSE.
    acum = 1
    DO WHILE(acum <= SizeOfSystem)
      aux = INT((L/VEHICLE_LENGTH)*RandomGenerator(ISEED))
      IF(.NOT.AuxiliaryVector(aux))THEN
        AuxiliaryVector(aux) = .TRUE.
        VEHICLES(acum)%POSITION = VEHICLE_LENGTH*aux
        VEHICLES(acum)%VELOCITY = MAXIMUM_VELOCITY
        acum = acum + 1
      END IF
    END DO
    DO i = 0, INT(L/VEHICLE_LENGTH)
      IF(AuxiliaryVector(i))THEN
        DO aux = 0, VEHICLE_LENGTH
          STREET(mod(VEHICLE_LENGTH*i+aux,L)) = i
        END DO
      END IF
    END DO
  END SUBROUTINE RandomPositionInitialization
! This routine reorder the vehicles using the bubble algorithm. This is important when overtaking is allowed. It takes the extremities of the vehicles' vector as inputs. 
  SUBROUTINE ParticlesReordering(Beginnig, Ending)
    INTEGER, INTENT(IN) :: Beginnig, Ending
    LOGICAL :: flag
    INTEGER :: i
    TYPE(PARTICLES) :: AuxiliaryVector
    DO
      flag=.TRUE.
      DO i=Beginnig , Ending-1
        IF(VEHICLES(i)%POSITION>VEHICLES(i+1)%POSITION)THEN
          AuxiliaryVector=VEHICLES(i)
          VEHICLES(i)=VEHICLES(i+1)
          VEHICLES(i+1)=AuxiliaryVector
          flag=.FALSE.
        END IF
      END DO
      IF(flag) EXIT
    END DO
  END SUBROUTINE ParticlesReordering
!As the algorithm presents a direction to update and non-physical featurs may appear if the updating of the overtaking algorithm always starts from the same particle, then this routine chooses regions where the initialization will cause the lowest disturbances possible. It gives two output corresponding to the limits of the vehicles' vector   
  SUBROUTINE ChooseTheFirst(Beginning, Ending)
    INTEGER, INTENT(OUT) :: Ending, Beginning
    INTEGER :: i, gap, aux, aux1, aux2, aux3, aux4, aux5 ! the number of aux parameter must be equal to the VEHICLE_LENGTH parameter plus one... use a vector form INTEGER, DIMENSION(0:VEHICLE_LENGTH+1) :: aux, if you will
    !This piece of code decides if there is a piece of the system in free state where will come the first vehicle
    !If there are a vehicle which velocity is inferior to the space in front of it and the vehicle behind can not overtake, then it must behave as a cooperator and is chosen as the first one
    DO i = N, 1, -1
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - VEHICLE_LENGTH
      aux = VEHICLES(i)%POSITION - VEHICLES(i - 1)%POSITION + VEHICLES(i)%VELOCITY
      IF(VEHICLES(i)%VELOCITY < gap .AND. VEHICLES(i - 1)%VELOCITY <= aux )THEN
        Ending = i
        Beginning = i-N+1
        Return
      END IF
      if(VEHICLES(i)%SPECIES == 0)then
        Ending = i
        Beginning = i-N+1
        Return
      End if
    END DO
    !If the piece of code above fails, thats it, if the all the system is in congested state the first vehicle is chosen using the piece of code below
    !This code is based on the fact that overtaking is not very commom even if all can overtake to inicialize the update from a particle that will neither overtake or be overtaken anyway. If there are a vehicle too slow to overtake the vehicles in front as well as the vehicles behind it, then it must behave as a cooperator and is chosen as the first one.
    DO i = N, 1, -1
      aux1 = VEHICLES(i)%POSITION + VEHICLES(i)%VELOCITY + VEHICLE_LENGTH
      aux2 = VEHICLES(i - 1)%POSITION + VEHICLES(i - 1)%VELOCITY + VEHICLE_LENGTH
      aux3 = VEHICLES(i - 2)%POSITION + VEHICLES(i - 2)%VELOCITY + VEHICLE_LENGTH
      aux4 = VEHICLES(i - 3)%POSITION + VEHICLES(i - 3)%VELOCITY + VEHICLE_LENGTH
      aux5 = VEHICLES(i + 1)%POSITION + VEHICLES(i + 1)%VELOCITY + VEHICLE_LENGTH
      IF(aux1 <= aux5)THEN
        IF(aux2 <= aux5)THEN
          IF(aux3 <= aux5)THEN
            IF(aux4 <= aux5)THEN
              Ending = i
              Beginning = i-N+1
              RETURN
            END IF
          END IF
        END IF
      END IF
    END DO
    Ending = N
    Beginning = 1
  END SUBROUTINE ChooseTheFirst
! This routine inicializes a mixed population randomly
  SUBROUTINE RandomSpeciesInitialization(Specie0, Specie1)
    INTEGER, INTENT(IN) :: Specie0, Specie1
    INTEGER :: aux, acum
    If(Specie0==N)then
      VEHICLES%SPECIES = 0
      Return
    Else if(Specie1==N)then
      VEHICLES%SPECIES = 1
      Return
    Else if(Specie0>Specie1)then
      VEHICLES%SPECIES = 0
      acum = 1
      Do While(acum <= Specie1)
        aux = INT(N*RandomGenerator(ISEED))+1
        If(VEHICLES(aux)%SPECIES == 0)then
          VEHICLES(aux)%SPECIES = 1
          acum = acum + 1
        End if
      End do
    Else
      VEHICLES%SPECIES = 1
      acum = 1
      Do While(acum <= Specie0)
        aux = INT(N*RandomGenerator(ISEED))+1
        If(VEHICLES(aux)%SPECIES == 1)then
          VEHICLES(aux)%SPECIES = 0
          acum = acum + 1
        End if
      End do
    End if
  END SUBROUTINE RandomSpeciesInitialization
! This subroutine implements the main algorithm one particle at the time, it takes the number of the particle and its species as inputs and updates the overtaking number 
  SUBROUTINE MainAlgorithm(j, especie, OvertakingNumber)
    INTEGER, INTENT(IN) :: j, especie
    INTEGER, INTENT(INOUT) :: OvertakingNumber
    INTEGER :: gap, ExpectedVelocity, EffectiveGap, interactionradius, OldVelocity, p_status, Aux,&
    &OverDistance, temp
    REAL(KIND=16) :: TimeHeadway, p
    LOGICAL :: Overtaking

    OldVelocity = VEHICLES(j)%VELOCITY
    gap = VEHICLES(j+1)%POSITION - VEHICLES(j)%POSITION - VEHICLE_LENGTH
    ! Evaluation of the cut off time
    if(VEHICLES(j)%VELOCITY > CUT_OFF_TIME)THEN
      interactionradius = CUT_OFF_TIME
    else
      interactionradius = VEHICLES(j)%VELOCITY
    end if
    ! evaluation of the deceleration probability (step 0 in Knospe's work)
    if(VEHICLES(j)%VELOCITY == 0)then
      p = P_0
      p_status = 1
      TimeHeadway = 1.0*gap
    else
      TimeHeadway = gap/VEHICLES(j)%VELOCITY
      if((VEHICLES(j+1)%BRAKE_LIGHT==1).AND.TimeHeadway < interactionradius)then
        p = P_B
        VEHICLES(j)%BRAKE_LIGHT = 1
        p_status = 2
      else
        p = P_f
        p_status = 3
      end if
    end if
    ! determination of the effective gap and the space available for overtaking
    ExpectedVelocity = VEHICLES(j+2)%POSITION - VEHICLES(j+1)%POSITION - VEHICLE_LENGTH
    if(VEHICLES(j+1)%VELOCITY < ExpectedVelocity)then
      ExpectedVelocity = VEHICLES(j+1)%VELOCITY
    end if
    EffectiveGap = ExpectedVelocity - SECURITY_GAP
    if(EffectiveGap<0)then
      EffectiveGap = 0
    end if
    EffectiveGap = EffectiveGap + gap                                ! Effective gap
    OverDistance = EffectiveGap + ExpectedVelocity + VEHICLE_LENGTH  ! Minimum distance to cover to overtake
    !Evaluation of the acceleration (step 1 in Knospe's work)
    if((VEHICLES(j+1)%BRAKE_LIGHT==0).AND.(VEHICLES(j)%BRAKE_LIGHT==0) &
    &.OR.TimeHeadway >= interactionradius.OR.VEHICLES(j)%VELOCITY>=OverDistance)then
      if(VEHICLES(j)%VELOCITY < MAXIMUM_VELOCITY)then
        VEHICLES(j)%VELOCITY = VEHICLES(j)%VELOCITY + 1
      end if
    end if
    ! Evaluation of the behavior according to the specie (step 2 in Knospe's work)
    If(VEHICLES(j)%VELOCITY > EffectiveGap)then
      If (especie == 0)then
        VEHICLES(j)%VELOCITY = EffectiveGap
        VEHICLES(j)%BRAKE_LIGHT = 1
      Else if(especie == 1)then
        If(VEHICLES(j)%VELOCITY < OverDistance) then
          VEHICLES(j)%VELOCITY = EffectiveGap
          VEHICLES(j)%BRAKE_LIGHT = 1
        Else
          ! If the minimum distance to overtake is lower than the velocity of the overtaking vehicle then this piece of code will find wheter there is space available to complete the maneuver. If it fail, the overtaking does not happen
          temp = VEHICLES(j)%VELOCITY
          Do while (temp>EffectiveGap)
            Overtaking = .true.
            Do Aux = -VEHICLE_LENGTH,1
              If(.NOT.STREET(MOD(VEHICLES(j)%POSITION+temp+Aux,L))==0)then
                Overtaking = .false.
                Exit
              End if
            End Do
            if(Overtaking)then
              VEHICLES(j)%VELOCITY = temp
              exit
            else 
              temp = temp - VEHICLE_LENGTH - 1
            end if   
          End do
          If(Overtaking)then
            p = p_f
            p_status = 3
            OvertakingNumber = OvertakingNumber + 1
            VEHICLES(j)%BRAKE_LIGHT = 0
          Else
            VEHICLES(j)%VELOCITY = EffectiveGap
            VEHICLES(j)%BRAKE_LIGHT = 1
          End If
        End if
      End if
    End if
    ! Finaly, the randomization (step 3 in Knospe's work)
    If(VEHICLES(j)%VELOCITY > 0)then
      If(RandomGenerator(ISEED) < p)then
        VEHICLES(j)%VELOCITY = VEHICLES(j)%VELOCITY - 1
        If(p_status==2 .and. TimeHeadway < interactionradius)then
          VEHICLES(j)%BRAKE_LIGHT = 1
        End if
      End if
    End if
    ! This deviates form Knospe's algorithm as in the presented in their paper, states of braking light turned on are asymptotically absorbing in any density 
    If(VEHICLES(j)%VELOCITY > OldVelocity .or. VEHICLES(j)%VELOCITY == MAXIMUM_VELOCITY)then
      VEHICLES(j)%BRAKE_LIGHT = 0
    End if
    ! The road updated position of the particle is marked
    Do Aux = -VEHICLE_LENGTH+1, 0
      STREET(MOD(VEHICLES(j)%POSITION+VEHICLES(j)%VELOCITY+Aux,L)) = j + 2*N
    End Do

  END SUBROUTINE MainAlgorithm
! This routine uses the previous routine to initialize the system. It receives the number of defectors as a parameter
  SUBROUTINE InitialConditions(DefNumb)
    INTEGER, INTENT(IN) :: DefNumb
    INTEGER :: i, gap
    CALL RandomPositionInitialization(N)
    CALL ParticlesReordering(1, N)
    CALL RandomSpeciesInitialization(N-DefNumb,DefNumb)
    VEHICLES%POSITION = VEHICLES%POSITION + L + VEHICLE_LENGTH
    VEHICLES%BRAKE_LIGHT = 0
    VEHICLES(N+1) = VEHICLES(1)
    VEHICLES(N+1)%POSITION = VEHICLES(N+1)%POSITION + L
    VEHICLES(N+2) = VEHICLES(2)
    VEHICLES(N+2)%POSITION = VEHICLES(N+2)%POSITION + L
    Do i = N, 1, -1
      gap = VEHICLES(i+1)%POSITION - VEHICLES(i)%POSITION - VEHICLE_LENGTH
      If(gap < MAXIMUM_VELOCITY)then
        VEHICLES(i)%VELOCITY = gap
        VEHICLES(i)%BRAKE_LIGHT = 1
      Else
        VEHICLES(i)%VELOCITY = MAXIMUM_VELOCITY
      End if
    End do
    Do i = 0, -N+1, -1
      VEHICLES(i) = VEHICLES(i+N)
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION - L
    End Do
    VEHICLES(N+1) = VEHICLES(1)
    VEHICLES(N+1)%POSITION = VEHICLES(N+1)%POSITION + L
    VEHICLES(N+2) = VEHICLES(2)
    VEHICLES(N+2)%POSITION = VEHICLES(N+2)%POSITION + L
  END SUBROUTINE InitialConditions

! ##############################################################
! Now we begin the measurements
! ##############################################################

! This routine measures the fundamental diagram measuring the 
  SUBROUTINE FundamentalDiagram
    INTEGER :: i, j, configurations, OvertakingNumber, Transient, Measuring, Beginning, Ending
    REAL(KIND=16) :: velocity, localconcentration, flux
    LOGICAL :: Crossing
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A22)")"FundamentalDiagram.dat"
    OPEN(1,FILE=TRIM(Arq1))
    ALLOCATE(STREET(0:L-1))
    Transient = 10000
    Measuring = 100
    Do N = INT(0.01*(L/VEHICLE_LENGTH)), INT(0.99*(L/VEHICLE_LENGTH)),INT(0.01*(L/VEHICLE_LENGTH))
      ALLOCATE(VEHICLES(-N+1:N+2))
      Do configurations = 1,1
        Call InitialConditions(0)
        OvertakingNumber = 0
        Do i = 1,Transient
          STREET = 0
          Call ChooseTheFirst(Beginning, Ending)
          Do j = Ending, Beginning, -1
            Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
          End do
          Call ParticlesReordering(Beginning, Ending)
          Do j = Beginning-1, -N+1, -1
            VEHICLES(j) = VEHICLES(j+N)
            VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
          End Do
          Do j = Ending+1, N+2
            VEHICLES(j) = VEHICLES(j-N)
            VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
          End Do
        End do

        velocity = 0
        flux = 0
        OvertakingNumber = 0
        Do i = 1, Measuring
          STREET = 0
          Call ChooseTheFirst(Beginning, Ending)
          Do j = Ending, Beginning, -1
            Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
          End do
          Do j = Beginning, Ending
            Crossing = .false.
            if(MOD(VEHICLES(j)%POSITION,L) < L/2)then
              Crossing = .true.
            end if
            VEHICLES(j)%POSITION=VEHICLES(j)%VELOCITY+VEHICLES(j)%POSITION
            if(MOD(VEHICLES(j)%POSITION,L) >= L/2 .and. Crossing)then
              Crossing = .true.
              velocity = velocity + VEHICLES(j)%VELOCITY
              flux = flux + 1.
            end if
          End Do
          Call ParticlesReordering(Beginning, Ending)
          Do j = Beginning-1, -N+1, -1
            VEHICLES(j) = VEHICLES(j+N)
            VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
          End Do
          Do j = Ending+1, N+2
            VEHICLES(j) = VEHICLES(j-N)
            VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
          End Do
        End do
        if(flux>0) then
          velocity = 5.4*velocity/flux
        else
          velocity = 0.0
        end if
        flux = flux*Measuring/3600
        localconcentration = 0
        do j = 1, N
          localconcentration = localconcentration + VEHICLES(j)%BRAKE_LIGHT
        end do
        localconcentration = localconcentration/N
        if(velocity>0) then
          !Simulational unities
          !write(1,*) flux/velocity, flux, velocity, OvertakingNumber
          !km/h unities
          write(1,*) 133.33*flux/velocity, 720*flux, 5.4*velocity, OvertakingNumber
        end if
      End do
      DEALLOCATE(VEHICLES)
    End do
  END SUBROUTINE FundamentalDiagram
! This routine measures the distance headway (gaps histogram), the time headway (time gaps histogram) and optimum velocity curve (average velocity as a function of the gap)   
  SUBROUTINE DistanceHeadway_TimeHeadway
    INTEGER :: i, j, configurations, OvertakingNumber, Transient, Measuring, Beginning, Ending, gap, Aux, ConfigNumber
    REAL(KIND=16) :: Cont1, Cont2, Cont3, Cont4
    LOGICAL :: Crossing
    REAL(KIND=16), DIMENSION(:,:,:), ALLOCATABLE :: Histogram
    CHARACTER(LEN=100) :: Arq1, Arq2, Arq3
    WRITE(Arq1,"(A19)")"DistanceHeadway.dax"
    OPEN(1,FILE=TRIM(Arq1))
    WRITE(Arq2,"(A15)")"TimeHeadway.dax"
    OPEN(2,FILE=TRIM(Arq2))
    WRITE(Arq3,"(A19)")"OptimumVelocity.dax"
    OPEN(3,FILE=TRIM(Arq3))

    ALLOCATE(STREET(0:L-1))
    Transient = 10000
    Measuring = 10000
    ConfigNumber = 10
    ALLOCATE(Histogram(4,0:L-1,ConfigNumber))
    N = INT(0.1*(L/VEHICLE_LENGTH))
    Histogram = 0
    ALLOCATE(VEHICLES(-N+1:N+2))
    Do configurations = 1,ConfigNumber
      Call InitialConditions(N)
      OvertakingNumber = 0
      Do i = 1,Transient
        STREET = 0
        Call ChooseTheFirst(Beginning, Ending)
        Do j = Ending, Beginning, -1
          Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
        End do
        Call ParticlesReordering(Beginning, Ending)
        Do j = Beginning-1, -N+1, -1
          VEHICLES(j) = VEHICLES(j+N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
        End Do
        Do j = Ending+1, N+2
          VEHICLES(j) = VEHICLES(j-N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
        End Do
      End do
      Do i = 1,Measuring
        STREET = 0
        Call ChooseTheFirst(Beginning, Ending)
        Do j = Ending, Beginning, -1
          Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
        End do
        Do j = Ending, Beginning, -1
          Crossing = .false.
          if(MOD(VEHICLES(j)%POSITION,L) < L/2)then
            Crossing = .true.
          end if
          VEHICLES(j)%POSITION=VEHICLES(j)%VELOCITY+VEHICLES(j)%POSITION
          if(MOD(VEHICLES(j)%POSITION,L) >= L/2 .and. Crossing)then
            Crossing = .true.
            gap = VEHICLES(j+1)%POSITION-VEHICLES(j)%POSITION-VEHICLE_LENGTH
            Histogram(1,gap,configurations)=Histogram(1,gap,configurations)+1.
            Histogram(2,gap,configurations)=Histogram(2,gap,configurations)+VEHICLES(j)%VELOCITY
            If(VEHICLES(j)%VELOCITY>0)then
              If((VEHICLES(j)%VELOCITY>=MAXIMUM_VELOCITY-5).or.gap>=MAXIMUM_VELOCITY-5)then
                Aux = 10*gap/VEHICLES(j)%VELOCITY
                If(Aux<L)Then
                  Histogram(3,Aux,configurations)=Histogram(3,Aux,configurations)+1.
                End If
              Else   
                Aux = 10*gap/VEHICLES(j)%VELOCITY
                If(Aux<L)Then
                  Histogram(4,Aux,configurations)=Histogram(4,Aux,configurations)+1.
                End If
              End if  
            End if
          else
            Crossing = .false.
          end if
        End Do
        Call ParticlesReordering(Beginning, Ending)
        Do j = Beginning-1, -N+1, -1
          VEHICLES(j) = VEHICLES(j+N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
        End Do
        Do j = Ending+1, N+2
          VEHICLES(j) = VEHICLES(j-N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
        End Do
      End do
    End do
    Do j = 1, ConfigNumber
      Cont1 = 0.
      Cont2 = 0.
      Cont3 = 0.
      Cont4 = 0.
      Do i = 0,L-1
        Cont1 = Cont1 + Histogram(1,i,j)
        Cont3 = Cont3 + Histogram(3,i,j)
        Cont4 = Cont4 + Histogram(4,i,j)
        If(Histogram(1,i,j)>0)then
          Histogram(2,i,j) = Histogram(2,i,j)/Histogram(1,i,j)
        Else
          Histogram(2,i,j) = 0.
        End If
      End do
      Histogram(1,0:L-1,j) = Histogram(1,0:L-1,j)/Cont1
      Histogram(3,0:L-1,j) = Histogram(3,0:L-1,j)/Cont3
      Histogram(4,0:L-1,j) = Histogram(4,0:L-1,j)/Cont4
    End do
    
    Do i = 0,300!L-1
      Cont1 = 0.
      Cont2 = 0.
      Cont3 = 0.
      Cont4 = 0.
      Do j = 1, ConfigNumber
! print the results of each configuration        
        write(1,"(f9.2,A1)",ADVANCE='NO')Histogram(1,i,j)," "
        write(2,"(f9.2,A1,f9.2)",ADVANCE='NO')Histogram(3,i,j)," ", Histogram(4,i,j)
        write(3,"(f9.2,A1)",ADVANCE='NO')Histogram(2,i,j)," "
        Cont1 = Cont1 + Histogram(1,i,j)
        Cont2 = Cont2 + Histogram(2,i,j)
        Cont3 = Cont3 + Histogram(3,i,j)
        Cont4 = Cont4 + Histogram(4,i,j)
      End do
      Cont1 = Cont1/ConfigNumber
      Cont2 = Cont2/ConfigNumber
      Cont3 = Cont3/ConfigNumber
      Cont4 = Cont4/ConfigNumber
      write(1,"(f9.2,A1)",ADVANCE='NO')Cont1," "
      write(2,"(f9.2,A1,f9.2)",ADVANCE='NO')Cont3," ",Cont4
      write(3,"(f9.2,A1)",ADVANCE='NO')Cont2," "
      write(1,*)" "
      write(2,*)" "
      write(3,*)" "
    End Do
    Close(1)
    Close(2)
    Close(3)
    DEALLOCATE(VEHICLES)
    DEALLOCATE(STREET)
    DEALLOCATE(Histogram)
  END SUBROUTINE DistanceHeadway_TimeHeadway
! This routine measures the temporal Crosscovariances and temporal Autocovariances of the 3 variables of interest
  SUBROUTINE TemporalCrossCovariances
    INTEGER :: i, j, k, configurations, OvertakingNumber, Transient, Measuring, Beginning, Ending, ConfigNumber
    INTEGER :: Minute, tau
    REAL(KIND=16) :: Cont1, Cont2, Cont3, Cont4, Cont5, Cont6
    LOGICAL :: Crossing
    REAL(KIND=16), DIMENSION(:,:), ALLOCATABLE :: Variables 
    REAL(KIND=16), DIMENSION(9):: CrossProducts
    REAL(KIND=16), DIMENSION(6):: Covariances
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A20)")"CrossCovariances.dat"
    OPEN(1,FILE=TRIM(Arq1))
    
    Transient = 1000
    Measuring = 10000
    ConfigNumber = 10
    Minute = 60
    tau=100
    N = INT(0.1*(L/VEHICLE_LENGTH))
    ALLOCATE(STREET(0:L-1))
    ALLOCATE(Variables(3,Measuring))
    ALLOCATE(VEHICLES(-N+1:N+2))
    CrossProducts = 0
    Variables = 0
    Do configurations = 1,ConfigNumber
      Call InitialConditions(0)
      OvertakingNumber = 0
      Do i = 1,Transient
        STREET = 0
        Call ChooseTheFirst(Beginning, Ending)
        Do j = Ending, Beginning, -1
          Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
        End do
        Call ParticlesReordering(Beginning, Ending)
        Do j = Beginning-1, -N+1, -1
          VEHICLES(j) = VEHICLES(j+N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
        End Do
        Do j = Ending+1, N+2
          VEHICLES(j) = VEHICLES(j-N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
        End Do
      End do
      Do i = 1, Measuring
        Cont1 = 0.0
        Cont2 = 0.0
        Cont3 = 0.0
        Do k = 1, Minute 
          STREET = 0
          Call ChooseTheFirst(Beginning, Ending)
          Do j = Ending, Beginning, -1
            Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
          End do
          Do j = Ending, Beginning, -1
            Crossing = .false.
            if(MOD(VEHICLES(j)%POSITION,L) < L/2)then
              Crossing = .true.
            end if
            VEHICLES(j)%POSITION=VEHICLES(j)%VELOCITY+VEHICLES(j)%POSITION
            if(MOD(VEHICLES(j)%POSITION,L) >= L/2 .and. Crossing)then
              Crossing = .true.
              Cont1 = Cont1 + 1.
              Cont2 = Cont2 + VEHICLES(j)%VELOCITY
            end if
          End do
        End Do
        if(Cont1>0)then
          Cont2 = Cont2/Cont1
          Cont3 = Cont1/Cont2
        end if 
        Variables(1,i) = Variables(1,i) + Cont3  ! density
        Variables(2,i) = Variables(2,i) + Cont2  ! velocity
        Variables(3,i) = Variables(3,i) + Cont1  ! flux
        Call ParticlesReordering(Beginning, Ending)
        Do j = Beginning-1, -N+1, -1
          VEHICLES(j) = VEHICLES(j+N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
        End Do
        Do j = Ending+1, N+2
          VEHICLES(j) = VEHICLES(j-N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
        End Do
      End do
    End do
    Variables = Variables/ConfigNumber   ! average over configurations
    Do i = 0,tau
      CrossProducts = 0.0
      Cont1 = 0.0
      Cont2 = 0.0
      Cont3 = 0.0
      Cont4 = 0.0
      Cont5 = 0.0
      Cont6 = 0.0
      Do j = 1, Measuring - i
        Cont3 = Cont3 + Variables(1,j)  ! density
        Cont2 = Cont2 + Variables(2,j)  ! velocity
        Cont1 = Cont1 + Variables(3,j)  ! flux
        Cont4 = Cont4 + Variables(1,j)*Variables(1,j) 
        Cont5 = Cont5 + Variables(2,j)*Variables(2,j) 
        Cont6 = Cont6 + Variables(3,j)*Variables(3,j) 
        CrossProducts(1) = CrossProducts(1) + Variables(1,j)*Variables(1,j+i)  !density x density 
        CrossProducts(2) = CrossProducts(2) + Variables(2,j)*Variables(2,j+i)  !velociy x velocity
        CrossProducts(3) = CrossProducts(3) + Variables(3,j)*Variables(3,j+i)  !flux x flux
        CrossProducts(4) = CrossProducts(4) + Variables(1,j)*Variables(2,j+i)  !density x velocity
        CrossProducts(5) = CrossProducts(5) + Variables(1,j)*Variables(3,j+i)  !density x flux
        CrossProducts(6) = CrossProducts(6) + Variables(2,j)*Variables(1,j+i)  !velociy x density
        CrossProducts(7) = CrossProducts(7) + Variables(2,j)*Variables(3,j+i)  !velociy x flux
        CrossProducts(8) = CrossProducts(8) + Variables(3,j)*Variables(1,j+i)  !flux x density
        CrossProducts(9) = CrossProducts(9) + Variables(3,j)*Variables(2,j+i)  !flux x velocity
      End Do
      Cont1 = Cont1/(Measuring - i)
      Cont2 = Cont2/(Measuring - i)
      Cont3 = Cont3/(Measuring - i)
      Cont4 = Cont4/(Measuring - i)
      Cont5 = Cont5/(Measuring - i)
      Cont6 = Cont6/(Measuring - i)
      CrossProducts = CrossProducts/(Measuring - i)
      ! Autocovariances
      Covariances(1) = (CrossProducts(1) - Cont3*Cont3)/(Cont4 - Cont3*Cont3) !density x density  covariance
      Covariances(2) = (CrossProducts(2) - Cont2*Cont2)/(Cont5 - Cont2*Cont2) !velocity x velocity  covariance
      Covariances(3) = (CrossProducts(3) - Cont1*Cont1)/(Cont6 - Cont1*Cont1) !flux x flux  covariance
      ! Crosscovariances
      Covariances(4) = (CrossProducts(4) -  Cont2*Cont3)/&
      &SQRT((Cont4 - Cont3*Cont3)*(Cont5 - Cont2*Cont2))     !density x velocity  covariance
      Covariances(5) = (CrossProducts(5) -  Cont1*Cont3)/&
      &SQRT((Cont4 - Cont3*Cont3)*(Cont6 - Cont1*Cont1))     !density x flux  covariance
      Covariances(6) = (CrossProducts(7) -  Cont2*Cont1)/&
      &SQRT((Cont6 - Cont1*Cont1)*(Cont5 - Cont2*Cont2))     !velocity x flux  covariance
      write(1,*)i, Covariances(1), Covariances(2), Covariances(3), Covariances(4), Covariances(5), Covariances(6)
    End do
    Close(1)
    DEALLOCATE(VEHICLES)
    DEALLOCATE(STREET)
    DEALLOCATE(Variables)
  END SUBROUTINE TemporalCrossCovariances
! This routine measures the Pearson correlation of the velocities of two neighbohrs
  SUBROUTINE SpacialCorrelation
    INTEGER :: i, j, configurations, OvertakingNumber, Transient, Measuring, Beginning, Ending, ConfigNumber
    INTEGER :: FocalVehicle
    REAL(KIND=16) :: Correlation
    REAL(KIND=16), DIMENSION(2) :: FocalVariables
    REAL(KIND=16), DIMENSION(:,:), ALLOCATABLE :: Variables 
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A21)")"CrossCorrelations.dat"
    OPEN(1,FILE=TRIM(Arq1))
    
    Transient = 1000
    Measuring = 10000
    ConfigNumber = 10
    N = INT(0.1*(L/VEHICLE_LENGTH))
    FocalVehicle = N/2             ! vehicle to be measured
    ALLOCATE(STREET(0:L-1))
    ALLOCATE(Variables(3,0:N-1))
    ALLOCATE(VEHICLES(-N+1:N+2))
    Variables = 0.0
    FocalVariables = 0.0
    Do configurations = 1,ConfigNumber
      Call InitialConditions(0)
      OvertakingNumber = 0
      Do i = 1,Transient
        STREET = 0
        Call ChooseTheFirst(Beginning, Ending)
        Do j = Ending, Beginning, -1
          Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
        End do
        Do j = Ending, Beginning, -1
          VEHICLES(j)%POSITION=VEHICLES(j)%VELOCITY+VEHICLES(j)%POSITION
        End do
        Call ParticlesReordering(Beginning, Ending)
        Do j = Beginning-1, -N+1, -1
          VEHICLES(j) = VEHICLES(j+N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
        End Do
        Do j = Ending+1, N+2
          VEHICLES(j) = VEHICLES(j-N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
        End Do
      End do
      Do i = 1, Measuring
        STREET = 0
        Call ChooseTheFirst(Beginning, Ending)
        Do j = Ending, Beginning, -1
          Call MainAlgorithm(j, VEHICLES(j)%SPECIES, OvertakingNumber)
        End do
        Do j = Ending, Beginning, -1
          VEHICLES(j)%POSITION=VEHICLES(j)%VELOCITY+VEHICLES(j)%POSITION
        End do
        Call ParticlesReordering(Beginning, Ending)
        Do j = Beginning-1, -N+1, -1
          VEHICLES(j) = VEHICLES(j+N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION - L
        End Do
        Do j = Ending+1, N+2
          VEHICLES(j) = VEHICLES(j-N)
          VEHICLES(j)%POSITION = VEHICLES(j)%POSITION + L
        End Do
        FocalVariables(1) = FocalVariables(1) + VEHICLES(FocalVehicle)%VELOCITY                                 ! averege of the focal vehicle' velocity
        FocalVariables(2) = FocalVariables(2) + VEHICLES(FocalVehicle)%VELOCITY*VEHICLES(FocalVehicle)%VELOCITY ! averege of the focal vehicle' squared velocity
        Do j = FocalVehicle, FocalVehicle + N - 1
          Variables(1,j-FocalVehicle) = Variables(1,j-FocalVehicle) + VEHICLES(j)%VELOCITY   ! averege of the target vehicle' velocity                           
          Variables(2,j-FocalVehicle) = Variables(2,j-FocalVehicle) + VEHICLES(j)%VELOCITY*VEHICLES(j)%VELOCITY ! averege of the target vehicle' squared velocity                           
          Variables(3,j-FocalVehicle) = Variables(3,j-FocalVehicle) + VEHICLES(j)%VELOCITY*VEHICLES(FocalVehicle)%VELOCITY ! average of the cross product
        End do
      End Do
    End do
    FocalVariables = FocalVariables/(Measuring*ConfigNumber)
    Variables = Variables/(Measuring*ConfigNumber)
    Do i = 0,N
      Correlation = (Variables(3,i) - FocalVariables(1)*Variables(1,i))/&
      &(SQRT(Variables(2,i)-Variables(1,i)*Variables(1,i))*SQRT(FocalVariables(2)-FocalVariables(1)*FocalVariables(1)))
      write(1,*)i, Correlation
    End do
    
    DEALLOCATE(VEHICLES)
    DEALLOCATE(STREET)
    DEALLOCATE(Variables)
  END SUBROUTINE SpacialCorrelation
  
END MODULE ROUTINES


PROGRAM MAIN
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  USE ROUTINES
  IMPLICIT NONE
  ISEED = 3117
  ! uncoment to run one of the routine below
  
  !Call FundamentalDiagram
  !Call DistanceHeadway_TimeHeadway
  !Call TemporalCrossCovariances
  !Call SpacialCorrelation
  
END PROGRAM MAIN

