! In this module we declare all constants used throughout the program as well as where we define abstract objects
! If one wish to vary some of these constants (the deceleration probability, for example) you may comment the respective line in this module and declare it in the GLOBAL_VARIABLES module  
MODULE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  TYPE :: PARTICLES 
    INTEGER :: POSITION
    INTEGER :: VELOCITY
    INTEGER :: SPECIES
    LOGICAL :: BRAKE_LIGHT
    REAL    :: DECELERATION_PROBABILITY 
  END TYPE PARTICLES
  INTEGER, PARAMETER :: rdp=selected_real_kind(10,300) 
  INTEGER, PARAMETER :: L = 50000                                                                                                                                                                                                                                                                                             
  INTEGER, PARAMETER :: MEASUREMENT_TIME = 100000
  INTEGER, PARAMETER :: TRANSIENT_TIME = 10000
  INTEGER, PARAMETER :: VEHICLE_LENGTH = 3
  INTEGER, PARAMETER :: MAXIMUM_VELOCITY = 20 
  INTEGER, PARAMETER :: SECURITY_GAP = 7 
  REAL, PARAMETER :: CUT_OFF_TIME = 4.0
  REAL, PARAMETER :: P = 0.1
  REAL, PARAMETER :: P_0 = .5
  REAL, PARAMETER :: P_B = .94
END MODULE EXTERNAL_PARAMETERS
! Here we declare all the global variables 
MODULE GLOBAL_VARIABLES
  USE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  TYPE(PARTICLES), DIMENSION(:), ALLOCATABLE :: VEHICLES
  LOGICAL, DIMENSION(0:L-1) :: STREET
  INTEGER :: ISEED, N, NUMBER_OF_DEFECTORS
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
  
MODULE SECONDARY_ROUTINES
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  IMPLICIT NONE
  CONTAINS
  
  ! This routine selects the position of a set (of size SizeOfSystem) of vehicles randomly without superposition 
  SUBROUTINE  RandomPositionInitialization(SizeOfSystem)
    INTEGER, INTENT(IN) :: SizeOfSystem
    INTEGER :: i, acum, aux
    LOGICAL, DIMENSION(0:INT(L/VEHICLE_LENGTH)) :: vetorauxiliar
    STREET = .FALSE.
    vetorauxiliar = .FALSE.
    acum = 1
    DO WHILE(acum <= SizeOfSystem)
      aux = INT((L/VEHICLE_LENGTH)*RandomGenerator(ISEED))      
      IF(.NOT.vetorauxiliar(aux))THEN 
        vetorauxiliar(aux) = .TRUE.   
        VEHICLES(acum)%POSITION = VEHICLE_LENGTH*aux
        VEHICLES(acum)%VELOCITY = MAXIMUM_VELOCITY     
        acum = acum + 1                       
      END IF
    END DO 
    DO i = 0, INT(L/VEHICLE_LENGTH)
      IF(vetorauxiliar(i))THEN
        DO aux = 0, VEHICLE_LENGTH
          STREET(VEHICLE_LENGTH*i+aux) = .TRUE.
        END DO
      END IF 
    END DO
  END SUBROUTINE RandomPositionInitialization
  !This subroutine organize the particles in positional order using the bubble algorithm
  !It requares the initial position and the final position of the piece of the VEHICLES object which one wish to order. 
  SUBROUTINE ParticlesReordering(Beginnig, Ending)
    INTEGER, INTENT(IN) :: Beginnig, Ending
    LOGICAL :: flag
    INTEGER :: i, aux
    TYPE(PARTICLES) :: auxiliaryvector      
    DO               
      flag=.TRUE.               
      DO i=Beginnig , Ending-1
        IF(VEHICLES(i)%POSITION>VEHICLES(i+1)%POSITION)THEN
          auxiliaryvector=VEHICLES(i)      
          VEHICLES(i)=VEHICLES(i+1)
          VEHICLES(i+1)=auxiliaryvector
          flag=.FALSE.
        END IF 
      END DO 
      IF(flag) EXIT 
    END DO
  END SUBROUTINE ParticlesReordering
  !This routine selects vehicles randomly and initialize then as belonging to a given species
  !It requires the number of particles and the the number of defectors
   SUBROUTINE RandomSpeciesInitialization(NumberOfParticles, NumberOfdefectors)
    INTEGER, INTENT(IN) :: NumberOfParticles, NumberOfdefectors
    INTEGER :: aux, acum, cooperators, defectors
    defectors = NumberOfdefectors 
    cooperators = NumberOfParticles - NumberOfdefectors
    IF(cooperators <= defectors)THEN  
      VEHICLES%SPECIES = 1
      acum = 1                            
      DO WHILE(acum <= cooperators)
        aux = INT(NumberOfParticles*RandomGenerator(ISEED))
        IF (aux > 0)THEN
          IF(VEHICLES(aux)%SPECIES == 1)THEN 
            VEHICLES(aux)%SPECIES = 0
            acum = acum + 1                       
          END IF 
        END IF 
      END DO 
    ELSE IF(cooperators > defectors)THEN  
      VEHICLES%SPECIES = 0
      acum = 1                            
      DO WHILE(acum <= defectors)
        aux = INT(NumberOfParticles*RandomGenerator(ISEED))
        IF (aux > 0)THEN
          IF (VEHICLES(aux)%SPECIES == 0)THEN 
            VEHICLES(aux)%SPECIES = 1
            acum = acum + 1                      
          END IF 
        END IF 
      END DO 
    END IF 
  END SUBROUTINE RandomSpeciesInitialization
  !This routine chooses the first vehicle to initialize the velocity update
  SUBROUTINE ChooseTheFirst(Ending)
    INTEGER, INTENT(OUT) :: Ending
    INTEGER :: i, gap, aux, aux1, aux2, aux3, aux4 ! the number of aux parameter must be equal to the VEHICLE_LENGTH parameter ... use a vector form INTEGER, DIMENSION(0:VEHICLE_LENGTH) :: aux, if you will 
    LOGICAL :: flag
    flag = .FALSE.
    Ending = N
    !This piece of code decides the first vehicle if there is a piece of the system in free state
    !If there are a vehicle which velocity is inferior to the space in front of it and the vehicle behind can not overtake, then it must behave as a cooperator and is chosen as the first one  
    DO i = N, 1, -1 
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - VEHICLE_LENGTH
      aux = VEHICLES(i)%POSITION - VEHICLES(i - 1)%POSITION + VEHICLES(i)%VELOCITY
      IF(VEHICLES(i)%VELOCITY < gap .AND. VEHICLES(i - 1)%VELOCITY <= aux )THEN
        Ending = i
        flag = .TRUE.
        EXIT
      END IF 
    END DO     
    !If the piece of code above fails, thats it, if the all the system is in congested state the first vehicle is chosen using the piece of code below
    !If there are a vehicle too slow to overtake the vehicles in front and the vehicles as well as the vehicles behind it, then it must behave as a cooperator and is chosen as the first one  
    IF(.NOT.flag)THEN
      DO i = N, 1, -1
        aux1 = VEHICLES(i)%POSITION + VEHICLES(i)%VELOCITY + VEHICLE_LENGTH 
        aux2 = VEHICLES(i - 1)%POSITION + VEHICLES(i - 1)%VELOCITY + VEHICLE_LENGTH
        aux3 = VEHICLES(i - 2)%POSITION + VEHICLES(i - 2)%VELOCITY + VEHICLE_LENGTH
        aux4 = VEHICLES(i - 3)%POSITION + VEHICLES(i - 3)%VELOCITY + VEHICLE_LENGTH
        IF(aux1 <= VEHICLES(i + 1)%POSITION)THEN
          IF(aux2 <= VEHICLES(i + 1)%POSITION)THEN
            IF(aux3 <= VEHICLES(i + 1)%POSITION)THEN
              IF(aux4 <= VEHICLES(i + 1)%POSITION)THEN
                Ending = i
                flag = .TRUE.
                EXIT
              END IF 
            END IF  
          END IF 
        END IF 
      END DO 
    END IF
  END SUBROUTINE ChooseTheFirst
  !This routine performs the second and tird steps of the  algorithm of the cooperators
  SUBROUTINE CooperatorsAlgorithm(i, gap, gapII, initialvelocity)
    INTEGER, INTENT(IN) :: i, gap, gapII, initialvelocity
    INTEGER :: aux, efectivegap
    ! aux=MIN(v_(n+1),gap(n+1))
    IF(VEHICLES(i + 1)%VELOCITY < gapII)THEN
      aux = VEHICLES(i + 1)%VELOCITY
    ELSE 
      aux = gapII
    END IF  
    ! aux=MAX(aux-SECURITY_GAP, 0)
    aux = aux - SECURITY_GAP
    IF(aux < 0)THEN
      aux = 0
    END IF  
    ! EfectiveGap= gap + aux = gap + MAX(MIN(v_(n+1),gap(n+1))-SECURITY_GAP, 0)
    efectivegap = gap + aux
    ! v_(n)(t+1)=MIN(EfectiveGap,v_(n)(t))
    IF(VEHICLES(i)%VELOCITY > efectivegap)THEN
      VEHICLES(i)%VELOCITY = efectivegap
    END IF   
    ! If the vehicle decelerates because of other vehicles its brake lights turn on
    IF(VEHICLES(i)%VELOCITY < initialvelocity)THEN
      VEHICLES(i)%BRAKE_LIGHT = .TRUE.
    ELSE
      VEHICLES(i)%BRAKE_LIGHT = .FALSE.
    END IF   
    ! NaSch tird step
    IF(VEHICLES(i)%VELOCITY > 0)THEN
      IF(RandomGenerator(ISEED) < VEHICLES(i)%DECELERATION_PROBABILITY)THEN
        VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1
      END IF   
    END IF
    !finaly we mark the position and trajectory of the vehicle
    DO aux = 0, VEHICLES(i)%VELOCITY
      STREET(MOD(VEHICLES(i)%POSITION+aux,L))= .TRUE.
    END DO          
    DO aux = 0, VEHICLE_LENGTH
      STREET(MOD(VEHICLES(i)%POSITION+VEHICLES(i)%VELOCITY+aux,L))= .TRUE.
    END DO          
  END SUBROUTINE CooperatorsAlgorithm
  !This routine performs the second and tird steps of the  algorithm of the defectors
  SUBROUTINE DefectorsAlgorithm(i, futuregap, gap, gapII, initialvelocity)
    INTEGER, INTENT(IN) :: i, futuregap, gap, gapII, initialvelocity
    INTEGER :: auxB, auxE, aux
    !The vehicle is a defector and is supposed to have enough velocity to overtake at least the front vehicle
    !We begin with the random deceleration 
    IF(RandomGenerator(ISEED) < VEHICLES(i)%DECELERATION_PROBABILITY)THEN
        VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1
    END IF 
    ! If the vehicle still cannot overtake then the CooperatorsAlgorithm routine is called, else we keep testing 
    IF(VEHICLES(i)%VELOCITY>futuregap)THEN
      DO
        !If the vehicle fails the test we call the CooperatorsAlgorithm routine. The overtaking is frustrated
        IF(VEHICLES(i)%VELOCITY <= gap + VEHICLE_LENGTH)THEN
          CALL CooperatorsAlgorithm(i, gap, gapII, initialvelocity)
          RETURN
        END IF   
        auxB=MOD(VEHICLES(i)%POSITION+VEHICLES(i)%VELOCITY,L)               !rear o the proposed posiion of the vehicle 
        auxE=MOD(VEHICLES(i)%POSITION+VEHICLES(i)%VELOCITY+VEHICLE_LENGTH,L)!front of the proposed position of the vehicle   
        IF(.NOT.STREET(auxB))THEN  
          IF(.NOT.STREET(auxE))THEN            !If both spots are free, the vehicle have space to make the maneuver
            VEHICLES(i+1)%BRAKE_LIGHT = .TRUE. !Turning the brake light of the overtaken vehicle on 
            VEHICLES(i)%BRAKE_LIGHT = .FALSE.  !Compactible with overtaking
            EXIT                               !We exit the tests
          END IF   
        END IF 
        VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1  !If the overtaking of the j-th vehicle fails, decrease the velocity and try again until either, find a suitble empty space, or fail the overtaking
      END DO
      !finaly we mark the position and trajectory of the vehicle
      DO aux = 0, VEHICLES(i)%VELOCITY
        STREET(MOD(VEHICLES(i)%POSITION+aux,L))= .TRUE.
      END DO          
      DO aux = 0, VEHICLE_LENGTH
        STREET(MOD(VEHICLES(i)%POSITION+VEHICLES(i)%VELOCITY+aux,L))= .TRUE.
      END DO          
    ELSE   
      CALL CooperatorsAlgorithm(i, gap, gapII, initialvelocity)
      RETURN
    END IF  
  END SUBROUTINE DefectorsAlgorithm
  
END MODULE SECONDARY_ROUTINES

MODULE PRIMARY_ROUTINES 
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  USE SECONDARY_ROUTINES
  IMPLICIT NONE
  CONTAINS
  !This routine initializes the system using the routineS of THE SECUNDARY_ROUTINES module
  SUBROUTINE InitialConditions
    INTEGER :: i, gap
    CALL RandomPositionInitialization(N)                      !Position initialization
    CALL RandomSpeciesInitialization(N, NUMBER_OF_DEFECTORS)  !Species initialization
    CALL ParticlesReordering(1, N)                            !Ordering of the vehicles 
    VEHICLES%BRAKE_LIGHT = .FALSE.
    !initialization of the velocities
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L 
    DO i = N, 1, -1
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - VEHICLE_LENGTH
      IF(gap < MAXIMUM_VELOCITY)THEN
        VEHICLES(i)%VELOCITY = gap
        VEHICLES(i)%BRAKE_LIGHT = .TRUE.
      ELSE
        VEHICLES(i)%VELOCITY = MAXIMUM_VELOCITY
      END IF   
    END DO      
    ! Aplication of the periodic boundary conditions TO THE VEHICLES. That it, initialization of the VEHICLE vector with 2N+2 entries
    VEHICLES%POSITION = VEHICLES%POSITION + L
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L
    VEHICLES(N + 2) = VEHICLES(2)
    VEHICLES(N + 2)%POSITION = VEHICLES(N + 2)%POSITION + L 
    DO i = -N, 0
     VEHICLES(i)%VELOCITY = VEHICLES(i + N)%VELOCITY 
     VEHICLES(i)%SPECIES= VEHICLES(i + N)%SPECIES
     VEHICLES(i)%BRAKE_LIGHT = VEHICLES(i + N)%BRAKE_LIGHT
     VEHICLES(i)%DECELERATION_PROBABILITY = VEHICLES(i + N)%DECELERATION_PROBABILITY
     VEHICLES(i)%POSITION = VEHICLES(i + N)%POSITION - L 
    END DO 
    VEHICLES%POSITION = VEHICLES%POSITION + L
  END SUBROUTINE InitialConditions
  !This routine perform the aceleration step, as well as the decision of the deceleration factor adopted by the vehicle
  SUBROUTINE AcelerationProcess
    INTEGER :: i, gap, aux
    REAL :: timeheadway, interactionradius
    !We perform the step to the 1 to N vehicles and copy the result to the rest
    DO i = N, 1, -1
      !All v=0 vehicles have the same velocity update rule than it its written here for computing time reasons
      IF(VEHICLES(i)%VELOCITY == 0)THEN
        VEHICLES(i)%DECELERATION_PROBABILITY = P_0
        VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY + 1
        CYCLE  
      END IF     
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - VEHICLE_LENGTH
      SELECT CASE (VEHICLES(i)%SPECIES) ! Diferent species have diferent update rules  
      CASE(0)                           ! Cooperators rules                         
        !InteractionRadius = MIN(v_n, CUT_OFF_TIME)
        IF(VEHICLES(i)%VELOCITY >= CUT_OFF_TIME)THEN
          interactionradius = CUT_OFF_TIME
        ELSE
          interactionradius = VEHICLES(i)%VELOCITY
        END IF 
        timeheadway = gap/VEHICLES(i)%VELOCITY       !t_h(n)=d_n/v_n  notice that v_n .ne. 0
        !if ((b_(n+1) = 0) and (b_n = 0)) or (t_h .ge. t_s ) then: v_n (t + 1) = min(v_n (t) + 1, v_max ).
        IF(VEHICLES(i + 1)%BRAKE_LIGHT .AND. timeheadway < interactionradius)THEN
          VEHICLES(i)%DECELERATION_PROBABILITY = P_B  !vehicle interacting does not have an increase in its velocty
          CYCLE
        ELSE 
          VEHICLES(i)%DECELERATION_PROBABILITY = P    ! free vehicle that have an increase in its velocty
          IF((.NOT.VEHICLES(i + 1)%BRAKE_LIGHT .AND. .NOT.VEHICLES(i)%BRAKE_LIGHT) &
          &.OR. timeheadway >= interactionradius)THEN
            IF(VEHICLES(i)%VELOCITY < MAXIMUM_VELOCITY)THEN
              VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY + 1
            END IF   
          END IF 
        END IF 
      CASE(1)                               !Defectors rule
        aux = gap + VEHICLES(i + 1)%VELOCITY
        !Does the vehicle have enough velocity to overtake? behave as a free vehicle!
        IF(VEHICLES(i)%VELOCITY > aux)THEN  
          VEHICLES(i)%DECELERATION_PROBABILITY = P
          IF(VEHICLES(i)%VELOCITY < MAXIMUM_VELOCITY)THEN 
            VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY + 1
          END IF   
          CYCLE
        !Else, behave as a cooperator   
        ELSE 
          IF(VEHICLES(i)%VELOCITY >= CUT_OFF_TIME)THEN
            interactionradius = CUT_OFF_TIME
          ELSE
            interactionradius = VEHICLES(i)%VELOCITY
          END IF 
          timeheadway = gap/VEHICLES(i)%VELOCITY
          IF(VEHICLES(i + 1)%BRAKE_LIGHT .AND. timeheadway < interactionradius )THEN
            VEHICLES(i)%DECELERATION_PROBABILITY = P_B
            CYCLE
          ELSE 
            VEHICLES(i)%DECELERATION_PROBABILITY=P
            IF((.NOT.VEHICLES(i + 1)%BRAKE_LIGHT .AND. .NOT.VEHICLES(i)%BRAKE_LIGHT)&
            &.OR. timeheadway >= interactionradius)THEN
              IF(VEHICLES(i)%VELOCITY<MAXIMUM_VELOCITY)THEN 
                VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY + 1
              END IF   
            END IF   
          END IF 
        END IF  
      END SELECT  
    END DO
    DO i = N + 1, N + 2
      VEHICLES(i) = VEHICLES(i - N)
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION + L 
    END DO 
    DO i = 0, -N, -1
     VEHICLES(i)=VEHICLES(i+N)
     VEHICLES(i)%POSITION=VEHICLES(i)%POSITION-L 
    END DO    
  END SUBROUTINE AcelerationProcess
  !This routine perform the response of the vehicles to interactions among thenselves sellecting the vehicles who perform a cooperator-like behavior or a defector-like behavior. 
  SUBROUTINE InteractionProcess    
    INTEGER ::i, j, ending, beginning, gap, gapII, aux, initialvelocity, futuregap
    !This is important for security reasons and is important only for the defectors
    STREET = .FALSE.    
    DO i = 1, N 
      DO j = 0, VEHICLE_LENGTH
        STREET(MOD(VEHICLES(i)%POSITION + j, L)) = .TRUE.
      END DO  
    END DO 
    !We choose the vehicle to begin the update of the velocities
    CALL ChooseTheFirst(ending)
    beginning = ending - N + 1
    DO i = ending, beginning, -1     
      initialvelocity = VEHICLES(i)%VELOCITY
      futuregap = VEHICLES(i + 1)%POSITION + VEHICLES(i + 1)%VELOCITY - VEHICLES(i)%POSITION + VEHICLE_LENGTH 
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - VEHICLE_LENGTH 
      gapII = VEHICLES(i + 2)%POSITION - VEHICLES(i + 1)%POSITION - VEHICLE_LENGTH 
      IF(VEHICLES(i)%VELOCITY <= gap + VEHICLE_LENGTH)THEN  !Overtaking is impossible: call the cooperators routine regardless of its strategy
        CALL CooperatorsAlgorithm(i, gap, gapII, initialvelocity)
        CYCLE
      END IF   
      IF(VEHICLES(i)%SPECIES== 0)THEN   !The vehicle is a cooperator: call the cooperators routine
        CALL CooperatorsAlgorithm(i, gap, gapII, initialvelocity)
        CYCLE
      ELSE                              !The vehicle is a defector and may overtake: call the defectors routine
        CALL DefectorsAlgorithm(i, futuregap, gap, gapII, initialvelocity)
        CYCLE
      END IF   
    END DO    
    !Update of the positions
    DO i = beginning, ending 
      VEHICLES(i)%POSITION=VEHICLES(i)%VELOCITY+VEHICLES(i)%POSITION
    END DO
    !Reordering of the vehicles 
    CALL ParticlesReordering(beginning, ending)
    !Updating all other vehicles 
    DO i = ending + 1, N + 2
      VEHICLES(i) = VEHICLES(i - N)
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION + L 
    END DO 
    DO i = beginning - 1, -N , -1
      VEHICLES(i) = VEHICLES(i + N)
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION - L 
    END DO 
  END SUBROUTINE InteractionProcess
      
  ! We initialize the system defining a global density c. 
  ! Then we divide the highway in 10 segments and measure the concentration and the flux in each segment, printing to the file  FundamentalDiagram.dat.
  ! The result is the fundamental diagram of the system as presented in Wolfgang Knospe, Ludger Santen, Andreas Schadschneider, and Michael Schreckenberg. Towards a realistic microscopic description ofhighway traffic.J. Phys. A-Math Gen., 33(48):L477, 2000 figure 1.
  SUBROUTINE FUNDAMENTAL_DIAGRAM
    INTEGER :: i, j
    REAL(KIND=RDP), DIMENSION(0:9) :: velocity, localconcentration, flux 
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A23)")"FundamentalDiagramD.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    N = INT(0.10*(L/VEHICLE_LENGTH)) !Global density= 10%
    ALLOCATE(VEHICLES(-N:N+2))
    NUMBER_OF_DEFECTORS=N
    
    CALL InitialConditions           !initialization of the objects  
    DO i = 1, TRANSIENT_TIME         !Transient
      CALL AcelerationProcess
      CALL InteractionProcess
    END DO
    velocity = 0
    localconcentration = 0
    DO i = 1, MEASUREMENT_TIME              !measurements
      CALL AcelerationProcess
      CALL InteractionProcess
      DO j = 1, N
        !print*,i,j,MOD(VEHICLES(j)%POSITION,L),L  
        SELECT CASE(MOD(VEHICLES(j)%POSITION,L))
        CASE(0:INT(0.1*L))
          velocity(0) = velocity(0) + VEHICLES(j)%VELOCITY
          localconcentration(0) = localconcentration(0) + 1
        CASE(INT(0.1*L) + 1:INT(0.2*L))
          velocity(1) = velocity(1) + VEHICLES(j)%VELOCITY
          localconcentration(1) = localconcentration(1) + 1
        CASE(INT(0.2*L) + 1: INT(0.3*L))
          velocity(2) = velocity(2) + VEHICLES(j)%VELOCITY
          localconcentration(2) = localconcentration(2) + 1
        CASE(INT(0.3*L) + 1: INT(0.4*L))
          velocity(3) = velocity(3) + VEHICLES(j)%VELOCITY
          localconcentration(3) = localconcentration(3) + 1
        CASE(INT(0.4*L) + 1: INT(0.5*L))
          velocity(4) = velocity(4) + VEHICLES(j)%VELOCITY
          localconcentration(4) = localconcentration(4) + 1
        CASE(INT(0.5*L) + 1: INT(0.6*L))
          velocity(5) = velocity(5) + VEHICLES(j)%VELOCITY
          localconcentration(5) = localconcentration(5) + 1
        CASE(INT(0.6*L) + 1: INT(0.7*L))
          velocity(6) = velocity(6) + VEHICLES(j)%VELOCITY
          localconcentration(6) = localconcentration(6) + 1
        CASE(INT(0.7*L) + 1: INT(0.8*L))
        velocity(7) = velocity(7) + VEHICLES(j)%VELOCITY
          localconcentration(7) = localconcentration(7) + 1
        CASE(INT(0.8*L) + 1: INT(0.9*L))
          velocity(8) = velocity(8) + VEHICLES(j)%VELOCITY
          localconcentration(8) = localconcentration(8) + 1
        CASE(INT(0.9*L) + 1: L-1)
          velocity(9) = velocity(9) + VEHICLES(j)%VELOCITY
          localconcentration(9) = localconcentration(9) + 1
        END SELECT
      END DO
      DO j = 0, 9
        IF(localconcentration(j) == 0)THEN
          flux(j) = 0
        ELSE
          velocity(j) = velocity(j)/localconcentration(j)
          localconcentration(j) = (localconcentration(j)*VEHICLE_LENGTH)/(0.1*L)
          flux(j) = localconcentration(j)*velocity(j)
        END IF 
        WRITE(1,"(f5.3,A1,f5.3)") localconcentration(j), " ", flux(j)
      END DO
    END DO
    DEALLOCATE(VEHICLES)
  END SUBROUTINE FUNDAMENTAL_DIAGRAM
END MODULE PRIMARY_ROUTINES 

PROGRAM MAIN
  USE EXTERNAL_PARAMETERS 
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  USE SECONDARY_ROUTINES
  USE PRIMARY_ROUTINES
  IMPLICIT NONE
  ISEED = 3577
  CALL FUNDAMENTAL_DIAGRAM
END PROGRAM MAIN

