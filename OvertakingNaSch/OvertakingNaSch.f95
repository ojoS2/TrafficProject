! In this module we declare all constants used throughout the program as well as where we define abstract objects
! If one wish to vary some of these constants (the deceleration probability, for example) you may comment the respective line in this module and declare it in the GLOBAL_VARIABLES module
MODULE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  ! This is the object we work with. It has tree basic properties: its position; velocity and species 
  TYPE :: PARTICLES
    INTEGER :: POSITION
    INTEGER :: VELOCITY
    CHARACTER :: SPECIES 
  END TYPE PARTICLES
  
  INTEGER, PARAMETER :: L = 10000                      ! system size --> default is L=10^4
  INTEGER, PARAMETER :: NUMBER_OF_CONFIGURATIONS = 100 ! number of different configurations --> default N_c=100
  INTEGER, PARAMETER :: TERMALIZATION_TIME= 100000     ! time waited to start measurements --> default T_t=10^5
  ! We have divided the evolution of the system in timescales (after the transient) for especific investigations. The total run time in the esperiment is given by the multiplication of both scales, but set one of those as 1 and increase the other if you preffer (or alter the respective DO loops).  
  INTEGER, PARAMETER :: FIRST_TIME_SCALE = 1000        ! Time-step              
  INTEGER, PARAMETER :: SECOND_TIME_SCALE = 1000       ! 10^3 time-steps scale
  ! Experiment time scale is FIRST_TIME_SCALE x SECOND_TIME_SCALE --> default=10^6
  INTEGER, PARAMETER :: MAXIMUM_VELOCITY = 5           ! Vehicle maximum velocity
  INTEGER, PARAMETER :: RDP = selected_real_kind(10,1000)! Parameter necessary to big float opertions 
  REAL, PARAMETER :: P = 0.0                             ! Deceleration probability -->default=0.0
END MODULE EXTERNAL_PARAMETERS
! Declaration of the global variables 
MODULE GLOBAL_VARIABLES
  USE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  TYPE(PARTICLES), DIMENSION(:), ALLOCATABLE :: VEHICLES         
  LOGICAL, DIMENSION(0:L-1) :: STREET  
  INTEGER :: ISEED, N, DEFECTOS_NUMBER
  REAL(KIND=RDP) :: OVERTAKING_NUMBER
!  REAL(KIND=RDP) :: P
END MODULE GLOBAL_VARIABLES
! Functions used throughout the program 
MODULE FUNCTIONS  
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  IMPLICIT NONE
  CONTAINS
  !This function generates random numbers in the interval [0:1) and random integer seeds in the interval [1:134456)
  FUNCTION RandomGenerator(Seed)       
    INTEGER, INTENT(INOUT) :: Seed
    REAL :: RandomGenerator
    ISEED = mod(8121*seed+28411, 134456) 
    RandomGenerator = real(ISEED)/134456.        
    RETURN
  END FUNCTION RandomGenerator

END MODULE FUNCTIONS
! Secondary importance routine 
MODULE SECUNDARY_ROUTINES
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  IMPLICIT NONE
  CONTAINS
  ! This routine selects the position of a set (of size SizeOfSystem) of vehicles randomly without superposition 
  SUBROUTINE RandomPositionInitialization(SizeOfSystem)
    INTEGER, INTENT(IN) :: SizeOfSystem
    INTEGER :: aux, acum
    acum = 1                            
    DO WHILE(acum <= SizeOfSystem)
      aux = INT(L*RandomGenerator(ISEED))      
      OCUPADO: IF(.NOT.STREET(aux))THEN 
        STREET(aux) = .TRUE.   
        VEHICLES(acum)%POSITION = aux
        VEHICLES(acum)%VELOCITY = MAXIMUM_VELOCITY     
        acum = acum + 1                       
      END IF OCUPADO
    END DO 
  END SUBROUTINE RandomPositionInitialization
  !This routine selects vehicles randomly and initialize then as belonging to a given species
  !It requires the number of particles and the the number of defectors
  SUBROUTINE RandomSpeciesInitialization(NumberOfParticles, NumberOfdefectors)
    INTEGER, INTENT(IN) :: NumberOfParticles, NumberOfdefectors
    INTEGER :: aux, acum, cooperators, defectors
    defectors = NumberOfdefectors
    cooperators = NumberOfParticles - NumberOfdefectors
    IF (cooperators<=defectors)THEN  
      VEHICLES%SPECIES = "D"   
      acum = 1                            
      DO WHILE(acum<=cooperators)
        aux = INT(NumberOfParticles*RandomGenerator(ISEED))
        IF (aux > 0)THEN
          IF(VEHICLES(aux)%SPECIES == "D")THEN 
            VEHICLES(aux)%SPECIES = "C" 
            acum = acum + 1                       
          END IF 
        END IF 
      END DO 
    ELSE  IF(cooperators>defectors)THEN  
      VEHICLES%SPECIES= "C"
      acum = 1                            
      DO WHILE(acum<=defectors)
        aux = INT(NumberOfParticles*RandomGenerator(ISEED))
        IF (aux>0)THEN
          IF (VEHICLES(aux)%SPECIES == "C") THEN 
            VEHICLES(aux)%SPECIES = "D" 
            acum = acum + 1                      
          END IF 
        END IF 
      END DO 
    END IF 
  END SUBROUTINE RandomSpeciesInitialization 
  !This subroutine organize the particles in positional order using the bubble algorithm
  !It requares the initial position and the final position of the piece of the VEHICLES object which one wish to order. 
  SUBROUTINE ParticlesReordering(Beginnig, Ending)
    INTEGER, INTENT(IN) :: Beginnig, Ending
    LOGICAL :: Flag
    INTEGER :: i    
    TYPE(PARTICLES) :: auxiliaryVector       
    DO               
      Flag = .TRUE.               
      DO i = Beginnig, Ending - 1
        IF( VEHICLES(i)%POSITION >  VEHICLES(i + 1)%POSITION)THEN
          auxiliaryVector             = VEHICLES(i)      
          VEHICLES(i)       = VEHICLES(i + 1)
          VEHICLES(i + 1)   = auxiliaryVector
          Flag = .FALSE.                          
        END IF    
      END DO 
      ! All the particles' position are tested and only when all are ordered the main DO loop ends 
      IF(Flag)THEN
        EXIT
      END IF  
    END DO 
  END SUBROUTINE ParticlesReordering
  
END MODULE SECUNDARY_ROUTINES
! Primary importance routines
MODULE PRIMARY_ROUTINES 
 USE EXTERNAL_PARAMETERS
 USE GLOBAL_VARIABLES
 USE FUNCTIONS
 USE SECUNDARY_ROUTINES
 IMPLICIT NONE
  CONTAINS 
  !This routine initializes the system using the routineS of THE SECUNDARY_ROUTINES module
  SUBROUTINE RandomInitialConditions
    INTEGER :: i, j, gap
    STREET = .FALSE.                                    ! Initialization of the street as empty
    CALL RandomPositionInitialization(N)                ! Initialization of positions and ocupation of the street
    CALL RandomSpeciesInitialization(N, DEFECTOS_NUMBER)          ! Initialization of the species       
    CALL ParticlesReordering( 1, N )                    ! Ordering of the vehicles
    
    ! This piece of code make the initial disposition of space and velocities of the vehicles compactible with the NaSch model implemented  with periodic boundary conditions
    
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L   
    DO i = N, 1, -1
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - 1
      IF(gap < MAXIMUM_VELOCITY)THEN
        VEHICLES(i)%VELOCITY = gap
      ELSE
        VEHICLES(i)%VELOCITY = MAXIMUM_VELOCITY
      END IF   
    END DO    
    
    ! Now comes the overall boundary conditions (Notice that the vehicles are already ordered)
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L 
    DO i = 0, -N, -1
      VEHICLES(i) = VEHICLES(N + i)
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION - L 
    END DO 
    DO i = -N, N + 1
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION + 2*L 
    END DO 
  END SUBROUTINE RandomInitialConditions
  ! This is the main routine in this code as it display a possible implementation of our algorithm. The algorithm is as following
  
  ! 0 >> The initial vehicle is chosen 
  
  ! 1 >> step all vehicles have its velocity increased by one if their velocity is not at the the maximum value :: v_n(t+1)=MIN(v_n(t)+1,V_MAX)
  
  ! 2 >> let the space available to the vehicle be  gap_1 = x_(n+1)-x_n-1  then one in two paths are chosen --> if the vehicle is a cooperator, then v_n(t+1)=MIN(v_n(t+1),gap);  if the vehicle is a defector, then the following inequality is tested :  x_(n+j) + v_(n+j)(t+1) < x_n + v_n < x_(n+j+1). If the first inequality is not obeyed, then we say that the overtaking is impossible, v_n(t+1)=MIN(v_n(t+1),gap) , and go to the step 3. Else we do the desceleration step here  v_n(t+1)=v_n(t+1)-1 with probability p  and test the both inequalities again. If the inequalities are both obeyed for some j>0 then the overtaking is succesfull. Else we make j=j-1 , update the velocity as v_n(t+1)=MIN(v_n(t+1), x_(n+j)-x_(n+j-1)-1), run the randomization step:v_n(t+1)=v_n(t+1)-1 with probability p  and try again. If j becomes j=0 then v_n(t+1)=MIN(v_n(t+1),gap), we say the overtaking is impossible and we go to the tird step
  
  ! 3 >> This step aplies only to cooperators and to defectors imcapable of overtaking in this round: v_n(t+1)=v_n(t+1)-1 with probability p
  
  ! 4 >> Update the position of the vehicles : X(t+1) = X(t) + V(t+1)
  SUBROUTINE OvertakingNaSch
    INTEGER :: i, gap, aux, beginning, ending
    ! The road is initialized as empty, but following we mark the position of all vehicles  (remembering that (j-N)-th vehicle is identical to the j-th vehicle, and  the (N+1)-th vehicle is identical to the 1). This way deos not matter wheter we update with the (j-N)-th or the j-th
    STREET = .FALSE.
    DO i = 1, N
      STREET(MOD((VEHICLES(i)%POSITION), L)) = .TRUE.
    END DO 
    ! Now we choose the first vehicle to be updated as the vehicle behind the fastest on the street, because the the fastest vehicle can not be overtaken 
    aux = VEHICLES(N)%VELOCITY
    ending = N
    DO i = N - 1, 1, -1
      IF(VEHICLES(i)%VELOCITY > aux)THEN
        aux = VEHICLES(i)%VELOCITY
        ending = i
      END IF 
      IF(aux == MAXIMUM_VELOCITY)EXIT 
    END DO
    ending = ending - 1
    beginning = ending - N + 1
    ! We have our segment, now we perform the first step
    DO i = beginning, ending
      IF(VEHICLES(i)%VELOCITY < MAXIMUM_VELOCITY)VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY + 1
    END DO
    !The second step is interwined with the tird step. Notice that we have named the first and the last vehicle in positional order and not in the flux order. With this we enphatsize the directional character of the algorithm 
    MainLoop: DO i = ending, beginning, -1
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - 1
      !This is useful for the evaluation of the first inequality x_(n+1) + v_(n+1)(t+1) < x_n + v_n 
      aux = gap + VEHICLES(i + 1)%VELOCITY + 1 
      ! The overtaking does not happen
      IF(VEHICLES(i)%VELOCITY <= aux)THEN  
        !Second NaSch step
        IF(VEHICLES(i)%VELOCITY>gap)THEN
          VEHICLES(i)%VELOCITY=gap
        END IF 
        !Tird NaSch step 
        IF(VEHICLES(i)%VELOCITY>0)THEN 
          IF(RandomGenerator(ISEED)<P)THEN   
            VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1 
          END IF   
        END IF 
        ! Now we mark the sites ocupied by the by the vehicle during its dynamics
        ! This sites cannot be crossed by other vehicles 
        DO aux=1,VEHICLES(i)%VELOCITY
          STREET(MOD((VEHICLES(i)%POSITION + aux), L)) = .TRUE.
        END DO  
        CYCLE MainLoop
      END IF           
      ! The first inequality is satisfied, now we verify the vehicle species
      IF(VEHICLES(i)%SPECIES == "C")THEN
        VEHICLES(i)%VELOCITY = gap
        IF(VEHICLES(i)%VELOCITY > 0)THEN 
          IF(RandomGenerator(ISEED)<P)THEN
            VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1 
          END IF   
        END IF 
        DO aux = 1,VEHICLES(i)%VELOCITY
          STREET(MOD((VEHICLES(i)%POSITION + aux), L)) = .TRUE.
        END DO  
        CYCLE MainLoop
      END IF 
      ! The first inequality is satisfied and it is a defector. We perform the tird NaSch step
      IF(VEHICLES(i)%VELOCITY > 0)THEN 
        IF(RandomGenerator(ISEED)<P)THEN
          VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1 
        END IF   
      END IF 
      aux = VEHICLES(i)%POSITION + VEHICLES(i)%VELOCITY
      
      ! We are now in position to verify the second inequality
      DO
        !The inequality is verified and street space is marked 
        IF(.NOT.STREET(MOD(aux, L)))THEN
          DO aux=1,VEHICLES(i)%VELOCITY
            STREET(MOD((VEHICLES(i)%POSITION + aux), L)) = .TRUE.
          END DO  
          OVERTAKING_NUMBER = OVERTAKING_NUMBER + 1
          CYCLE MainLoop
        END IF   
        !The overtaking fails completly and the vehicle must behave as a cooperator
        IF(VEHICLES(i)%VELOCITY <= gap)THEN
          VEHICLES(i)%VELOCITY = gap
          IF(VEHICLES(i)%VELOCITY > 0)THEN 
            IF(RandomGenerator(ISEED)<p)THEN
              VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1 
            END IF   
          END IF 
          DO aux = 1,VEHICLES(i)%VELOCITY
            STREET(MOD((VEHICLES(i)%POSITION + aux), L)) = .TRUE.
          END DO  
          CYCLE MainLoop
        END IF   
        ! The inequality fails, but we have enough velocity to try again. Notice the re-run of the tird NaSch step
        VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1
        IF(VEHICLES(i)%VELOCITY > 0)THEN 
          IF(RandomGenerator(ISEED)<P)THEN
            VEHICLES(i)%VELOCITY = VEHICLES(i)%VELOCITY - 1 
          END IF   
        END IF 
        aux = VEHICLES(i)%POSITION + VEHICLES(i)%VELOCITY
      END DO
    END DO MainLoop
    ! Once updated the velocities we perform the updated of the velocities only on the selected section of the VHICLES object
    DO i = beginning, ending
      VEHICLES(i)%POSITION = VEHICLES(i)%VELOCITY + VEHICLES(i)%POSITION
    END DO
    ! We reorder 
    CALL ParticlesReordering(beginning, ending) 
    ! Finaly we update the other vehicles 
    DO i = ending + 1, N+1
      VEHICLES(i) = VEHICLES(i - N)
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION + L
    END DO 
    DO i = beginning -1, - N, -1
      VEHICLES(i) = VEHICLES(i + N )
      VEHICLES(i)%POSITION = VEHICLES(i)%POSITION - L
    END DO 
  END SUBROUTINE OvertakingNaSch

  ! This routine prints the fundamental diagram    
  SUBROUTINE FundamentalDiagram
    INTEGER :: i, j, k, configurations
    REAL(KIND = RDP) :: velocity, aux
    CHARACTER(LEN=100) :: Arq1, Arq2
    WRITE(Arq2,"(A20)")"OvertakingN0.000.dat"
    OPEN(2,FILE=TRIM(Arq2))   
    WRITE(Arq1,"(A27)")"FundamentalDiagram0.000.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    DO N=INT(0.01*L),INT(0.99*L),INT(0.01*L) ! Fundamental Diagram with 100 points 
      ALLOCATE(VEHICLES(-N : N + 1)) 
      DEFECTOS_NUMBER = N                     ! ALLD Population, It may be changed
      velocity = 0.
      OVERTAKING_NUMBER = 0
      DO configurations = 1, NUMBER_OF_CONFIGURATIONS
        CALL RandomInitialConditions
        ! Termalization of the system
        aux = OVERTAKING_NUMBER
        DO i = 1,TERMALIZATION_TIME
          CALL OvertakingNaSch
        END DO
        OVERTAKING_NUMBER = aux
        DO k = 1, SECOND_TIME_SCALE
          DO j = 1,FIRST_TIME_SCALE  
            CALL OvertakingNaSch
            DO i = 1, N
              velocity = velocity + VEHICLES(i)%VELOCITY
            END DO 
          END DO  
        END DO   
      END DO 
      ! the averaging over the velocity occurs in two steps because of float point errors 
      velocity = velocity/(N*SECOND_TIME_SCALE)
      velocity = velocity/(FIRST_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
      OVERTAKING_NUMBER = OVERTAKING_NUMBER/(SECOND_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
      OVERTAKING_NUMBER = OVERTAKING_NUMBER/(FIRST_TIME_SCALE*N)
      DEALLOCATE(VEHICLES)
      WRITE(1,*)REAL(N)/L,(REAL(N)/L)*velocity,velocity
      WRITE(2,*)REAL(N)/L,OVERTAKING_NUMBER
      print*, velocity,real(N)/L, OVERTAKING_NUMBER
    END DO  
    CLOSE(1)
    CLOSE(2)
  END SUBROUTINE FundamentalDiagram
  
  ! This routine is a visual aid. It prints a spatio-temporal pattern of the system in a matrix form that can be visualized in gnuplot with the following code, for example 
  ! reset
  ! set terminal qt 
  ! set pm3d map
  ! splot "matrix_file.dat" matrix
  
  ! This code may present problems when the matrix is very big, we recomend a size of L=1.000.  
  SUBROUTINE SpatioTemporalPattern
    INTEGER :: i, j, k
    INTEGER, DIMENSION(0:L-1) :: streetocupation
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A24)")"SpatioTemporalMatrix.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    ! specify the density
    N = 0.2*L
    ! specify the number of defectors
    DEFECTOS_NUMBER=N
    ALLOCATE(VEHICLES(-N : N + 1)) 
    CALL RandomInitialConditions
    ! We work with a square matrix. If you preffer other formats you may change the range of the loop below
    DO j=1,L
      ! this will make the empty space yellow and the occupied space purple in the gnuplot code presented above, you may invert the output colors inverting streetocupation=0 below or even grad the pattern according to the velocity as streetocupation(MOD(VEHICLES(i)%POSITION,L))=VEHICLES(i)%VELOCITY
      streetocupation=9
      DO i=1,N
        ! Here we mark the position of the vehicles
        streetocupation(MOD(VEHICLES(i)%POSITION,L))=0
      END DO  
      ! impression of a line in the matrix
      DO i = 0,L-1
        WRITE(1,"(I1,A1)",ADVANCE='NO')streetocupation(i)," "
      END DO
      WRITE(1,"(A1)")" "
      !update of the system
      CALL OvertakingNaSch
    END DO 
    CLOSE(1)
    DEALLOCATE(VEHICLES)
  END SUBROUTINE SpatioTemporalPattern

 
END MODULE PRIMARY_ROUTINES 

PROGRAM MAIN
  USE EXTERNAL_PARAMETERS 
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  USE SECUNDARY_ROUTINES
  USE PRIMARY_ROUTINES
  IMPLICIT NONE
  ! Initiate the random seed 
  ISEED = 1313
  
  ! Uncoment one of this lines to call the main functions of this code
  ! but beware. The first SpatioTemporalPattern routine is fast. Nevertheless, for the very modest parameters L=1000 and NUMBER_OF_CONFIGURATIONS=10 the FundamentalDiagram routine cost 3.5 hours in my personal computer which attributes are (output from $ inxi -Fxz )
  ! Processes: 253 Uptime: 5d 21h 12m Memory: 5.72 GiB used: 4.60 GiB (80.5%) Init: systemd runlevel: 5 Compilers: 
               !gcc: 9.2.1 Shell: bash v: 5.0.3 inxi: 3.0.36
  !CPU:       Topology: Dual Core model: Intel Core i5-2410M bits: 64 type: MT MCP arch: Sandy Bridge rev: 7 L2                    cache: 3072 KiB 
               !flags: avx lm nx pae sse sse2 sse3 sse4_1 sse4_2 ssse3 vmx bogomips: 18359 
               ! Speed: 1796 MHz min/max: 800/2300 MHz Core speeds (MHz): 1: 1796 2: 1796 3: 1796 4: 1796              
  !CALL SpatioTemporalPattern
  CALL FundamentalDiagram
END PROGRAM MAIN



