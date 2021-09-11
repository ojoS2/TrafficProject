!********************************************************************************************************************************************
!********************************************************************************************************************************************
! This is a FORTRAN program used to obtain the results in the work entitled "Evolution of behaviors in heterogeneous traffic models as driven annealed disorders and its relation to the n-vector model". Its language was made to simple to conceptualize and not to minimize computational costs. Feel free to modify it to your heart's content, but keep it open so everybody may contribute to a fast scientific advancement 
!********************************************************************************************************************************************
!********************************************************************************************************************************************!********************************************************************************************************************************************



! In this module we declare all constants used throughout the program as well and define complex objects
! If one wish to vary some of these constants in the routines to make some graphs, you may comment the respective line in this module and declare it in the GLOBAL_VARIABLES module

MODULE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  ! This is the object we work with in multibehavioral population analysis. It has two basic properties regarding variables: its position and, velocity. Also  10 properties dealing with internal parameters of the algorithms  that are changed by imitation and mutation.
TYPE :: PARTICLES
    ! Dynamical quantities
    INTEGER :: POSITION
    INTEGER :: VELOCITY
    ! Behavioral parameters
    INTEGER :: MAXIMUM_VELOCITY 
    INTEGER :: TEMPORAL_STATUS
    REAL :: P
    REAL :: P_F
    REAL :: P_S
    REAL :: P_0
    REAL :: P_ACC
    REAL :: P_DCC
    REAL :: X_T
    REAL :: X_s
    INTEGER :: SP
  END TYPE PARTICLES
    
  INTEGER, PARAMETER :: L = 10000                                          ! system size --> default is L=10^4
  INTEGER, PARAMETER :: NUMBER_OF_CONFIGURATIONS = 1  ! number of different configurations --> default N_c=100
  INTEGER, PARAMETER :: TERMALIZATION_TIME= 100000           ! time waited to start measurements --> default T_t=10^5
  
  ! I have divided the evolution of the system in timescales (after the transient) for especific investigations. The total run time in the esperiment is given by the multiplication of both scales, but set one of those as 1 and increase the other if you preffer (or alter the respective DO loops).  
  
  INTEGER, PARAMETER :: FIRST_TIME_SCALE = 1000                  !  Time-step  scales            
  INTEGER, PARAMETER :: SECOND_TIME_SCALE = 1000             ! next time-scale
  INTEGER, PARAMETER :: IMITATION_TIME_SCALE = 1000           ! Time-scale related to the evaluation of strategies in the imitation process 
  INTEGER, PARAMETER :: TIME_SERIES_SCALE = 10000              ! Number of points taken in a time series analysis
  INTEGER, PARAMETER :: TRAJECTORIES_SCALE = 10000           ! Number of points on trajectories in the the space of parameters 
  INTEGER, PARAMETER :: EVOLUTION_SCALE = 1000000            ! Duration of evolution considered 
  INTEGER, PARAMETER :: RDP = selected_real_kind(10,1000) ! Parameter necessary to big float opertions 
  REAL, PARAMETER :: CONTACT_FRQUENCY_MINIMUM = 0.8      ! Related to the imitation process
  
  ! Now I will initialize the limits of the parameters
  INTEGER, PARAMETER :: MAX_VEL_MAX_LIMIT = 9
  INTEGER, PARAMETER :: MAX_VEL_MIN_LIMIT = 4
  REAL, PARAMETER :: P_MAX_LIMIT = .5
  REAL, PARAMETER :: P_MIN_LIMIT = .1
  REAL, PARAMETER :: P_F_MAX_LIMIT = .5
  REAL, PARAMETER :: P_F_MIN_LIMIT = .0
  REAL, PARAMETER :: P_S_MAX_LIMIT = .95
  REAL, PARAMETER :: P_S_MIN_LIMIT = .9
  REAL, PARAMETER :: P_0_MAX_LIMIT = .7
  REAL, PARAMETER :: P_0_MIN_LIMIT = .4
  REAL, PARAMETER :: P_ACC_MAX_LIMIT = .95
  REAL, PARAMETER :: P_ACC_MIN_LIMIT = .9
  REAL, PARAMETER :: P_DCC_MAX_LIMIT = .95
  REAL, PARAMETER :: P_DCC_MIN_LIMIT = .8
  REAL, PARAMETER :: X_T_MAX_LIMIT = 1.3
  REAL, PARAMETER :: X_T_MIN_LIMIT = 1.2
  REAL, PARAMETER :: X_s_MAX_LIMIT = 2.
  REAL, PARAMETER :: X_s_MIN_LIMIT = 0.
  
  
END MODULE EXTERNAL_PARAMETERS

!Declaration of the global variables. These variables are used and modified by routines through the program so they should be used with caution 

MODULE GLOBAL_VARIABLES
  USE EXTERNAL_PARAMETERS
  IMPLICIT NONE
  TYPE(PARTICLES), DIMENSION(:), ALLOCATABLE :: VEHICLES    ! Vector of particles
  TYPE(PARTICLES) :: TEMPVECTOR, AVTEMPVECTOR                   ! Used in [VectorBuilder] and [OrderParameter] routines
  LOGICAL, DIMENSION(0:L-1) :: STREET                                       ! Environment 
  INTEGER :: ISEED                                                                        ! Randon generator seed
  INTEGER :: N                                                                               ! Number of particles
  
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

! Secondary importance routines, generally related to the inicial conditions 

MODULE SECUNDARY_ROUTINES
  USE EXTERNAL_PARAMETERS
  USE GLOBAL_VARIABLES
  USE FUNCTIONS
  IMPLICIT NONE
  CONTAINS
  
  ! This routine selects the position of [SizeOfSystem] vehicles randomly without superposition.  
  
  SUBROUTINE RandomPositionInitialization(SizeOfSystem)
    INTEGER, INTENT(IN) :: SizeOfSystem
    INTEGER :: aux, acum
    acum = 1                            
    DO WHILE(acum <= SizeOfSystem)
      aux = INT(L*RandomGenerator(ISEED))      
      BUSY: IF(.NOT.STREET(aux))THEN 
        STREET(aux) = .TRUE.   
        VEHICLES(acum)%POSITION = aux
        VEHICLES(acum)%VELOCITY = VEHICLES(acum)%VELOCITY
        acum = acum + 1                       
      END IF BUSY
    END DO 
  END SUBROUTINE RandomPositionInitialization
  
  !This routine selects vehicles randomly and initialize then as belonging to a given species with their respective parameters
  !It requires the number of particles of each specie
  
  SUBROUTINE RandomSpeciesInitialization(N_NaSch, N_NaSchCC, N_SFI, N_TT,N_BJH,N_VDR,N_TOCA)
    INTEGER, INTENT(IN) :: N_NaSch, N_NaSchCC, N_SFI, N_TT, N_BJH, N_VDR, N_TOCA
    INTEGER :: Indice, aux
    VEHICLES%SP=0                                   ! All species initialized as a flag
    aux = 0                                                                    ! I add the species all at a time, notice that the sum of all entry parameters must be N  
    DO
       IF(aux==N_NaSch)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%SP==0)THEN
          VEHICLES(Indice)%SP = 1
          aux = aux+1
       END IF   
    END DO
    aux = 0
    DO
       IF(aux==N_NaSchCC)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%SP==0)THEN
          VEHICLES(Indice)%SP = 2
          aux = aux+1
       END IF   
    END DO
    aux = 0
    DO
       IF(aux==N_SFI)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%SP==0)THEN
          VEHICLES(Indice)%SP = 3
          aux = aux+1
       END IF   
    END DO
    aux = 0
    DO
       IF(aux==N_TT)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%SP==0)THEN
          VEHICLES(Indice)%SP = 4
          aux = aux+1
       END IF   
    END DO
    aux = 0
    DO
       IF(aux==N_BJH)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%SP==0)THEN
          VEHICLES(Indice)%SP = 5
          aux = aux+1
       END IF   
    END DO
    aux = 0
    DO
       IF(aux==N_VDR)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%SP==0)THEN
          VEHICLES(Indice)%SP = 6
          aux = aux+1
       END IF   
    END DO
    aux = 0
    DO
       IF(aux==N_TOCA)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%SP==0)THEN
          VEHICLES(Indice)%SP = 7
          aux = aux+1
       END IF   
    END DO

  END SUBROUTINE RandomSpeciesInitialization 
  
  
  !This subroutine organize the particles in positional order using the bubble algorithm. I made it keeping overtaking in mind, but one may change it for a les expensive code such as numbering them according to their positions.
  !It requares the initial position and the final position of the piece of the VEHICLES object which one wish to order. 
  
  
  SUBROUTINE ParticlesReordering(Beginnig, Ending)
    INTEGER, INTENT(IN) :: Beginnig, Ending
    LOGICAL :: Flag
    INTEGER :: i    
    TYPE(PARTICLES) :: auxiliaryVector       
    DO               
      Flag = .TRUE.               
      DO i = Beginnig, Ending - 1
        IF(VEHICLES(i)%POSITION >  VEHICLES(i + 1)%POSITION)THEN
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
  
  !  Inicialization of the raw parameters. Homogeneous populations
  
  SUBROUTINE HomogeneousParametersInicialization
    VEHICLES%MAXIMUM_VELOCITY=5
    VEHICLES%P=0.1
    VEHICLES%P_F=0.0
    VEHICLES%P_S=0.9
    VEHICLES%P_0=0.5
    VEHICLES%P_ACC=0.9
    VEHICLES%P_DCC=0.9
    VEHICLES%X_T=1.2
    VEHICLES%X_s=1
    VEHICLES%TEMPORAL_STATUS=1
  END SUBROUTINE HomogeneousParametersInicialization  
 
  !  Inicialization of the raw parameters. Heterogeneous populations
 
  SUBROUTINE HeterogeneousParametersInicialization
    INTEGER :: Indice
    DO Indice=1,N
      VEHICLES(Indice)%MAXIMUM_VELOCITY=INT(5*RandomGenerator(ISEED))+4
      VEHICLES(Indice)%P= 0.1+RandomGenerator(ISEED)/4
      VEHICLES(Indice)%P_F=RandomGenerator(ISEED)/2
      VEHICLES(Indice)%P_S=0.9+RandomGenerator(ISEED)/20
      VEHICLES(Indice)%P_0=0.4+RandomGenerator(ISEED)/3.5
      VEHICLES(Indice)%P_ACC=0.9+RandomGenerator(ISEED)/20
      VEHICLES(Indice)%P_DCC=0.8+RandomGenerator(ISEED)/8
      VEHICLES(Indice)%X_T= 1.2+RandomGenerator(ISEED)/10
      VEHICLES(Indice)%X_s=2*RandomGenerator(ISEED)
      VEHICLES(Indice)%TEMPORAL_STATUS=1
    END DO
  END SUBROUTINE HeterogeneousParametersInicialization  

  ! This routine is for analysis of the Heterogeneous NaSch model presnted in the section Ex 1. in the main text
  
  SUBROUTINE HeteroNaSchInicializationOfSpecies(c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8)
    REAL, INTENT(IN) :: c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8
    INTEGER :: Indice, aux
    VEHICLES%SP=1                                   ! All particles are NaSch particles
    VEHICLES%MAXIMUM_VELOCITY=0      ! The flag is in the maximum velocity
    aux = INT(c_1*N)                                                       
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 2
          VEHICLES(Indice)%P = 0.18
          aux = aux-1
       END IF   
    END DO
    
    aux = INT(c_2*N)                                                       
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 3
          VEHICLES(Indice)%P = 0.17
          aux = aux-1
       END IF   
    END DO
    
    aux = INT(c_3*N)                                                       
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 4
          VEHICLES(Indice)%P = 0.16
          aux = aux-1
       END IF   
    END DO
    
    aux = INT(c_4*N)                                                       
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 5
          VEHICLES(Indice)%P = 0.15
          aux = aux-1
       END IF   
    END DO
    
    aux = INT(c_5*N)                                                       
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 6
          VEHICLES(Indice)%P = 0.14
          aux = aux-1
       END IF   
    END DO
    
    aux = INT(c_6*N)                                                       
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 7
          VEHICLES(Indice)%P = 0.13
          aux = aux-1
       END IF   
    END DO
    
    aux = INT(c_7*N)                                                       
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 8
          VEHICLES(Indice)%P = 0.12
          aux = aux-1
       END IF   
    END DO
    
    ! This below is necessary for \sum (c_i) = 1  for small N
    aux = N - INT(c_1*N) - INT(c_2*N)- INT(c_3*N)- INT(c_4*N)- INT(c_5*N)- INT(c_6*N)- INT(c_7*N) 
    DO
       IF(aux==0)EXIT
       Indice=INT(N*RandomGenerator(ISEED))+1
       IF(VEHICLES(Indice)%MAXIMUM_VELOCITY==0)THEN
          VEHICLES(Indice)%MAXIMUM_VELOCITY = 9
          VEHICLES(Indice)%P = 0.11
          aux = aux-1
       END IF   
    END DO
  END SUBROUTINE HeteroNaSchInicializationOfSpecies

  !This routine builds auxiliary vectors used in the routine OrderParameter
  
  SUBROUTINE VectorBuilder(indice)
    INTEGER, INTENT(IN) :: Indice
    
    SELECT CASE (VEHICLES(indice)%SP)
    !(v_{max}, p, p_f, p_s, p_0, p_{acc}, p_{dec}, X_t, X_s, SP)
    CASE(1)  ! NaSch (v_max, p, 0, 0, 0, 0, 0, 0, 0, 1).  All unused parameters in this model get null values
    
      TEMPVECTOR%MAXIMUM_VELOCITY=VEHICLES(indice)%MAXIMUM_VELOCITY
      TEMPVECTOR%P=VEHICLES(indice)%P
      TEMPVECTOR%P_F=0
      TEMPVECTOR%P_S=0
      TEMPVECTOR%P_0=0
      TEMPVECTOR%P_ACC=0
      TEMPVECTOR%P_DCC=0
      TEMPVECTOR%X_T=0
      TEMPVECTOR%X_S=0
      TEMPVECTOR%SP=1
      
    CASE(2) ! NaSchCC (v_max, p, p_f, 0, 0, 0, 0, 0, 0, 2).  All unused parameters in this model get null values
    
      TEMPVECTOR%MAXIMUM_VELOCITY=VEHICLES(indice)%MAXIMUM_VELOCITY
      TEMPVECTOR%P=VEHICLES(indice)%P
      TEMPVECTOR%P_F=VEHICLES(indice)%P_F
      TEMPVECTOR%P_S=0
      TEMPVECTOR%P_0=0
      TEMPVECTOR%P_ACC=0
      TEMPVECTOR%P_DCC=0
      TEMPVECTOR%X_T=0
      TEMPVECTOR%X_S=0
      TEMPVECTOR%SP=2
    
    CASE(3)  !SFI (v_max, p, p_f, 0, 0, 0, 0, 0, 0, 3).  All unused parameters in this model get null values
      
      TEMPVECTOR%MAXIMUM_VELOCITY=VEHICLES(indice)%MAXIMUM_VELOCITY
      TEMPVECTOR%P=VEHICLES(indice)%P
      TEMPVECTOR%P_F=VEHICLES(indice)%P_F 
      TEMPVECTOR%P_S=0
      TEMPVECTOR%P_0=0
      TEMPVECTOR%P_ACC=0
      TEMPVECTOR%P_DCC=0
      TEMPVECTOR%X_T=0
      TEMPVECTOR%X_S=0
      TEMPVECTOR%SP=3
      
    CASE(4) !TT (v_max, 0, 0, 0, 0, 0, 0, 0, X_s, 4).  All unused parameters in this model get null values
    
      TEMPVECTOR%MAXIMUM_VELOCITY=VEHICLES(indice)%MAXIMUM_VELOCITY
      TEMPVECTOR%P=0
      TEMPVECTOR%P_F=0 
      TEMPVECTOR%P_S=0
      TEMPVECTOR%P_0=0
      TEMPVECTOR%P_ACC=0
      TEMPVECTOR%P_DCC=0
      TEMPVECTOR%X_T=0
      TEMPVECTOR%X_S=VEHICLES(indice)%X_S
      TEMPVECTOR%SP=4
      
    CASE(5) !BJH (v_max, p, 0, p_s, 0, 0, 0, 0, 0, 5).  All unused parameters in this model get null values
 
      TEMPVECTOR%MAXIMUM_VELOCITY=VEHICLES(indice)%MAXIMUM_VELOCITY
      TEMPVECTOR%P=VEHICLES(indice)%P
      TEMPVECTOR%P_F=0
      TEMPVECTOR%P_S=VEHICLES(indice)%P_S
      TEMPVECTOR%P_0=0
      TEMPVECTOR%P_ACC=0
      TEMPVECTOR%P_DCC=0
      TEMPVECTOR%X_T=0
      TEMPVECTOR%X_S=0
      TEMPVECTOR%SP=5
      
    CASE(6) !VDR (v_max, p, 0, 0, p_0, 0, 0, 0, 0, 6),  All unused parameters in this model get null values
      
      TEMPVECTOR%MAXIMUM_VELOCITY=VEHICLES(indice)%MAXIMUM_VELOCITY
      TEMPVECTOR%P=VEHICLES(indice)%P
      TEMPVECTOR%P_F=0 
      TEMPVECTOR%P_S=0
      TEMPVECTOR%P_0=VEHICLES(indice)%P_0
      TEMPVECTOR%P_ACC=0
      TEMPVECTOR%P_DCC=0
      TEMPVECTOR%X_T=0
      TEMPVECTOR%X_S=0
      TEMPVECTOR%SP=6
      
    CASE(7)  !TO (v_max, 0, 0, 0, 0, p_acc, p_dcc, X_t, 0, 7),  All unused parameters in this model get null values
      
      TEMPVECTOR%MAXIMUM_VELOCITY=VEHICLES(indice)%MAXIMUM_VELOCITY
      TEMPVECTOR%P=VEHICLES(indice)%P
      TEMPVECTOR%P_F=0 
      TEMPVECTOR%P_S=0
      TEMPVECTOR%P_0=0
      TEMPVECTOR%P_ACC=VEHICLES(indice)%P_ACC
      TEMPVECTOR%P_DCC=VEHICLES(indice)%P_DCC
      TEMPVECTOR%X_T=VEHICLES(indice)%X_T
      TEMPVECTOR%X_S=0
      TEMPVECTOR%SP=7
      
    CASE DEFAULT
      PRINT*,"Error, invalid selection[VectorBuilder]"
      STOP
    END SELECT
  END SUBROUTINE VectorBuilder

  
END MODULE SECUNDARY_ROUTINES

! Primary importance routines


MODULE PRIMARY_ROUTINES 
 USE EXTERNAL_PARAMETERS
 USE GLOBAL_VARIABLES
 USE FUNCTIONS
 USE SECUNDARY_ROUTINES
 IMPLICIT NONE
  CONTAINS 
  
  
  !Initial conditions of the Heterogeneou NaSch model
  
  SUBROUTINE HeteroNaSchInitialConditions
    INTEGER :: i, gap
    STREET = .FALSE.                                             ! Initialization of the environment as empty (see {RandomPositionInitialization} routine)
    CALL RandomPositionInitialization(N)                ! Initialization of positions and ocupation of the street
                                                     
    CALL HeteroNaSchInicializationOfSpecies(.01,0.0,0.0,0.0,0.0,0.0,0.0,0.99)!
                                                                               ! Inicialization of the NaSch species in the set (v_max=1+n, p = 0.19-0.01n), n in [1,8]
                                                                               ! each entry represents the n-th position in the set and all entries must sum to 1
                                                                               ! put 1 in one entry and 0 in the rest to have a homogeneous population
                                                                               
    CALL ParticlesReordering( 1, N )                        ! Ordering of the vehicles
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L   
    DO i = N, 1, -1
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - 1
      IF(gap < VEHICLES(i)%MAXIMUM_VELOCITY)THEN
        VEHICLES(i)%VELOCITY = gap
      ELSE
        VEHICLES(i)%VELOCITY = VEHICLES(i)%MAXIMUM_VELOCITY
      END IF   
    END DO    
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L 
    
  END SUBROUTINE HeteroNaSchInitialConditions

  
  ! Initial conditions of multibehavioral populations 
  
  SUBROUTINE MultiBehavioralRandomInitialConditions
    INTEGER :: i, gap
    STREET = .FALSE.                                             ! Initialization of the environment as empty (see {RandomPositionInitialization} routine)
    CALL HeterogeneousParametersInicialization   ! Inicialization of the intrinsic parameters of all models
    !CALL HomogeneousParametersInicialization
    CALL RandomPositionInitialization(N)                ! Initialization of positions and ocupation of the street
    CALL RandomSpeciesInitialization(INT(0.1429*N),INT(0.1429*N),INT(0.1429*N),INT(0.1429*N),INT(0.1429*N),&
                                                          &INT(0.1429*N),N-6*INT(0.1429*N))      ! Initialization of the species. The sum of the arguments in this routine must be N
                                                          !(0,0,0,0,0,0,N)!
                                                                               ! Put N in any entry and 0 in the rest to simulate a single model
                                                                               
    CALL ParticlesReordering( 1, N )                        ! Ordering of the vehicles
    
    ! This piece of code make the initial disposition of space and velocities of the vehicles compactible with the NaSch model implemented  with periodic boundary conditions
    
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L   
    DO i = N, 1, -1
      gap = VEHICLES(i + 1)%POSITION - VEHICLES(i)%POSITION - 1
      IF(gap < VEHICLES(i)%MAXIMUM_VELOCITY)THEN
        VEHICLES(i)%VELOCITY = gap
      ELSE
        VEHICLES(i)%VELOCITY = VEHICLES(i)%MAXIMUM_VELOCITY
      END IF   
    END DO    
    
    ! Now comes the overall boundary conditions (Notice that the vehicles are already ordered)
    
    VEHICLES(N + 1) = VEHICLES(1)
    VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L 
    
  END SUBROUTINE MultiBehavioralRandomInitialConditions
  
  ! Now the algorithm of the models. They process their velocities particle by particle 
  
  SUBROUTINE NaSch(Indice)
    INTEGER, INTENT(IN) :: Indice
    INTEGER :: Aux
    IF(VEHICLES(Indice)%VELOCITY+1<=VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
      VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY+1
    END IF
    Aux=VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1
    IF(VEHICLES(Indice)%VELOCITY>Aux)THEN
      VEHICLES(Indice)%VELOCITY=Aux
    END IF 
    IF(VEHICLES(Indice)%VELOCITY>0)THEN
      IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P)THEN
        VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY-1  
      END IF
    END IF   
  END SUBROUTINE NaSch
  
  SUBROUTINE NaSchCC(Indice)
    INTEGER, INTENT(IN) :: Indice
    INTEGER :: Aux
    IF(VEHICLES(Indice)%VELOCITY+1<=VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
      VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY+1
    END IF
    Aux=VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1
    IF(VEHICLES(Indice)%VELOCITY>Aux)THEN
      VEHICLES(Indice)%VELOCITY=Aux
    END IF 
    IF(VEHICLES(Indice)%VELOCITY==VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
       IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P_F)THEN
        VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY-1  
      END IF
    ELSE
      IF(VEHICLES(Indice)%VELOCITY>0)THEN
        IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P)THEN
          VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY-1  
        END IF
      END IF
    END IF  
  END SUBROUTINE NaSchCC

  SUBROUTINE SFI(Indice)
    INTEGER, INTENT(IN) :: Indice
    INTEGER :: Aux
    VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY
    Aux=VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1
    IF(VEHICLES(Indice)%VELOCITY>Aux)THEN
      VEHICLES(Indice)%VELOCITY=Aux
    END IF 
    IF(VEHICLES(Indice)%VELOCITY==VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
      IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P)THEN
        VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY-1  
      END IF
    END IF   
  END SUBROUTINE SFI
  
  SUBROUTINE TT(Indice)
    INTEGER, INTENT(IN) :: Indice
    INTEGER :: Aux
    Aux=VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1
    IF(Aux>INT(VEHICLES(Indice)%X_s*VEHICLES(Indice)%MAXIMUM_VELOCITY   ))THEN
      VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY
    END IF   
    IF(VEHICLES(Indice)%VELOCITY>Aux)THEN
       VEHICLES(Indice)%VELOCITY=Aux
    END IF 
  END SUBROUTINE TT
  
  SUBROUTINE BJH(Indice)
    INTEGER, INTENT(IN) :: Indice
    INTEGER :: Aux
    IF(VEHICLES(Indice)%VELOCITY+1<=VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
      IF(VEHICLES(Indice)%VELOCITY==0)THEN
        IF(VEHICLES(Indice)%TEMPORAL_STATUS==1)THEN
          IF(RandomGenerator(ISEED)<(1-VEHICLES(Indice)%P_S))THEN
            VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY+1  
            VEHICLES(Indice)%TEMPORAL_STATUS=1
          ELSE 
            VEHICLES(Indice)%TEMPORAL_STATUS=0
          END IF
        ELSE IF(VEHICLES(Indice)%TEMPORAL_STATUS==0)THEN 
          IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P_S)THEN
            VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY+1  
            VEHICLES(Indice)%TEMPORAL_STATUS=1
          ELSE 
            VEHICLES(Indice)%TEMPORAL_STATUS=0
          END IF
        END IF   
      ELSE 
        VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY+1
      END IF 
    END IF
    
    Aux=VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1
    IF(VEHICLES(Indice)%VELOCITY>Aux)THEN   
      VEHICLES(Indice)%VELOCITY=Aux
    END IF 
    IF(VEHICLES(Indice)%VELOCITY>0)THEN
      IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P)THEN
        VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY-1  
      END IF
    END IF   
  END SUBROUTINE BJH

  SUBROUTINE VDR(Indice)
    INTEGER, INTENT(IN) :: Indice
    INTEGER :: Aux
    REAL :: p_eff
    IF(VEHICLES(Indice)%VELOCITY==0)THEN
       p_eff=VEHICLES(Indice)%P_0
    ELSE 
       p_eff=VEHICLES(Indice)%P
    END IF 
    IF(VEHICLES(Indice)%VELOCITY+1<=VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
      VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY+1
    END IF
    Aux=VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1
    IF(VEHICLES(Indice)%VELOCITY>Aux)THEN
      VEHICLES(Indice)%VELOCITY=Aux
    END IF 
    IF(VEHICLES(Indice)%VELOCITY>0)THEN
      IF(RandomGenerator(ISEED)<p_eff)THEN
        VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY-1  
      END IF
    END IF   
  END SUBROUTINE VDR

  SUBROUTINE TOCA(Indice)
    INTEGER, INTENT(IN) :: Indice
    INTEGER :: Aux
    
    Aux=VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1
    IF(VEHICLES(Indice)%VELOCITY<VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
      IF(Aux>VEHICLES(Indice)%X_T*VEHICLES(Indice)%VELOCITY)THEN
        IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P_ACC)THEN
          VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY+1  
        END IF 
      END IF 
    END IF 
    IF(VEHICLES(Indice)%VELOCITY>Aux)THEN
      VEHICLES(Indice)%VELOCITY=Aux
    END IF 
    IF(VEHICLES(Indice)%VELOCITY>0)THEN
      IF(Aux<VEHICLES(Indice)%X_T*VEHICLES(Indice)%VELOCITY)THEN
        IF(RandomGenerator(ISEED)<VEHICLES(Indice)%P_DCC)THEN
          VEHICLES(Indice)%VELOCITY=VEHICLES(Indice)%VELOCITY-1  
        END IF 
      END IF 
    END IF   
  END SUBROUTINE TOCA

 ! Here we perform the steps to any model and update their positions and velocities. In the next version of this algorithm, I present a variation where %SP is a real number related to probabilities. 
  
  SUBROUTINE PositionUpdate
     INTEGER :: Indice, Aux
     DO Indice = 1,N
       CALL NaSch(Indice) ! If it is only one strategy, uncoment this line and (change the model) coment the rest to save time 
       Aux=VEHICLES(Indice)%SP 
       SELECT CASE(INT(VEHICLES(Indice)%SP))
       CASE(1)
            CALL NaSch(Indice)
       CASE(2)
            CALL NaSchCC(Indice)
       CASE(3)
            CALL SFI(Indice)
       CASE(4)
            CALL TT(Indice)
       CASE(5)
            CALL BJH(Indice)
       CASE(6)
            CALL VDR(Indice)
       CASE(7)
            CALL TOCA(Indice)
       CASE DEFAULT
         PRINT*,"Selection out of Scope",  Aux,VEHICLES(Indice)%SP 
       END SELECT
     END DO 
     
     DO Indice = 1,N
       VEHICLES(Indice)%POSITION=VEHICLES(Indice)%POSITION+VEHICLES(Indice)%VELOCITY
     END DO
     
     VEHICLES(N + 1) = VEHICLES(1)
     VEHICLES(N + 1)%POSITION = VEHICLES(N + 1)%POSITION + L 
     
  END SUBROUTINE PositionUpdate

 ! Imitation dynamics. This is the cern of the algorithm.
 
  SUBROUTINE ImitationDynamics(Indice, VelFocal, VelTarget, ContactFrequency)
    INTEGER, INTENT(IN) :: Indice
    REAL, INTENT(OUT) :: VelFocal, VelTarget, ContactFrequency
    INTEGER  :: I
    VelFocal = 0
    VelTarget = 0
    ContactFrequency=0
    DO I=1,IMITATION_TIME_SCALE    ! The dynamics of the particles goes on as the imitation parameters are measured
       CALL PositionUpdate
       VelFocal=VelFocal+VEHICLES(Indice)%VELOCITY
       VelTarget=VelTarget+VEHICLES(Indice+1)%VELOCITY
       IF(VEHICLES(Indice+1)%POSITION-VEHICLES(Indice)%POSITION-1<=VEHICLES(Indice)%MAXIMUM_VELOCITY)THEN
         ContactFrequency=ContactFrequency+1
       END IF
   END DO
   ContactFrequency=ContactFrequency/IMITATION_TIME_SCALE
   VelFocal=VelFocal/IMITATION_TIME_SCALE
   VelTarget=VelTarget/IMITATION_TIME_SCALE
  END SUBROUTINE ImitationDynamics  

! Here we have the selection of the particle, the calculation of the probability to change and the process of changing by imitation

  SUBROUTINE Imitation
    INTEGER  :: I, Indice
    REAL :: VelFocal, VelTarget, ContactFrequency
    DO I=1,INT(0.05*L)
       Indice=INT(RandomGenerator(ISEED)*N) + 1
       CALL ImitationDynamics(Indice, VelFocal, VelTarget, ContactFrequency)
       IF(ContactFrequency>= CONTACT_FRQUENCY_MINIMUM)THEN                  ! minimum connectivity criterion
          IF(RandomGenerator(ISEED)<(VelTarget/(VelFocal+VelTarget)))THEN    ! increasing in payoff criterion
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice+1)%MAXIMUM_VELOCITY
             VEHICLES(Indice)%P=VEHICLES(Indice+1)%P
             VEHICLES(Indice)%P_0=VEHICLES(Indice+1)%P_0
             VEHICLES(Indice)%P_S=VEHICLES(Indice+1)%P_S
             VEHICLES(Indice)%P_F=VEHICLES(Indice+1)%P_F
             VEHICLES(Indice)%P_ACC=VEHICLES(Indice+1)%P_ACC
             VEHICLES(Indice)%P_DCC=VEHICLES(Indice+1)%P_DCC
             VEHICLES(Indice)%X_T=VEHICLES(Indice+1)%X_T
             VEHICLES(Indice)%X_S=VEHICLES(Indice+1)%X_S
             VEHICLES(Indice)%SP=VEHICLES(Indice+1)%SP
          END IF
       END IF 
    END DO
   
  END SUBROUTINE Imitation
  
! Here I present the mutation algorithm  
! The changes must be done in the apropriated parameters of each model, otherwise many changes will not present measurable differences in  bahavior  
  SUBROUTINE Mutation
    INTEGER :: I,  Indice
    REAL :: Aux
    DO I=1, INT(.05*N)
      Indice = INT(RandomGenerator(ISEED)*N)+1
      Aux = RandomGenerator(ISEED)
      SELECT CASE(VEHICLES(Indice)%SP)
      CASE(1)  
        IF(Aux<0.25)THEN
           IF(VEHICLES(Indice)%P+0.01<=P_MAX_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P+0.01
           END IF
        ELSE IF(Aux<0.5)THEN
           IF(VEHICLES(Indice)%P-0.01>=P_MIN_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P-0.01
           END IF   
        ELSE IF(Aux<0.75)THEN   
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY+1<MAX_VEL_MAX_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY+1
           END IF
        ELSE    
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY-1>MAX_VEL_MIN_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY-1
           END IF
        END IF  
        
      CASE(2:3)  
      
        IF(Aux<0.16667)THEN
           IF(VEHICLES(Indice)%P+0.01<=P_MAX_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P+0.01
           END IF
        ELSE IF(Aux<0.33334)THEN
           IF(VEHICLES(Indice)%P-0.01>=P_MIN_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P-0.01
           END IF   
        ELSE IF(Aux<0.5)THEN   
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY+1<MAX_VEL_MAX_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY+1
           END IF
        ELSE  IF(Aux<0.66668)THEN  
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY-1>MAX_VEL_MIN_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY-1
           END IF
        ELSE IF(Aux<0.83335)THEN      
           IF(VEHICLES(Indice)%P_F+0.01<=P_F_MAX_LIMIT)THEN
             VEHICLES(Indice)%P_F=VEHICLES(Indice)%P_F+0.01
           END IF   
        ELSE   
         IF(VEHICLES(Indice)%P_F-0.01>=P_F_MIN_LIMIT)THEN
             VEHICLES(Indice)%P_F=VEHICLES(Indice)%P_F-0.01
           END IF   
        END IF  
        
      CASE(4)  
        IF(Aux<0.25)THEN
          IF(VEHICLES(Indice)%X_s+0.01<=X_s_MAX_LIMIT)THEN
             VEHICLES(Indice)%X_s=VEHICLES(Indice)%X_s+0.01
           END IF
        ELSE IF(Aux<0.5)THEN
           IF(VEHICLES(Indice)%X_s-0.01>=X_s_MIN_LIMIT)THEN
             VEHICLES(Indice)%X_s=VEHICLES(Indice)%X_s-0.01
           END IF   
        ELSE IF(Aux<0.75)THEN   
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY+1<MAX_VEL_MAX_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY+1
           END IF
        ELSE    
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY-1>MAX_VEL_MIN_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY-1
           END IF
        END IF  
      
      CASE(5)  
      
        IF(Aux<0.25)THEN
           IF(VEHICLES(Indice)%P_S+0.01<=P_S_MAX_LIMIT)THEN
             VEHICLES(Indice)%P_S=VEHICLES(Indice)%P_S+0.01
           END IF
        ELSE IF(Aux<0.5)THEN
           IF(VEHICLES(Indice)%P_S-0.01>=P_S_MIN_LIMIT)THEN
             VEHICLES(Indice)%P_S=VEHICLES(Indice)%P_S-0.01
           END IF   
        ELSE IF(Aux<0.75)THEN   
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY+1<MAX_VEL_MAX_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY+1
           END IF
        ELSE    
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY-1>MAX_VEL_MIN_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY-1
           END IF
        END IF  
            
      CASE(6)  
       
        IF(Aux<0.16667)THEN
           IF(VEHICLES(Indice)%P+0.01<=P_MAX_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P+0.01
           END IF
        ELSE IF(Aux<0.33334)THEN
           IF(VEHICLES(Indice)%P-0.01>=P_MIN_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P-0.01
           END IF   
        ELSE IF(Aux<0.5)THEN   
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY+1<MAX_VEL_MAX_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY+1
           END IF
        ELSE  IF(Aux<0.66668)THEN  
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY-1>MAX_VEL_MIN_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY-1
           END IF
        ELSE IF(Aux<0.83335)THEN      
           IF(VEHICLES(Indice)%P_0+0.01<=P_0_MAX_LIMIT)THEN
             VEHICLES(Indice)%P_0=VEHICLES(Indice)%P_0+0.01
           END IF   
        ELSE   
         IF(VEHICLES(Indice)%P_0-0.01>=P_0_MIN_LIMIT)THEN
             VEHICLES(Indice)%P_0=VEHICLES(Indice)%P_0-0.01
           END IF   
        END IF  
        
      CASE(7) 
      
        IF(Aux<0.1)THEN
           IF(VEHICLES(Indice)%P+0.01<=P_MAX_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P+0.01
           END IF
        ELSE IF(Aux<0.2)THEN
           IF(VEHICLES(Indice)%P-0.01>=P_MIN_LIMIT)THEN
             VEHICLES(Indice)%P=VEHICLES(Indice)%P-0.01
           END IF   
        ELSE IF(Aux<0.3)THEN   
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY+1<MAX_VEL_MAX_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY+1
           END IF
        ELSE  IF(Aux<0.4)THEN  
           IF(VEHICLES(Indice)%MAXIMUM_VELOCITY-1>MAX_VEL_MIN_LIMIT)THEN
             VEHICLES(Indice)%MAXIMUM_VELOCITY=VEHICLES(Indice)%MAXIMUM_VELOCITY-1
           END IF
        ELSE IF(Aux<0.5)THEN      
           IF(VEHICLES(Indice)%X_T+0.01<=X_T_MAX_LIMIT)THEN
             VEHICLES(Indice)%X_T=VEHICLES(Indice)%X_T+0.01
           END IF   
        ELSE IF(Aux<0.6)THEN  
          IF(VEHICLES(Indice)%X_T-0.01>=X_T_MIN_LIMIT)THEN
             VEHICLES(Indice)%X_T=VEHICLES(Indice)%X_T-0.01
           END IF   
        ELSE IF(Aux<0.7)THEN
           IF(VEHICLES(Indice)%P_ACC+0.01<=P_ACC_MAX_LIMIT)THEN
             VEHICLES(Indice)%P_ACC=VEHICLES(Indice)%P_ACC+0.01
           END IF
        ELSE IF(Aux<0.8)THEN
           IF(VEHICLES(Indice)%P_ACC-0.01>=P_ACC_MIN_LIMIT)THEN
             VEHICLES(Indice)%P_ACC=VEHICLES(Indice)%P_ACC-0.01
           END IF   
        ELSE IF(Aux<0.9)THEN
           IF(VEHICLES(Indice)%P_DCC+0.01<=P_DCC_MAX_LIMIT)THEN
             VEHICLES(Indice)%P_DCC=VEHICLES(Indice)%P_DCC+0.01
           END IF
        ELSE 
           IF(VEHICLES(Indice)%P_DCC-0.01>=P_DCC_MIN_LIMIT)THEN
             VEHICLES(Indice)%P_DCC=VEHICLES(Indice)%P_DCC-0.01
           END IF   
        END IF  
    
      CASE DEFAULT
        print*,"Invalid selection (mutation)"
      END SELECT
    END DO
  END SUBROUTINE Mutation
   
  
  
! *********************************************************************************************************************************************
! I begin now the measurements routines
!*********************************************************************************************************************************************  
  
  
  ! This routine prints the fundamental diagram   of the model  
  
  SUBROUTINE FundamentalDiagram
    INTEGER :: i, j, k, configurations
    REAL(KIND = RDP) :: velocity
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A23)")"FundamentalDiagram7.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    DO N=INT(0.01*L),INT(0.99*L),INT(0.01*L) ! Fundamental Diagram with 100 points 
      ALLOCATE(VEHICLES(-N : N + 1)) 
      velocity = 0.
      DO configurations = 1, NUMBER_OF_CONFIGURATIONS
        CALL MultiBehavioralRandomInitialConditions
        !CALL HeteroNaSchInitialConditions
        DO i = 1,TERMALIZATION_TIME    ! Termalization of the system
          CALL PositionUpdate
        END DO
        DO k = 1, SECOND_TIME_SCALE
          DO j = 1,FIRST_TIME_SCALE  
            CALL PositionUpdate
            DO i = 1, N
              velocity = velocity + VEHICLES(i)%VELOCITY
            END DO 
          END DO  
        END DO   
      END DO 
      ! the averaging over the velocity occurs in two steps because of float point errors 
      velocity = velocity/(N*SECOND_TIME_SCALE)
      velocity = velocity/(FIRST_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
      WRITE(1,*)REAL(N)/L,(REAL(N)/L)*velocity,velocity
      DEALLOCATE(VEHICLES)
      print *,REAL(N)/L,(REAL(N)/L)*velocity,velocity
    END DO  
    CLOSE(1)
  END SUBROUTINE FundamentalDiagram

  
  
  
 SUBROUTINE FundamentalDiagramLocalPlots
    INTEGER :: i, j, k, aux, configurations
    REAL(KIND = RDP) :: velocity
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A33)")"FundamentalDiagramLocalPlots1.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    DO N=INT(0.01*L),INT(0.99*L),INT(0.01*L) ! Fundamental Diagram with 100 points 
      ALLOCATE(VEHICLES(-N : N + 1)) 
      DO configurations = 1, NUMBER_OF_CONFIGURATIONS
        CALL MultiBehavioralRandomInitialConditions
        !CALL HeteroNaSchInitialConditions
        DO i = 1,TERMALIZATION_TIME    ! Termalization of the system
          CALL PositionUpdate
        END DO
        DO k = 1, SECOND_TIME_SCALE
          DO j = 1,FIRST_TIME_SCALE  
            CALL PositionUpdate
            aux = 0
            velocity = 0
            DO i = 1,N
              IF(MOD(VEHICLES(i)%POSITION,L)>INT(0.45*L) .AND.MOD(VEHICLES(i)%POSITION,L)<INT(0.55*L)) THEN
                velocity = velocity + 1.0*VEHICLES(i)%VELOCITY
                aux=aux+1
              END IF 
            END DO 
            IF(aux>0)THEN
               velocity=velocity/aux
               WRITE(1,*)1.0*aux/(0.1*L),velocity*aux/(0.1*L),velocity
            END IF 
          END DO            
        END DO   
      END DO 
      ! the averaging over the velocity occurs in two steps because of float point errors 
      !velocity = velocity/(N*SECOND_TIME_SCALE)
      !velocity = velocity/(FIRST_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
     ! WRITE(1,*)REAL(N)/L,(REAL(N)/L)*velocity,velocity
      DEALLOCATE(VEHICLES)
      print *,REAL(N)/L!,(REAL(N)/L)*velocity,velocity
    END DO  
    CLOSE(1)
 END SUBROUTINE FundamentalDiagramLocalPlots

  
  
  ! This routine is a visual aid. It prints a spatio-temporal pattern of the system in a matrix form that can be visualized in gnuplot with the following code
  ! reset
  ! set terminal qt 
  ! set pm3d map
  ! splot "SpatioTemporalMatrix.dat" matrix
  
  ! The data processing may present problems when the matrix is very big and I recomend a size of L=1.000 but just because my machine is not very powerful.  
  SUBROUTINE SpatioTemporalPattern
    INTEGER :: i, j
    INTEGER, DIMENSION(0:L-1) :: streetocupation
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A25)")"SpatioTemporalMatrix9.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    ! specify the density
    N = 0.9*L
    ALLOCATE(VEHICLES(-N : N + 1)) 
    CALL MultiBehavioralRandomInitialConditions
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
      CALL PositionUpdate
    END DO 
    CLOSE(1)
    DEALLOCATE(VEHICLES)
  END SUBROUTINE SpatioTemporalPattern

  
  ! This routine measures the distance headway  of the particles and exibts the results in a histogram
  
  SUBROUTINE DistanceHeadway
    INTEGER :: i, j , k, ii, Aux
    REAL(KIND=rdp), DIMENSION(0:L) :: DH
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A9)")"DHMP3.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    N=INT(0.16*L)                             ! Initialize the density of particles
    DH=0
    ALLOCATE(VEHICLES(1 : N + 1)) 
    DO ii = 1, NUMBER_OF_CONFIGURATIONS
      CALL MultiBehavioralRandomInitialConditions
      !CALL HeteroNaSchInitialConditions
      DO i = 1,TERMALIZATION_TIME
        CALL PositionUpdate
      END DO
      DO k = 1, SECOND_TIME_SCALE
        DO j = 1,FIRST_TIME_SCALE  
          CALL PositionUpdate
          DO i = 1, N
            Aux = VEHICLES(i+1)%POSITION-VEHICLES(i)%POSITION-1
            DH(Aux)=DH(Aux)+1
          END DO 
        END DO  
      END DO 
      print*,ii
    END DO
    DH=DH/(N*FIRST_TIME_SCALE)
    DH=DH/(SECOND_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
    DO i=0,L
      WRITE(1,*)i,DH(i)
    END DO
    DEALLOCATE(VEHICLES)  
    CLOSE(1)
  END SUBROUTINE DistanceHeadway

  SUBROUTINE DistanceHeadway2
    INTEGER :: i, j , k, ii, jj, Aux
    REAL(KIND=rdp), DIMENSION(0:L) :: DH
    CHARACTER(LEN=100) :: Arq1
    jj=0
    DO N=INT(0.10*L),INT(0.25*L),INT(0.01*L)
      DH=0
      ALLOCATE(VEHICLES(1 : N + 1)) 
      DO ii = 1, NUMBER_OF_CONFIGURATIONS
        CALL MultiBehavioralRandomInitialConditions
        DO i = 1,TERMALIZATION_TIME
          CALL PositionUpdate
        END DO
        DO k = 1, SECOND_TIME_SCALE
          DO j = 1,FIRST_TIME_SCALE  
            CALL PositionUpdate
            DO i = 1, N
              Aux = VEHICLES(i+1)%POSITION-VEHICLES(i)%POSITION-1
              DH(Aux)=DH(Aux)+1
            END DO 
          END DO  
        END DO 
        print*,ii
      END DO
      DH=DH/(N*FIRST_TIME_SCALE)
      DH=DH/(SECOND_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
      WRITE(Arq1,"(A5,I2.2,A4)")"7DHMP",jj,".dat"
      OPEN(1,FILE=TRIM(Arq1))   
      DO i=0,L
        WRITE(1,*)i,DH(i)
      END DO
      DEALLOCATE(VEHICLES)  
      jj=jj+1
      CLOSE(1)
    END DO  
  END SUBROUTINE DistanceHeadway2

  
  
  SUBROUTINE DesiredSpeedxHeadway
    INTEGER :: i, j , k, ii, Aux
    REAL(KIND=rdp), DIMENSION(0:L,2) :: DH
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A10)")"05DSxH.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    N=INT(0.05*L)                             ! Initialize the density of particles
    DH=0
    ALLOCATE(VEHICLES(1 : N + 1)) 
    DO ii = 1, NUMBER_OF_CONFIGURATIONS
      CALL MultiBehavioralRandomInitialConditions
      !CALL HeteroNaSchInitialConditions
      DO i = 1,TERMALIZATION_TIME
        CALL PositionUpdate
      END DO
      DO k = 1, SECOND_TIME_SCALE
        DO j = 1,FIRST_TIME_SCALE  
          CALL PositionUpdate
          DO i = 1, N
            Aux = VEHICLES(i+1)%POSITION-VEHICLES(i)%POSITION-1
            DH(Aux,1)=DH(Aux,1)+1
            DH(Aux,2)=DH(Aux,2)+VEHICLES(i)%VELOCITY
          END DO 
        END DO  
      END DO 
      print*,ii
    END DO
    DO i=0,L
       IF(INT(DH(i,1))==0)THEN
          DH(i,2)=0
       ELSE
          DH(i,2)=DH(i,2)/DH(i,1)
       END IF 
       DH(i,1)=DH(i,1)/(N*FIRST_TIME_SCALE)
       DH(i,1)=DH(i,1)/(SECOND_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
    END DO
    DO i=0,L
      WRITE(1,*)i,DH(i,1),DH(i,2)
    END DO
    DEALLOCATE(VEHICLES)  
    CLOSE(1)
  END SUBROUTINE DesiredSpeedxHeadway

  
  
  ! This routine measures the discretized time headway between two particles, given the that the velocity of the back particle is not null. The results is presented in a histogram. 
  
  SUBROUTINE TimeHeadway
    INTEGER :: i, j , k, ii, Aux
    REAL(KIND=rdp), DIMENSION(0:L) :: TH
    CHARACTER(LEN=100) :: Arq1
    WRITE(Arq1,"(A9)")"7TH05.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    N=INT(0.05*L)                             ! Initialize the density of particles
    TH=0
    ALLOCATE(VEHICLES(1 : N + 1)) 
    DO ii = 1, NUMBER_OF_CONFIGURATIONS
      CALL MultiBehavioralRandomInitialConditions
      DO i = 1,TERMALIZATION_TIME
        CALL PositionUpdate
      END DO
      DO k = 1, SECOND_TIME_SCALE
        DO j = 1,FIRST_TIME_SCALE  
          CALL PositionUpdate
          DO i = 1, N
            Aux = VEHICLES(i+1)%POSITION-VEHICLES(i)%POSITION-1
            IF (VEHICLES(i)%VELOCITY>0)THEN
              TH(Aux)=TH(Aux)+10.0*Aux/VEHICLES(i)%VELOCITY
            END IF 
          END DO 
        END DO  
      END DO 
      print*,ii
    END DO
    TH=TH/(N*FIRST_TIME_SCALE)
    TH=TH/(SECOND_TIME_SCALE*NUMBER_OF_CONFIGURATIONS)
    DO i=0,L
      WRITE(1,*)i,TH(i)
    END DO
    DEALLOCATE(VEHICLES)  
    CLOSE(1)
  END SUBROUTINE TimeHeadway
  
  
  ! This routine produces the trajectory of the heterogeneous NaSch model in the space of parameters (v_max,p). It is possible to mark a specific particle or plot the averages. I plot the averages in the following. 
    
  
  SUBROUTINE ParametersSpaceTrajectory
    INTEGER :: i, j 
    REAL(KIND=rdp), DIMENSION(2,TRAJECTORIES_SCALE) :: Aux
    REAL :: TempVel, TempP
    CHARACTER(LEN=100) :: Arq1
    N=INT(0.05*L)
    ALLOCATE(VEHICLES(1 : N+1)) 
    CALL HeteroNaSchInitialConditions
    Aux=0          
    DO j=1,TERMALIZATION_TIME
      CALL PositionUpdate
    END DO
    DO j=1,TRAJECTORIES_SCALE
        DO i=1,FIRST_TIME_SCALE  
          CALL PositionUpdate
          !CALL Mutation
          CALL Imitation
        END DO  
        TempVel=0
        TempP=0
        DO i=1, N
          TempVel=TempVel+VEHICLES(i)%MAXIMUM_VELOCITY
          TempP=TempP+VEHICLES(i)%P
        END DO 
        TempVel=TempVel/N
        TempP=TempP/N
        Aux(1,j)=TempVel
        Aux(2,j)=TempP
        print*,j
    END DO
    DEALLOCATE(VEHICLES)
    WRITE(Arq1,"(A28)")"ParameterSpaceTrajectory.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    DO j = 1, TIME_SERIES_SCALE
      WRITE(1,*)Aux(1,j),Aux(2,j)
    END DO   
    CLOSE(1)
  END SUBROUTINE ParametersSpaceTrajectory
  
  ! In systems with many variables, its analysis in the space of parameters is not feasible. Time series, on the other hand, offer a great simplification of this analysis at a cost of less information per graphic. The routine below make time-series of the parameters, the models fraction, or even the species (behaviors), under evolutionary dynamics.  It can produce data relative to individual runs, many runs togheter or their average.  In this case, I made a individual time serie of the models in the multibehavioral case    
  
  SUBROUTINE TimeSeries
    INTEGER :: i, j , k
    REAL(KIND=rdp), DIMENSION(7,TIME_SERIES_SCALE) :: Aux
    CHARACTER(LEN=100) :: Arq1
    N=INT(0.05*L)
    ALLOCATE(VEHICLES(1 : N+1)) 
    Aux=0
    DO k=1,NUMBER_OF_CONFIGURATIONS
      CALL MultiBehavioralRandomInitialConditions
      DO j=1,TERMALIZATION_TIME
        CALL PositionUpdate
      END DO
      TimeSerie: DO j=1,TIME_SERIES_SCALE
        DO i=1,FIRST_TIME_SCALE  
          CALL PositionUpdate
          !CALL Imitation
          !CALL Mutation
        END DO  
        DO i=1, N
          Aux(VEHICLES(i)%SP,j) = Aux(VEHICLES(i)%SP,j)+1.0/N
        END DO 
      END DO  TimeSerie
    END DO
    DEALLOCATE(VEHICLES)
    WRITE(Arq1,"(A14)")"TimeSeries.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    DO k = 1, TIME_SERIES_SCALE
      WRITE(1,*)k,Aux(1,k),Aux(2,k),Aux(3,k),Aux(4,k),Aux(5,k),Aux(6,k),Aux(7,k)
    END DO   
    CLOSE(1)
  END SUBROUTINE TimeSeries
  

  
  SUBROUTINE DensityVelocityFluxHeadwayTimeSeries
    INTEGER :: i, j , k, Aux
    REAL(KIND=rdp), DIMENSION(TIME_SERIES_SCALE) :: Auxp, Auxv, AuxJ, AuxH
    CHARACTER(LEN=100) :: Arq1
    N=INT(0.05*L)
    ALLOCATE(VEHICLES(1 : N+1)) 
    Auxp=0
    Auxv=0
    AuxJ=0
    AuxH=0
    DO k=1,NUMBER_OF_CONFIGURATIONS
      CALL MultiBehavioralRandomInitialConditions
      DO j=1,TERMALIZATION_TIME
        CALL PositionUpdate
      END DO
      TimeSerie: DO j=1,TIME_SERIES_SCALE
        DO i=1,FIRST_TIME_SCALE  
          CALL PositionUpdate
          !CALL Imitation
          !CALL Mutation
        END DO  
        Aux = 0
        DO i=1, N
          IF(MOD(VEHICLES(i)%POSITION,L)>INT(0.45*L).AND.MOD(VEHICLES(i)%POSITION,L)<INT(0.55*L))THEN 
            Aux = Aux + 1          
            Auxp(j) = Auxp(j) + 1.0/(0.1*L)
            Auxv(j) = Auxv(j) + VEHICLES(i)%VELOCITY
            AuxH(j) = AuxH(j) + (VEHICLES(i+1)%POSITION - VEHICLES(i)%POSITION - 1)
          END IF   
        END DO 
        IF(Aux==0)THEN
            Auxv(j) = 0
            AuxJ(j) = 0
            AuxH(j) = 0
        ELSE 
            Auxv(j) = Auxv(j)/Aux
            AuxH(j) = AuxH(j)/Aux
            AuxH(j) = Auxv(j)*Auxp(j)
        END IF 
      END DO  TimeSerie
    END DO
    DEALLOCATE(VEHICLES)
    WRITE(Arq1,"(A40)")"DensityVelocityFluxHeadwayTimeSeries.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    DO k = 1, TIME_SERIES_SCALE
      WRITE(1,*)k,Auxp(k),Auxv(k),AuxJ(k),AuxH(k)
    END DO   
    CLOSE(1)
  END SUBROUTINE DensityVelocityFluxHeadwayTimeSeries

  
  
   ! Now the measurements of the order parameter in the next subroutine. If this was only the measurements of real and unlimited 'big' integers we could have the average orientation of a of a normal vector as order parameter. For example we could fix v_max in the interval [1,99] and p in the interval [0.01,0.99] and measure the variable x_1 and x_2 which have 99 possible values each and a one to one correspondence to v_max and p.
  ! This way we could vary x_1 and x_2 in the interval [-1,1] in 99 points (49 negatives, 49 positives and 0) and we could define theta=artg(x_1,x_2) in such way that the representative vector in the parameter's space would be the normal vector (cos(theta),sen(theta)) and the order parameter would be identical to the used in the n-vector model.
  !Nevertheless, some parâmeters, such as beta, do not have unlimited variation, so we cannot  parametrize it to a variable  to all others as I have made above so I abandoned this idea. If you are treating a model or set of  parameters in which the idea presented above is feasible, The code below is easy to manipulate.
  !Let m_(j,beta) be the vector in the space of parameters related to the particle j following the algorithm beta, the order parameter theta is defined as teta = [sum  m_(j,beta)E(m)/(|E(m)|.|m_(j,beta)|)]/N, which maximum value is 1 which appears in homogeneous populations.
  !Aditionally, In this routine I also present Binder cumulant measurements for ensemble averages. 
  
  SUBROUTINE OrderParameter
    INTEGER :: i, j, k, ii
    REAL(KIND=rdp), DIMENSION(1:10) :: AverageVector
    REAL(KIND=rdp), DIMENSION(1:10, N) :: EffectiveVector
    REAL(KIND=rdp), DIMENSION(1:NUMBER_OF_CONFIGURATIONS) ::  Results
    REAL(KIND=rdp) :: Theta, ModEffectiveVector, ModAverageVector,VectorProduct, Cont,&
                      &Cont2, Cont4, BinderCumulant, ProbMut     
    CHARACTER(LEN=  100) :: Arq1
    
    WRITE(Arq1,"(A18)")"OrderParameter.dat"
    OPEN(1,FILE=TRIM(Arq1))   
    N=INT(0.05*L)
    DO ii=0,NUMBER_OF_CONFIGURATIONS
      ProbMut=1.0*ii/(ii+1)  ! probability to engange in the mutation process...1 is the fixed number of imitation process that can be changed to vary 
                                          ! tau_i/tau_m in a certain range. Here we go from no mutation to mutation 20 times more frequent then imitation
      Results=0
      EffectiveVector=0
      AverageVector=0
      N=INT(0.05*L)
      ALLOCATE(VEHICLES(1 : N+1)) 
      DO k=1,NUMBER_OF_CONFIGURATIONS
        CALL MultiBehavioralRandomInitialConditions
        DO j=1,TERMALIZATION_TIME
          CALL PositionUpdate
        END DO    
                                                                                  
        DO j=1,EVOLUTION_SCALE                               ! Evolutionary procosses
          CALL PositionUpdate
          IF(RandomGenerator(ISEED)<ProbMut)THEN ! Either mutation or Imitation happens
            CALL Mutation  
          ELSE   
            CALL Imitation
          END IF  
        END DO 
        
        !Measurements
        DO j=1,EVOLUTION_SCALE
          CALL PositionUpdate
          IF(RandomGenerator(ISEED)<ProbMut)THEN
            CALL Mutation  
          ELSE   
            CALL Imitation
          END IF  
          DO i=1,N            ! I have made the transformation to the auxiliar vetor in a separete routine
            CALL VectorBuilder(i)          
            AverageVector(1)=AverageVector(1)+1.0*TEMPVECTOR%MAXIMUM_VELOCITY
            EffectiveVector(1,i)=1.0*TEMPVECTOR%MAXIMUM_VELOCITY
            AverageVector(2)=AverageVector(2)+TEMPVECTOR%P
            EffectiveVector(2,i)=TEMPVECTOR%P
            AverageVector(3)=AverageVector(3)+TEMPVECTOR%P_F
            EffectiveVector(3,i)=TEMPVECTOR%P_F
            AverageVector(4)=AverageVector(4)+TEMPVECTOR%P_S
            EffectiveVector(4,i)=TEMPVECTOR%P_S
            AverageVector(5)=AverageVector(5)+TEMPVECTOR%P_0
            EffectiveVector(5,i)=TEMPVECTOR%P_0
            AverageVector(6)=AverageVector(6)+TEMPVECTOR%P_ACC
            EffectiveVector(6,i)=TEMPVECTOR%P_ACC
            AverageVector(7)=AverageVector(7)+TEMPVECTOR%P_DCC
            EffectiveVector(7,i)=TEMPVECTOR%P_DCC
            AverageVector(8)=AverageVector(8)+TEMPVECTOR%X_T
            EffectiveVector(8,i)=TEMPVECTOR%X_T
            AverageVector(9)=AverageVector(9)+TEMPVECTOR%X_S
            EffectiveVector(9,i)=TEMPVECTOR%X_S
            AverageVector(10)=AverageVector(10)+1.0*TEMPVECTOR%SP
            EffectiveVector(10,i)=TEMPVECTOR%SP
          END DO          
          AverageVector=AverageVector/N       
          Theta=0
          DO i=1,N !m_i*m_i and m_i*<m>
            ModEffectiveVector=EffectiveVector(1,i)*EffectiveVector(1,i)+EffectiveVector(2,i)*EffectiveVector(2,i)+&
            &EffectiveVector(3,i)*EffectiveVector(3,i)+EffectiveVector(4,i)*EffectiveVector(4,i)+&
            &EffectiveVector(5,i)*EffectiveVector(5,i)+EffectiveVector(6,i)*EffectiveVector(6,i)+&
            &EffectiveVector(7,i)*EffectiveVector(7,i)+EffectiveVector(8,i)*EffectiveVector(8,i)+&
            &EffectiveVector(9,i)*EffectiveVector(9,i)+EffectiveVector(10,i)*EffectiveVector(10,i)
            VectorProduct=EffectiveVector(1,i)*AverageVector(1)+EffectiveVector(2,i)*AverageVector(2)+&
            &EffectiveVector(3,i)*AverageVector(3)+EffectiveVector(4,i)*AverageVector(4)+&
            &EffectiveVector(5,i)*AverageVector(5)+EffectiveVector(6,i)*AverageVector(6)+&
            &EffectiveVector(7,i)*AverageVector(7)+EffectiveVector(8,i)*AverageVector(8)+&
            &EffectiveVector(9,i)*AverageVector(9)+EffectiveVector(10,i)*AverageVector(10)
            Theta=Theta+(VectorProduct/(SQRT(ModAverageVector*ModEffectiveVector))) 
          END DO        
          Theta=Theta/N
          Results(k)=Results(k)+Theta
        END DO    
        Results(k)=Results(k)/EVOLUTION_SCALE
        PRINT*,k,Results(k)
      END DO 
      DEALLOCATE(VEHICLES) 
      Cont=0
      Cont2=0
      Cont4=0
      DO k=1, NUMBER_OF_CONFIGURATIONS
        Cont=Cont+Results(k)
        Cont2=Cont2+Results(k)*Results(k)
        Cont4=Cont4+Results(k)*Results(k)*Results(k)*Results(k)
      END DO 
      Cont=Cont/NUMBER_OF_CONFIGURATIONS
      Cont2=Cont2/NUMBER_OF_CONFIGURATIONS
      Cont4=Cont4/NUMBER_OF_CONFIGURATIONS
      BinderCumulant=1-Cont4/(3*Cont2*Cont2)
      PRINT*,ii,Cont,BinderCumulant
      WRITE(1,*)ProbMut,Theta,BinderCumulant
    END DO
    CLOSE(1)
  END SUBROUTINE OrderParameter

  
  
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
  !CALL FundamentalDiagramLocalPlots
 ! CALL DistanceHeadway2
  !CALL ParametersSpaceTrajectory
  !CALL TimeHeadway
  !CALL DesiredSpeedxHeadway
  CALL DensityVelocityFluxHeadwayTimeSeries
END PROGRAM MAIN
