!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmg2dlib (basic use)


!> @param x x coordinate of the mesh node
!> @param y y coordinate of the mesh node
!>
!> @return the wanted edge length at node
!>
!> Compute the scalar edge length that we want to apply to the node of
!> coor \a x \a y
!>
FUNCTION scalar_size(x, y) RESULT(siz)

  REAL(KIND=8), INTENT(IN) :: x,y
  REAL(KIND=8)             :: siz

  !! TO FILL WITH THE WANTED SCALAR SIZE.

END FUNCTION



PROGRAM main

  IMPLICIT NONE

!! To build this application:
! 1 you must clone the Mmg repo and build Mmg by hand:
! git clone https://github.com/MmgTools/ */

! $ gfortran firstSizeMap.F90 -o firstSizeMap -L ~/mmg/build/lib/ -lmmg2d -lm -I ~/mmg/build/include/
! with ~/mmg/build/lib/ that must be replace with the path toward the Mmg library
! ~/mmg/build/include/ the path toward the Mmg include directory.

!! Include the mmg2d library hader file
#include "mmg/mmg2d/libmmg2df.h"

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol
  INTEGER            :: ier,argc,k
  CHARACTER(len=300) :: exec_name,filein

  INTEGER            :: ne,ned,np
  REAL(KIND=8)       :: scalar_size,x,y,s

  PRINT*,"  -- MMG2DLIB: mesh adaptation"

  argc =  COMMAND_ARGUMENT_COUNT();
  CALL get_command_argument(0, exec_name)

  IF ( argc /=1 ) THEN
     PRINT*," Usage: ",TRIM(ADJUSTL(exec_name))," filein"
     PRINT*,"        filein: name of your input mesh"
     PRINT*,""
     PRINT*," This application saves the"
     PRINT*,"        init.mesh and init.sol files : initial mesh and metric for&
          & vizualization purpose"
     PRINT*,"        firstSizeMap.mesh/sol file : the adapted mesh"


     CALL EXIT(1);
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, filein)

  !> ------------------------------ STEP   I: load the data  -------------------
  !! 1) Initialisation of mesh and sol structures
  !!   args of InitMesh:
  !! MMG5_ARG_start: we start to give the args of a variadic func
  !! MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
  !! mmgMesh: your MMG5_pMesh (that store your mesh)
  !! MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
  !! mmgSol: your MMG5_pSol (that store your metric) */

  mmgMesh = 0
  mmgSol  = 0

  CALL MMG2D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  !> 2) Build mesh in MMG5 format
  !! Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
  !! file formatted or manually set your mesh using the MMG2D_Set* functions

  !> with MMG2D_loadMesh function
  CALL MMG2D_loadMesh(mmgMesh,TRIM(ADJUSTL(filein)),&
       LEN(TRIM(ADJUSTL(filein))),ier)
  IF ( ier /= 1 )  CALL EXIT(102)

  !> 3) Build sol in MMG5 format
  !! Two solutions: just use the MMG2D_loadMet function that will read a .sol(b)
  !!    file formatted or manually set your sol using the MMG2D_Set* functions

  !! Get mesh size
  !! ne: Element number
  !! ned: Edge number
  !! np: Point number

  ! CALL THE Get_meshSize FUNCTION HERE

  !! Set the solution/metric: sol applied on vertex entities, number of
  !! vertices=np, in this example, the solution/metric is scalar

  ! CALL THE Set_solSize FUNCTION HERE

  !! Computation of the size map over the mesh
  DO k=1, np
     ! Get mesh coordinates

     ! CALL THE Get_vertex FUNCTION HERE

     ! Computation of the size at mesh node
     s = 0.; !scalar_size(x,y);

     ! Give solution value at position k

     ! CALL THE Set_scalarSol FUNCTION HERE

  ENDDO

  ! Save mesh for vizualization purpose
  CALL  MMG2D_saveMesh(mmgMesh,"init.mesh",LEN("init.mesh"),ier)
  IF (  ier/= 1 ) THEN
     PRINT*,"Unable to save initial mesh."
     CALL EXIT(105)
  ENDIF

  ! Save solution for vizualization
  CALL MMG2D_saveSol(mmgMesh,mmgSol,"init.sol",LEN("init.sol"),ier)
  IF ( ier /= 1 ) THEN
     PRINT*,"Unable to save initial metric."
     CALL EXIT(106)
  ENDIF

  !> ------------------------------ STEP  II: call the remesher ----------------
  !! Remesher options
  ! Global hausdorff value (default value = 0.01) applied on the whole boundary
  CALL  MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hausd, 0.0001D0,ier)
  IF ( ier/= 1 ) THEN
     PRINT*,"Unable to set the hausdorff parameter."
     CALL EXIT(107)
  ENDIF
  ! Minimal edge size parameter
  CALL  MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmin, 0.00001D0,ier)
  IF ( ier/= 1 ) THEN
     PRINT*,"Unable to set the minimal edge size parameter."
     CALL EXIT(107)
  ENDIF

  !! remesh function
  CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)
  IF ( ier == MMG5_STRONGFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH"
    STOP MMG5_STRONGFAILURE
  ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
     PRINT*,"BAD ENDING OF MMG2DLIB"
  ELSE
     PRINT*,"MMG2DLIB SUCCEED"
  ENDIF

  !> ------------------------------ STEP III: Save the result  -----------------
  !! get results
  !! Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
  !!    that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
  !!    using the MMG2D_getMesh/MMG2D_getSol functions

  !> 1) Automatically save the mesh
  CALL MMG2D_saveMesh(mmgMesh,"firstSizeMap.mesh",LEN("firstSizeMap.mesh"),ier)
  IF ( ier /= 1 ) CALL EXIT(108)

  !> 2) Automatically save the solution
  CALL MMG2D_saveSol(mmgMesh,mmgSol,"firstSizeMap.sol",LEN("firstSizeMap.sol"),ier)
  IF ( ier /= 1 ) CALL EXIT(109)

  !> 3) Free the MMG2D5 structures
  CALL MMG2D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main
