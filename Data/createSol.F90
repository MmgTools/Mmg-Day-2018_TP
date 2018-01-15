!> @author
!> Cecile Dobrzynski, Charles Dapogny, Pascal Frey and Algiane Froehly
!> @brief
!>  Example for using mmg2dlib (basic use)


!> @param x x coordinate of the mesh node
!> @param y y coordinate of the mesh node
!>
!> @return the value of the analytic function at the node of coor \a x and
!> \a y of the mesh
!>
!> Compute the nodal value of an analytic function
!>
FUNCTION f(x, y) RESULT(sol)

  REAL(KIND=8), INTENT(IN) :: x,y
  REAL(KIND=8)             :: sol

  !! TO FILL WITH THE WANTED ANALYTIC FUNCTION.

END FUNCTION

!> @param x x coordinate of the mesh node
!> @param y y coordinate of the mesh node
!> @param hess the computed hessian value at coor \a x and \a y
!>
!> Compute the nodal value of the hessian of the analytic function.
!>
!> \remark: the Hessian is symetric definite positive so we store only h_11,h_12,h22.
!>
SUBROUTINE hessian(x, y, hess)

  REAL(KIND=8), INTENT(IN)                  :: x,y
  REAL(KIND=8), DIMENSION(0:2), INTENT(OUT) :: hess

  !! HESSIAN COMPUTATION

END SUBROUTINE

!> @param hess the hessian
!> @param lambda the computed hessian eigenvalues
!> @param vp eigen vectors
!>
!> Compute the eigenvalues and eigenvectors of the hessian \a hess.
!>
SUBROUTINE eigenvals(hess, lambda,vp)

  REAL(KIND=8), DIMENSION(0:2),     INTENT(IN)   :: hess
  REAL(KIND=8), DIMENSION(0:1)     ,INTENT(OUT ) :: lambda
  REAL(KIND=8), DIMENSION(0:1,0:1) ,INTENT(OUT ) :: vp

  REAL(KIND=8) :: sqDelta,dd,trm,vnorm
  REAL(KIND=8) :: epsd2 = 10e-30

  dd  = hess(0)-hess(2)
  trm = hess(0)+hess(2)
  sqDelta = sqrt(dd*dd + 4.0*hess(1)*hess(1))
  lambda(0) = 0.5*(trm - sqDelta)

  ! Case when m = lambda(0)*I
  IF ( sqDelta < epsd2) THEN
    lambda(1) = lambda(0)
    vp(0,0) = 1.0
    vp(0,1) = 0.0

    vp(1,0) = 0.0
    vp(1,1) = 1.0
    ier = 0
 ELSE
    vp(0,0) = hess(1)
    vp(0,1) = (lambda(0) - hess(0))
    vnorm = sqrt(vp(0,0)*vp(0,0) + vp(0,1)*vp(0,1))

    IF ( vnorm < epsd2 ) THEN
       vp(0,0) = (lambda(0) - hess(2))
       vp(0,1) = hess(1)
       vnorm = sqrt(vp(0,0)*vp(0,0) + vp(0,1)*vp(0,1))
    ENDIF

    IF ( vnorm <= epsd2 )THEN
       PRINT*, "Error in the eigenvector computation."
       CALL EXIT(200)
    ENDIF

    vnorm = 1.0/vnorm
    vp(0,0) = vp(0,0) * vnorm
    vp(0,1) = vp(0,1) * vnorm

    vp(1,0) = -vp(0,1)
    vp(1,1) = vp(0,0)

    lambda(1) = hess(0)*vp(1,0)*vp(1,0) + 2.0*hess(1)*vp(1,0)*vp(1,1) &
         + hess(2)*vp(1,1)*vp(1,1)
 ENDIF

END SUBROUTINE


!> @param x x coordinate of the mesh node
!> @param y y coordinate of the mesh node
!> @param epsilon maximal thershold for the error of interpolation
!>
!> @return the wanted edge length at node
!>
!> Compute the scalar edge length that we want to apply to the node of
!> coor \a x \a y
!>
FUNCTION scalar_size(x, y, epsilon) RESULT(siz)

  REAL(KIND=8), INTENT(IN) :: x,y,epsilon
  REAL(KIND=8)             :: siz

  !!  TO FILL WITH THE WANTED SCALAR SIZE.

END FUNCTION


!> @param x x coordinate of the mesh node
!> @param y y coordinate of the mesh node
!> @param epsilon maximal thershold for the error of interpolation
!> @param siz computed metric
!>
!> @return the wanted edge length at node
!>
!> Compute the metric tensor at the node of
!> coor \a x \a y
!>
SUBROUTINE tensor_size(x, y, epsilon, siz)

  REAL(KIND=8), INTENT(IN)                  :: x,y,epsilon
  REAL(KIND=8), DIMENSION(0:2), INTENT(OUT) :: siz

  !!  TO FILL WITH THE WANTED TENSORIAL SIZE.

END SUBROUTINE



PROGRAM main

  IMPLICIT NONE

!! To build this application:
! 1 you must clone the Mmg repo and build Mmg by hand:
! git clone https://github.com/MmgTools/ */

! $ gfortran createSol.F90 -o createSol -L ~/mmg/build/lib/ -lmmg2d -lm -I ~/mmg/build/include/
! with ~/mmg/build/lib/ that must be replace with the path toward the Mmg library
! ~/mmg/build/include/ the path toward the Mmg include directory.

!! Include the mmg2d library hader file
#include "mmg/mmg2d/libmmg2df.h"

  MMG5_DATA_PTR_T    :: mmgMesh
  MMG5_DATA_PTR_T    :: mmgSol
  INTEGER            :: ier,metType,argc,k
  CHARACTER(len=300) :: exec_name,filein,buffer

  INTEGER            :: ne,ned,np
  REAL(KIND=8)       :: scalar_size,f
  REAL(KIND=8)       :: scalarSol,tensorSol(0:2),epsilon,x,y

  PRINT*,"  -- MMG2DLIB: mesh adaptation"

  argc =  COMMAND_ARGUMENT_COUNT()
  CALL get_command_argument(0, exec_name)

  IF ( argc /=3 ) THEN
     PRINT*," Usage: ",TRIM(ADJUSTL(exec_name))," filein epsilon met_type"
     PRINT*,"        filein: name of your input mesh"
     PRINT*,"        epsilon: maximal wanted error of interpolation"
     PRINT*,"        met_type: wanted type for the metric: 0=scalar, 1=tensorial"
     PRINT*,""
     PRINT*," This application saves the"
     PRINT*,"        vizuSolution.mesh file: the initial mesh and analytic function"
     PRINT*,"        vizuMet.mesh file  : the initial mesh and computed metric"
     PRINT*,"        adaptedMesh.mesh file  : the final mesh"


     CALL EXIT(1)
  ENDIF

  ! Name and path of the mesh file
  CALL get_command_argument(1, filein)
  ! epsilon
  CALL get_command_argument(2, buffer)
  READ(buffer,*) epsilon
  ! metric type
  CALL get_command_argument(3, buffer)
  READ(buffer,*) metType

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
  CALL MMG2D_Get_meshSize(mmgMesh,np,ne,ned,ier)
  IF ( ier/=1 ) THEN
     PRINT*,"Unable to get the mesh size."
     CALL EXIT(103)
  ENDIF

  !! Set the solution/metric: sol applied on vertex entities, number of
  !! vertices=np, in this example, the solution/metric is scalar
  CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar,ier)
  IF ( ier/= 1 ) THEN
     PRINT*,"Unable to set the solution size."
     CALL EXIT(104)
  ENDIF

  !! Computation of the size map over the mesh
  DO k=1, np
     ! Get mesh coordinates
     CALL MMG2D_Get_vertex(mmgMesh,x,y,%val(0),%val(0),%val(0),ier)
     IF ( ier/= 1 ) THEN
        PRINT*,"Unable to get the vertex coor at position ",k
        CALL EXIT(104)
     ENDIF

     ! Computation of the analytical function at mesh node
     scalarSol = f(x,y)

     ! Give solution value at position k
     CALL MMG2D_Set_scalarSol(mmgSol,scalarSol,k,ier)
     IF ( ier/= 1 ) THEN
        PRINT*,"Unable to set the solution at position ",k
        CALL EXIT(105)
     ENDIF
  ENDDO

  ! Save mesh for vizualization purpose
  CALL  MMG2D_saveMesh(mmgMesh,"vizuSolution.mesh",LEN("vizuSolution.mesh"),ier)
  IF (  ier/= 1 ) THEN
     PRINT*,"Unable to save initial mesh for vizualization."
     CALL EXIT(105)
  ENDIF

  ! Save solution for vizualization
  CALL MMG2D_saveSol(mmgMesh,mmgSol,"vizuSolution.sol",LEN("vizuSolution.sol"),ier)
  IF ( ier /= 1 ) THEN
     PRINT*,"Unable to save analytic function for vizualization."
     CALL EXIT(106)
  ENDIF

  !! Compute the metric at mesh nodes
  IF (metType==0) THEN
     ! Scalar metric
     CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar,ier)
     IF  (ier /= 1) THEN
        PRINT*,"Unable to set the scalar metric size."
        CALL EXIT(107);
     ENDIF

    DO k=1, np
       ! Get mesh coordinates
       CALL MMG2D_Get_vertex(mmgMesh,x,y,%val(0),%val(0),%val(0),ier)
       IF ( ier/= 1 ) THEN
          PRINT*,"Unable to get the vertex coor at position: ",k
          CALL EXIT(108);
       ENDIF

       ! Computation of the metric at mesh node
       scalarSol = scalar_size(x,y,epsilon);

       ! Give metric value at position k
       CALL MMG2D_Set_scalarSol(mmgSol,scalarSol,k,ier)
       IF ( ier/= 1 ) THEN
          print*,"Unable to set the scalar metric at position: ",k
          CALL EXIT(109);
       ENDIF
    ENDDO

  ELSE
     ! Tensorial metric
     CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Tensor,ier)
     IF  (ier /= 1) THEN
        PRINT*,"Unable to set the tensorial metric size."
        CALL EXIT(110);
     ENDIF

     DO k=1, np
       ! Get mesh coordinates
       CALL MMG2D_Get_vertex(mmgMesh,x,y,%val(0),%val(0),%val(0),ier)
       IF ( ier/= 1 ) THEN
          PRINT*,"Unable to get the vertex coor at position: ",k
          CALL EXIT(111);
       ENDIF

       ! Computation of the metric at mesh node
       call tensor_size(x,y,epsilon,tensorSol);

       ! Give metric value at position k
       CALL MMG2D_Set_tensorSol(mmgSol,tensorSol(0),tensorSol(1),tensorSol(2),k,ier)
       IF ( ier/= 1 ) THEN
          print*,"Unable to set the tensorial metric at position: ",k
          CALL EXIT(112);
       ENDIF
    ENDDO
 ENDIF

 ! Save the initial mesh and associated metric
  CALL  MMG2D_saveMesh(mmgMesh,"vizuMet.mesh",LEN("vizuMet.mesh"),ier)
  IF (  ier/= 1 ) THEN
     PRINT*,"Unable to save initial mesh for metric vizualization."
     CALL EXIT(113)
  ENDIF

  ! Save solution for vizualization
  CALL MMG2D_saveSol(mmgMesh,mmgSol,"vizuMet.sol",LEN("vizuMet.sol"),ier)
  IF ( ier /= 1 ) THEN
     PRINT*,"Unable to save the metric for vizualization."
     CALL EXIT(114)
  ENDIF



  !> ------------------------------ STEP  II: call the remesher ----------------
  !! Remesher options
  ! Global hausdorff value (default value = 0.01) applied on the whole boundary
  CALL  MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hausd, 0.0001D0,ier)
  IF ( ier/= 1 ) THEN
     PRINT*,"Unable to set the hausdorff parameter."
     CALL EXIT(115)
  ENDIF
  ! Minimal edge size parameter
  CALL  MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmin, 0.00001D0,ier)
  IF ( ier/= 1 ) THEN
     PRINT*,"Unable to set the minimal edge size parameter."
     CALL EXIT(116)
  ENDIF

  !! remesh function
  ier = MMG5_SUCCESS !! CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)
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
  CALL MMG2D_saveMesh(mmgMesh,"adaptedMesh.mesh",LEN("adaptedMesh.mesh"),ier)
  IF ( ier /= 1 ) CALL EXIT(117)

  !> 2) Automatically save the solution
  CALL MMG2D_saveSol(mmgMesh,mmgSol,"adaptedMesh.sol",LEN("adaptedMesh.sol"),ier)
  IF ( ier /= 1 ) CALL EXIT(118)

  !> 3) Free the MMG2D5 structures
  CALL MMG2D_Free_all(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

END PROGRAM main
