/**
 * Example of use of the mmg2d library (basic use of mesh adaptation)
 *
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/** To build this application: */
/* 1 you must clone the Mmg repo and build Mmg by hand:
git clone https://github.com/MmgTools/ */

// $ gcc createSol.c -L ~/mmg/build/lib/ -lmmg2d -lm -I ~/mmg/build/include/
// with ~/mmg/build/lib/ that must be replace with the path toward the Mmg library
// ~/mmg/build/include/ the path toward the Mmg include directory.


/** Include the mmg2d library hader file */
#include "mmg/mmg2d/libmmg2d.h"

/**
 * \param x x coordinate of the mesh node
 * \param y y coordinate of the mesh node
 *
 * \return the value of the analytic function at the node of coor \a x
 * and \a y of the mesh
 *
 * Compute the nodal value of an analytic function
 *
 */
double f(double x, double y) {
  double sol;

  // TO FILL WITH THE WANTED ANALYTIC FUNCTION

  /* Sinus function */
  sol = sin(x/2+y/2);

  return sol;
}

/**
 * \param x x coordinate of the mesh node
 * \param y y coordinate of the mesh node
 * \param hess the computed hessian value at coor \a x and \a y
 *
 * Compute the nodal value of the hessian of the analytic function.
 *
 * \remark: the Hessian is symetric definite positive so we store only h_11,h_12,h22
 *
 */
void hessian(double x, double y, double hess[3] ) {

  // HESSIAN COMPUTATION
  /* For the sinus function */
  // dx f = cos(x/2.+y/2.)/2.
  // dy f = cos(x/2.+y/2.)/2.
  hess[0] = -sin(x/2.+y/2.)/4.;
  hess[1] = -sin(x/2.+y/2.)/4.;
  hess[2] = -sin(x/2.+y/2.)/4.;

  return;
}

/**
 * \param hess the hessian
 * \param lambda the computed hessian eigenvalues
 * \param vp eigen vectors
 *
 * Compute the eigenvalues and eigenvectors of the hessian \a hess.
 *
 */
void eigenvals(double hess[3], double lambda[2],double vp[2][2]) {
  double       sqDelta,dd,trm,vnorm;
  const double epsd2 = 10e-30;

  dd  = hess[0]-hess[2];
  trm = hess[0]+hess[2];
  sqDelta = sqrt(dd*dd + 4.0*hess[1]*hess[1]);
  lambda[0] = 0.5*(trm - sqDelta);

  /* Case when m = lambda[0]*I */
  if ( sqDelta < epsd2) {
    lambda[1] = lambda[0];
    vp[0][0] = 1.0;
    vp[0][1] = 0.0;

    vp[1][0] = 0.0;
    vp[1][1] = 1.0;
    return;
  }
  vp[0][0] = hess[1];
  vp[0][1] = (lambda[0] - hess[0]);
  vnorm = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);

  if ( vnorm < epsd2 ) {
    vp[0][0] = (lambda[0] - hess[2]);
    vp[0][1] = hess[1];
    vnorm = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);
  }
  assert(vnorm > epsd2);

  vnorm = 1.0/vnorm;
  vp[0][0] *= vnorm;
  vp[0][1] *= vnorm;

  vp[1][0] = -vp[0][1];
  vp[1][1] = vp[0][0];

  lambda[1] = hess[0]*vp[1][0]*vp[1][0] + 2.0*hess[1]*vp[1][0]*vp[1][1]
    + hess[2]*vp[1][1]*vp[1][1];

  return;
}

/**
 * \param x x coordinate of the mesh node
 * \param y y coordinate of the mesh node
 * \param epsilon maximal thershold for the error of interpolation
 *
 * \return the wanted edge length at node
 *
 * Compute the scalar edge length that we want to apply to the node of
 * coor \a x \a y
 *
 */
double scalar_size(double x, double y, double epsilon) {
  double       lambda[2],vp[2][2],isqhmax,hess[3],siz;
  const double hmax = 10.;   // maximal edge size
  int          k;

  /* Hessian computation */
  hessian(x,y,hess);
  for ( k=0; k<3; ++k )
    hess[k] *= 2./(9.*epsilon);

  /* Eigenvalues/vectors */
  eigenvals(hess,lambda,vp);

  lambda[0] = fabs(lambda[0]);
  lambda[1] = fabs(lambda[1]);

  /* Size truncature */
  isqhmax = 1./(hmax*hmax);
  if ( isqhmax > lambda[0] ) lambda[0] = isqhmax;
  if ( isqhmax > lambda[1] ) lambda[1] = isqhmax;

  /* Metric computation */
  if ( lambda[0] >= lambda[1] ) siz = 1./sqrt(lambda[0]);
  else if ( lambda[1] > lambda[0] ) siz = 1./sqrt(lambda[1]);

  return siz;
}

/**
 * \param x x coordinate of the mesh node
 * \param y y coordinate of the mesh node
 * \param epsilon maximal thershold for the error of interpolation
 * \param siz computed metric
 *
 * Compute the metric tensor to prescribe to the node of
 * coor \a x \a y
 *
 */
void tensor_size(double x, double y, double epsilon, double siz[3]) {
  double       lambda[2],vp[2][2],isqhmax;
  const double hmax = 10.;   // maximal edge size
  int          k;

  /* Hessian computation */
  hessian(x,y,siz);
  for ( k=0; k<3; ++k )
    siz[k] *= 2./(9.*epsilon);

  /* Eigenvalues/vectors */
  eigenvals(siz,lambda,vp);

  lambda[0] = fabs(lambda[0]);
  lambda[1] = fabs(lambda[1]);

  /* Size truncature */
  isqhmax = 1./(hmax*hmax);
  if ( isqhmax > lambda[0] ) lambda[0] = isqhmax;
  if ( isqhmax > lambda[1] ) lambda[1] = isqhmax;

  /* Compute M = R D Rt */
  siz[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
  siz[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
  siz[2] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];

  return;
}


int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  double          epsilon,x,y,scalarSol,tensorSol[3];
  int             ier,metType,k;
  char            *filename, vizuSol[17],vizuMet[15],*saveMet;

  fprintf(stdout,"  -- TEST MMG2DLIB \n");

  if ( argc != 4 ) {
    printf(" Usage: %s filein epsilon met_type\n",argv[0]);
    printf("        filein: name of your input mesh\n");
    printf("        epsilon: maximal wanted error of interpolation\n");
    printf("        met_type: wanted type for the metric: 0=scalar, 1=tensorial\n\n");
    printf(" This application saves the:\n");
    printf("        vizuSolution.mesh file: the initial mesh and analytic function \n");
    printf("        vizuMet.mesh file  : the initial mesh and computed metric\n");
    printf("        adaptedMesh.mesh file  : the final mesh\n");

    return(1);
  }

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  epsilon = atof(argv[2]);

  metType = atoi(argv[3]);


  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh:
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */

  mmgMesh = NULL;
  mmgSol  = NULL;

  MMG2D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

 /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG2D_Set* functions */

  /** with MMG2D_loadMesh function */
  if ( MMG2D_loadMesh(mmgMesh,filename) != 1 ) {
    printf("Unable to read the input mesh: %s\n",filename);
    exit(EXIT_FAILURE);
  }

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG2D_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG2D_Set* functions */

  /** Get mesh size */
  int ne; // Element number
  int ned; // Edge number
  int np; // Point number
  if ( MMG2D_Get_meshSize(mmgMesh,&np,&ne,&ned) !=1 ) {
    printf("Unable to get the mesh size.\n");
    exit(EXIT_FAILURE);
  }

  /** Set the solution/metric: sol applied on vertex
      entities, number of vertices=np, in this example, the
      solution/metric is scalar */
  if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar) != 1 ) {
    printf("Unable to set the solution size.\n");
    exit(EXIT_FAILURE);
  }

  /** Computation of the analytical fonction over the mesh */
  for ( k=1; k<=np; ++k ) {
    /* Get mesh coordinates */
    if ( MMG2D_Get_vertex(mmgMesh,&x,&y,NULL,NULL,NULL) != 1 ) {
      printf("Unable to get the vertex coor at position: %d.\n",k);
      exit(EXIT_FAILURE);
    }

    /* Computation of the analytical function at mesh node */
    scalarSol = f(x,y);
    /* Give solution value at position k */
    if ( MMG2D_Set_scalarSol(mmgSol,scalarSol,k) != 1 ) {
      printf("Unable to set the solution at position %d.\n",k);
      exit(EXIT_FAILURE);
    }
  }

  /* Save mesh and analytic function for vizualisation */
  if ( MMG2D_saveMesh(mmgMesh,"vizuSolution.mesh") != 1 ) {
    printf("Unable to save the initial mesh for vizualization.\n");
    exit(EXIT_FAILURE);
  }
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"vizuSolution.sol") != 1 ) {
    printf("Unable to save the analytic function for vizualization.\n");
    exit(EXIT_FAILURE);
  }

  /** Compute the metric at mesh nodes */
  if ( metType==0 ) {
    /* Scalar metric */
    if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar) != 1 ) {
      printf("Unable to set the scalar metric size.\n");
      exit(EXIT_FAILURE);
    }

    for ( k=1; k<=np; ++k ) {
      /* Get mesh coordinates */
      if ( MMG2D_Get_vertex(mmgMesh,&x,&y,NULL,NULL,NULL) != 1 ) {
        printf("Unable to get the vertex coor at position: %d.\n",k);
        exit(EXIT_FAILURE);
      }

      /* Computation of the metric at mesh node */
      scalarSol = scalar_size(x,y,epsilon);

      /* Give metric value at position k */
      if ( MMG2D_Set_scalarSol(mmgSol,scalarSol,k) != 1 ) {
        printf("Unable to set the scalar metric at position %d.\n",k);
        exit(EXIT_FAILURE);
      }
    }
  }
  else {
    /* Tensorial metric */
    if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Tensor) != 1 ) {
      printf("Unable to set the tensorial metric size.\n");
      exit(EXIT_FAILURE);
    }

    for ( k=1; k<=np; ++k ) {
      /* Get mesh coordinates */
      if ( MMG2D_Get_vertex(mmgMesh,&x,&y,NULL,NULL,NULL) != 1 ) {
        printf("Unable to get the vertex coor at position: %d.\n",k);
        exit(EXIT_FAILURE);
      }

      /* Computation of the metric at mesh node */
      tensor_size(x,y,epsilon,tensorSol);

      /* Give metric value at position k */
      if ( MMG2D_Set_tensorSol(mmgSol,tensorSol[0],
                               tensorSol[1],tensorSol[2],k) != 1 ) {
        printf("Unable to set the tensorial metric at position %d.\n",k);
        exit(EXIT_FAILURE);
      }
    }
  }

  /* Save the initial mesh and the associated metric */
  if ( MMG2D_saveMesh(mmgMesh,"vizuMet.mesh") != 1 ) {
    printf("Unable to save the mesh for metric vizu.\n");
    exit(EXIT_FAILURE);
  }
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"vizuMet.sol") != 1 ) {
    printf("Unable to save the metric for vizu.\n");
    exit(EXIT_FAILURE);
  }

  /** ------------------------------ STEP   II: call the remesher ----------- */
  /** Remesher options */
  /* Global hausdorff value (default value = 0.01) applied on the whole boundary */
  if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hausd, 0.0001) != 1 ) {
    printf("Unable to set the hausdorff parameter.\n");
    exit(EXIT_FAILURE);
  }
  /* Minimal edge size */
  if ( MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmin, 0.00001) != 1 ) {
    printf("Unable to set the minimal edge size param.\n");
    exit(EXIT_FAILURE);
  }

  /** Remesher */
  ier = MMG2D_mmg2dlib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG2DLIB\n");


  /** ------------------------------ STEP   III: Save the result ------------ */
  /** Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
      that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
      using the MMG2D_getMesh/MMG2D_getSol functions */

  /** with MMG2D_saveMesh function */
  if ( MMG2D_saveMesh(mmgMesh,"adaptedMesh.mesh") != 1 ) {
    printf("Unable to save adapted mesh.\n");
    exit(EXIT_FAILURE);
  }

  /* Save solution for vizualization */
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"adaptedMesh.sol") != 1 ) {
    printf("Unable to save final metric.\n");
    exit(EXIT_FAILURE);
  }

  /** 3) Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free ( filename );
  filename = NULL;

  return(0);
}
