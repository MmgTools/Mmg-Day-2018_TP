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

// $ gcc firstSizeMap.c -o firstSizeMap -L ~/mmg/build/lib/ -lmmg2d -lm -I ~/mmg/build/include/
// with ~/mmg/build/lib/ that must be replace with the path toward the Mmg library
// ~/mmg/build/include/ the path toward the Mmg include directory.


/** Include the mmg2d library hader file */
#include "mmg/mmg2d/libmmg2d.h"

/**
 * \param x x coordinate of the mesh node
 * \param y y coordinate of the mesh node
 *
 * \return the wanted edge length at node
 *
 * Compute the scalar edge length that we want to apply to the node of
 * coor \a x \a y
 *
 */
double scalar_size(double x, double y) {
  double siz;

  // TO FILL WITH THE WANTED SCALAR SIZE.
  if ( x*x+y*y <= 3*3 ) {
    siz = 0.05;
  }
  else
    siz = 0.2;

  return siz;
}

/**
 * \param x x coordinate of the mesh node
 * \param y y coordinate of the mesh node
 * \param siz computed metric
 *
 * Compute the metric tensor to prescribe to the node of
 * coor \a x \a y
 *
 */
void tensor_size(double x, double y, double siz[3]) {

  // TO FILL WITH THE WANTED TENSORIAL SIZE.

  return;
}


int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  double          x,y,s;
  int             ier,k;
  char            *filein;

  fprintf(stdout,"  -- MMG2DLIB: mesh adaptation \n");

  if ( argc != 2 ) {
    printf(" Usage: %s filein\n",argv[0]);
    printf("        filein: name of your input mesh\n\n");
    printf(" This application saves the:\n");
    printf("        init.mesh and init.sol files : initial mesh and metric for"
           " vizualization purpose. \n");
    printf("        firstSizeMap.mesh file : the adapted mesh \n");

    return(1);
  }

  /* Name and path of the mesh file */
  filein = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filein == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filein,argv[1]);


  /** ------------------------------ STEP   I: load the data ---------------- */
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
  if ( MMG2D_loadMesh(mmgMesh,filein) != 1 ) {
    printf("Unable to read the input mesh: %s\n",filein);
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

  /** Computation of the size map over the mesh */
  for ( k=1; k<=np; ++k ) {
    /* Get mesh coordinates */
    if ( MMG2D_Get_vertex(mmgMesh,&x,&y,NULL,NULL,NULL) != 1 ) {
      printf("Unable to get the vertex coor at position: %d.\n",k);
      exit(EXIT_FAILURE);
    }

    /* Computation of the analytical function at mesh node */
    s = scalar_size(x,y);

    /* Give solution value at position k */
    if ( MMG2D_Set_scalarSol(mmgSol,s,k) != 1 ) {
      printf("Unable to set the metric at position %d.\n",k);
      exit(EXIT_FAILURE);
    }
  }

  /* Save mesh for vizualization purpose */
  if ( MMG2D_saveMesh(mmgMesh,"init.mesh") != 1 ) {
    printf("Unable to save initial mesh.\n");
    exit(EXIT_FAILURE);
  }

  /* Save solution for vizualization */
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"init.sol") != 1 ) {
    printf("Unable to save initial metric.\n");
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
  if ( MMG2D_saveMesh(mmgMesh,"firstSizeMap.mesh") != 1 ) {
    printf("Unable to save adapted mesh.\n");
    exit(EXIT_FAILURE);
  }

  /* Save solution for vizualization */
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"firstSizeMap.sol") != 1 ) {
    printf("Unable to save final metric.\n");
    exit(EXIT_FAILURE);
  }

  /** Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free ( filein );
  filein = NULL;

  return(0);
}
