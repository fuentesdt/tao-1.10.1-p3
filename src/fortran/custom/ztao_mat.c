/*$Id$*/

#include "private/fortranimpl.h"
/* #include "mat.h" */
#include "petscmat.h"

#ifdef PETSC_HAVE_FORTRAN_CAPS
#define matfdcoloringsetfunctiontao_    MATFDCOLORINGSETFUNCTIONTAO
#define matcreateada_                    MATCREATEADA

#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define matfdcoloringsetfunctiontao_    matfdcoloringsetfunctiontao
#define matcreateada_                    matcreateada
#endif

#include "tao_solver.h"

EXTERN_C_BEGIN

void PETSC_STDCALL matcreateada_(Mat mat, Vec D1, Vec D2,Mat *J,int *info)
{
  *info = MatCreateADA(
	(Mat)PetscToPointer((mat)),
	(Vec)PetscToPointer((D1)),
	(Vec)PetscToPointer((D2)),J);
}

/*
    MatFDColoringSetFunction sticks the Fortran function into the fortran_func_pointers;
    this function is then accessed by ourmatfdcoloringfunction()

    NOTE: FORTRAN USER CANNOT PUT IN A NEW J OR B currently.

    USER CAN HAVE ONLY ONE MatFDColoring in code Because there is no place to hang f7!
*/

static void (*f8)(TAO_SOLVER*,Vec*,Vec*,void*,int*);

static int ourmatfdcoloringfunctiontao(TAO_SOLVER tao,Vec x,Vec y,void *ctx)
{
  int info = 0;
  (*f8)(&tao,&x,&y,ctx,&info);
  return info;
}

void PETSC_STDCALL matfdcoloringsetfunctiontao_(MatFDColoring *fd,void (*f)(TAO_SOLVER*,Vec*,Vec*,void*,int*),
                                 void *ctx,int *info)
{
  f8 = f;
  *info = MatFDColoringSetFunction(*fd,(int (*)(void))ourmatfdcoloringfunctiontao,ctx);
}

EXTERN_C_END
