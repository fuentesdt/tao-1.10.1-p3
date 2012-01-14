!
!  This program uses CPP for preprocessing, as indicated by the use of
!  TAO include files in the directories $TAO_DIR/include/finclude and
!  $PETSC_DIR/include/finclude.  This convention enables use of the CPP
!  preprocessor, which allows the use of the #include statements that
!  define TAO objects and variables.
!
!  Since one must be very careful to include each file no more than once
!  in a Fortran routine, application programmers must explicitly list
!  each file needed for the various TAO and PETSc components within their
!  program (unlike the C/C++ interface).
!
!  See the Fortran section of the PETSc users manual for details.
!
!  The following include statements are generally used in TAO programs:
!     tao_solver.h - TAO solvers
!     petscksp.h   - Krylov subspace methods
!     petscpc.h    - preconditioners
!     petscmat.h   - matrices
!     petscvec.h   - vectors
!     petsc.h      - basic PETSc routines
!  In addition, we need the following for use of distributed arrays
!     petscda.h    - distributed arrays (DAs)

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscda.h"
#include "finclude/tao_solver.h"

! Common blocks:
! In this example we use common blocks to store data needed by the 
! application-provided call-back routines, FormFunction(), FormGradient(),
! and FormHessian().  Note that we can store (pointers to) TAO objects
! within these common blocks.
!
! common /params/ - contains parameters for the global application
!   ecc         - eccentricity parameter
!   b           - .5 * upper limit for 2nd variable
!   hx, hy      - increment size in both directions
!   area        - area of triangles
!   mx, my      - discretization including boundaries
!
! common /mpi_info/
!   rank        - mpi rank (0...nproc-1)
!

      double precision ecc,b,hx,hy,area
      integer          mx,my
      double precision wq(0:1024), wl(0:1024)

      integer          rank


      common /params/  ecc,b,hx,hy,area,mx,my,wq,wl
      common /mpi_info/ rank

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
