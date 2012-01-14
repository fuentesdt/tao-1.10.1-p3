!
!  "$Id$";
!
!  Utility routine used to set some values into a Fortran common block

      subroutine TaoSetCommonBlock()
      implicit none

#include "finclude/petsc.h"
#include "finclude/tao_solver.h"

      TAO_NULL_CHARACTER = PETSC_NULL_CHARACTER
      TAO_NULL_INTEGER   = PETSC_NULL_INTEGER
      TAO_NULL_SCALAR    = PETSC_NULL_SCALAR
      TAO_NULL_DOUBLE    = PETSC_NULL_DOUBLE

      return
      end

