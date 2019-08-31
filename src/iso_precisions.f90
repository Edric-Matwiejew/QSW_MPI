!> @brief Imports numerical precisions.
!> @details Defines single precision (sp), double precision (dp) and quadrupule presision (qp) kind types using the 
!> instrinsic *iso_fortran_env* module types real32, real64 and real128 respectively.

module ISO_Precisions

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, qp => real128

end module ISO_Precisions


