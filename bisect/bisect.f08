! ---------------------------------------------------------------
! Simple program to demonstrate the bisection numerical technique
! for finding the solution to f(x) = 0 where f(x) is continuous
! on [a, b] where f(a) and f(b) have opposite signs
! ---------------------------------------------------------------

MODULE mod_bisect
IMPLICIT NONE
CONTAINS

  SUBROUTINE bisect_sub (func, endpts, maxitr, tol, output, err)
    ! Bisection subroutine
    ! Find a solution to f(x) = 0 given the continuous function
    ! f on the interval [a,b], where f(a) and f(b) have opposite
    ! signs.
    !
    ! Parameters
    ! ----------
    ! endPts: array of real, length 2 for [a,b]
    ! tol   : real; tolerance for sufficient accuracy
    ! maxitr: integer, max number of iterations
    ! output: real, variable to hold output
    ! func  : procedure pointer to function being evaluted
    ! err   : integer variable to hold error code
    real, dimension(2), intent (in) :: endpts
    real, intent (in) :: tol
    integer, intent (in) :: maxitr
    real, intent (out) :: output
    integer, intent (out) :: err

    ! define the interface for the passed function
    INTERFACE
      FUNCTION func (y)
        real, intent (in) :: y
        real :: func
      END FUNCTION func
    END INTERFACE

    ! local variables
    integer :: i = 1 ! iterator
    real :: a
    real :: b
    real :: midpoint ! unset to start
    real :: FA ! function F evaluated at a, unset to start
    real :: FP ! function F evaluated at midpoint

    a = endpts(1)
    b = endpts(2)

    ! Evaluate the function at a. While the iterator < maxitr,
    ! compute the midpoint between a and b and evaulate the function
    ! at that midpoint (result -> FP)
    FA = func(a)

    DO WHILE (i < maxitr)
      midpoint = a + (b - a) / 2.0
      FP = func(midpoint)
      WRITE (*, '(/ "Number iterations: ", I5)') i

      IF ( (FP == 0.0) .OR. ((b-a)/2.0 < tol) ) THEN
        output = midpoint
        err = 0
        RETURN
      END IF

      i = i + 1

      IF ( (FA * FP) > 0.0 ) THEN
        ! set the midpoint as the new lower boundary;
        ! we've narrowed the range of our solutions to one
        ! that is closer to b
        a = midpoint
      ELSE
        ! set midpoint as the new upper boundary; we've narrowed
        ! the range of solutions to one closer to a
        b = midpoint
      END IF

    END DO

    ! at this point, the procedure has exhausted maximum allowed
    ! iterations and should terminate unsuccessfully

    err = 1
    RETURN
    END SUBROUTINE bisect_sub

END MODULE mod_bisect

MODULE mod_example_funcs
CONTAINS
  ! sample functions to compute the root for using the
  ! bisection method

  FUNCTION func1 (x)
    ! f(x) = 3(x+1)(x-0.5)(x-1)

    ! parameters
    real, intent (in) :: x
    real :: func1

    func1 = 3.0 * (x + 1.0) * (x - 0.5) * (x - 1.0)

    return

  END FUNCTION func1
END MODULE mod_example_funcs


PROGRAM main

  USE mod_example_funcs
  USE mod_bisect
  IMPLICIT NONE

  ! inputs to bisect_sub
  real, dimension(2) :: interval
  real :: tol
  real :: output
  integer :: maxiter
  integer :: errcode

  ! define interface for generic functions f(x) which we
  ! can then pass to bisect()
  ABSTRACT INTERFACE
    FUNCTION func (x)
      ! define parameters
      real, intent (in) :: x
      ! define type of return value
      real :: func
    END FUNCTION func
  END INTERFACE

  ! use a procedure pointer to be able to pass different funcs around
  PROCEDURE (func), POINTER :: f_p => NULL ()

  ! point at function
  f_p => func1

  ! assign interval, tolerance, maximum iterations
  ! TODO CLI input
  interval = (/ -2.0, 1.5 /)
  tol = 10E-20
  maxiter = 50

  ! call subroutine to compute solution
  CALL bisect_sub(f_p, interval, maxiter, tol, output, errcode)

  IF (errcode /= 0) THEN
    WRITE (*, '(/ "ERROR COMPUTING ROOT: ", I5)') errcode
  ELSE
    WRITE(*, '(/ "computed roots: ", ES14.4 )' )  output
  END IF

END PROGRAM main
