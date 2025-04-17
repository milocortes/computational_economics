module mod_julfort

    use iso_fortran_env

    implicit none

    public :: linint_Equi
    !##############################################################################
    !##############################################################################
    ! Interface declarations
    !##############################################################################
    !##############################################################################


    !##############################################################################
    ! INTERFACE assert_eq
    !
    ! Interface for equality assertions by assert_eqx functions.
    !##############################################################################
    interface assert_eq

    module procedure assert_eq2, assert_eq3, assert_eq4, assert_eq5, &
        assert_eqn

    end interface

    !##############################################################################
    ! INTERFACE linint_Equi
    !
    ! Interface for inverting gridpoints.
    !##############################################################################
    interface linint_Equi

    module procedure linint_Equi_1, linint_Equi_m

    end interface

contains




!############################################################################## 
!##############################################################################
! MODULE assertions
!##############################################################################
!##############################################################################
 
 
    !##############################################################################
    ! FUNCTION assert_eq2
    !
    ! Checks equality for two integers.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
    !     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
    !     Computing, 2nd edition. Cambridge: Cambridge University Press.
    !##############################################################################
function assert_eq2(n1, n2, string)
    
    
    !##### INPUT/OUTPUT VARIABLES #############################################
    
    ! integers to compare for equality
    integer, intent(in) :: n1, n2
    
    ! routine from which error should be thrown
    character(len=*), intent(in) :: string
    
    ! return value
    integer :: assert_eq2
    
    
    !##### ROUTINE CODE #######################################################
    
    ! if equality, set return value to n1
    if (n1 == n2)then
        assert_eq2 = n1
    
    ! else throw error message
    else
        call error(string, 'an assertion failed in assert_eq2')
    end if

end function assert_eq2


!##############################################################################
! FUNCTION assert_eq3
!
! Checks equality for three integers.
!
! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
!     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
!     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
!     Computing, 2nd edition. Cambridge: Cambridge University Press.
!##############################################################################
function assert_eq3(n1, n2, n3, string)


    !##### INPUT/OUTPUT VARIABLES #############################################
    
    ! integers to compare for equality
    integer, intent(in) :: n1, n2, n3
    
    ! routine from which error should be thrown
    character(len=*), intent(in) :: string
    
    ! return value
    integer :: assert_eq3
    
    
    !##### ROUTINE CODE #######################################################
    
    ! if equality, set return value to n1
    if (n1 == n2 .and. n2 == n3)then
        assert_eq3 = n1
    
    ! else throw error message
    else
        call error(string, 'an assertion failed in assert_eq3')
    end if

end function assert_eq3


!##############################################################################
! FUNCTION assert_eq4
!
! Checks equality for four integers.
!
! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
!     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
!     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
!     Computing, 2nd edition. Cambridge: Cambridge University Press.
!##############################################################################
function assert_eq4(n1, n2, n3, n4, string)


    !##### INPUT/OUTPUT VARIABLES #############################################
    
    ! integers to compare for equality
    integer, intent(in) :: n1, n2, n3, n4
    
    ! routine from which error should be thrown
    character(len=*), intent(in) :: string
    
    ! return value
    integer :: assert_eq4
    
    
    !##### ROUTINE CODE #######################################################
    
    ! if equality, set return value to n1
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4)then
        assert_eq4 = n1
    
    ! else throw error message
    else
        call error(string, 'an assertion failed in assert_eq4')
    end if

end function assert_eq4


!##############################################################################
! FUNCTION assert_eq5
!
! Checks equality for five integers.
!
! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
!     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
!     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
!     Computing, 2nd edition. Cambridge: Cambridge University Press.
!##############################################################################
function assert_eq5(n1, n2, n3, n4, n5, string)


    !##### INPUT/OUTPUT VARIABLES #############################################
    
    ! integers to compare for equality
    integer, intent(in) :: n1, n2, n3, n4, n5
    
    ! routine from which error should be thrown
    character(len=*), intent(in) :: string
    
    ! return value
    integer :: assert_eq5
    
    
    !##### ROUTINE CODE #######################################################
    
    ! if equality, set return value to n1
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4 .and. n4 == n5)then
        assert_eq5 = n1
    
    ! else throw error message
    else
        call error(string, 'an assertion failed in assert_eq5')
    end if

end function assert_eq5


!##############################################################################
! FUNCTION assert_eqn
!
! Checks equality for n integers.
!
! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
!     Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P. (1992). 
!     Numerical Recipes in Fortran 90: The Art of Parallel Scientific
!     Computing, 2nd edition. Cambridge: Cambridge University Press.
!##############################################################################
function assert_eqn(nn, string)


    !##### INPUT/OUTPUT VARIABLES #############################################
    
    ! integers to compare for equality
    integer, intent(in) :: nn(:)
    
    ! routine from which error should be thrown
    character(len=*), intent(in) :: string
    
    ! return value
    integer :: assert_eqn
    
    
    !##### ROUTINE CODE #######################################################
    
    ! if equality, set return value to n1
    if (all(nn(2:) == nn(1)))then
        assert_eqn = nn(1)
    
    ! else throw error message
    else
        call error(string, 'an assertion failed in assert_eqn')
    end if

end function assert_eqn

    
    !##############################################################################
    ! FUNCTION grid_Inv_Equi_m
    !
    ! Calculates inverse of gridpoints of an equidistant grid.
    !##############################################################################
function grid_Inv_Equi_m(x, left, right, n)
    
    
    !##### INPUT/OUTPUT VARIABLES #############################################
 
    ! point that shall be calculated
    real*8, intent(in) :: x(:)
 
    ! left and right interval point
    real*8, intent(in) :: left, right
 
    ! last grid point: 0,1,...,n
    integer, intent(in) :: n     
 
    ! value of the inverse of the gridpoint x \in [left, right]
    real*8 :: grid_Inv_Equi_m(size(x, 1))
 
 
    !##### OTHER VARIABLES ####################################################
 
    real*8 :: h
 
 
    !##### ROUTINE CODE #######################################################
 
    ! check for left <= right
    if(left >= right)call error('grid_Inv_Equi', &
        'left interval point greater than right point')
 
    ! calculate distance between grid points
    h = (right-left)/n

    ! calculate grid value
    grid_Inv_Equi_m = (x-left)/h

end function grid_Inv_Equi_m

    !##############################################################################
    ! FUNCTION grid_Inv_Equi_1
    !
    ! Calculates inverse of gridpoints of an equidistant grid.
    !##############################################################################
function grid_Inv_Equi_1(x, left, right, n)
    
    
    !##### INPUT/OUTPUT VARIABLES #############################################
 
    ! point that shall be calculated
    real*8, intent(in) :: x
 
    ! left and right interval point
    real*8, intent(in) :: left, right
 
    ! last grid point: 0,1,...,n
    !integer, intent(in) :: n
    integer(int64), intent(in) :: n

    ! value of the inverse of the gridpoint x \in [left, right]
    real*8 :: grid_Inv_Equi_1
 
 
    !##### OTHER VARIABLES ####################################################
 
    real*8 :: h
 
 
    !##### ROUTINE CODE #######################################################
 
    ! check for left <= right
    if(left >= right)call error('grid_Inv_Equi', &
        'left interval point greater than right point')
 
    ! calculate distance between grid points
    h = (right-left)/n
 
    ! calculate grid value
    grid_Inv_Equi_1 = (x-left)/h

end function grid_Inv_Equi_1

    !##############################################################################
    ! subroutine linint_Equi_1
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
subroutine linint_Equi_1(x, left, right, n, il, ir, phi)
    
    
    !##### INPUT/OUTPUT VARIABLES #############################################
 
    ! point that shall be calculated
    real*8, intent(in) :: x
 
    ! left and right interval point
    real*8, intent(in) :: left, right
 
    ! last grid point: 0,1,...,n
    !integer, intent(in) :: n
    integer(int64), intent(in) :: n

    ! left interpolation point
    !integer, intent(out) :: il
    integer(int64), intent(out) :: il(1)

    ! right interpolation point
    !integer, intent(out) :: ir
    integer(int64), intent(out) :: ir(1)

    ! interpolation fraction
    real*8, intent(out) :: phi(1)
 
    !##### OTHER VARIABLES ####################################################
 
    real*8 :: h, xinv, xl, xr
 
 
    !##### ROUTINE CODE #######################################################
         
    ! invert the grid to get point
    xinv = grid_Inv_Equi_1(min(max(x, left), right), left, right, n)
 
    ! get left and right gridpoint
    il = min(max(floor(xinv), 0), n-1)
    ir = il+1

    ! determine left and right gridpoint
    h  = (right-left)/n
    xl = h*dble(il(1))+left
    xr = h*dble(ir(1))+left

    ! get share on the left point
    phi = (xr-x)/(xr-xl)

end subroutine linint_Equi_1


!##############################################################################
! subroutine linint_Equi_m
!
! Calculates linear interpolant on an equidistant grid.
!##############################################################################
subroutine linint_Equi_m(x, left, right, n, il, ir, phi)


    !##### INPUT/OUTPUT VARIABLES #############################################
 
    ! point that shall be calculated
    real*8, intent(in) :: x(:)
 
    ! left and right interval point
    real*8, intent(in) :: left, right
 
    ! last grid point: 0,1,...,n
    integer, intent(in) :: n
 
    ! left interpolation point
    integer, intent(out) :: il(:)
    
    ! right interpolation point
    integer, intent(out) :: ir(:)

    ! interpolation fraction
    real*8, intent(out) :: phi(:)
 
    !##### OTHER VARIABLES ####################################################
 
    integer :: m
    real*8 :: h, xinv(size(x, 1)), xl(size(x, 1)), xr(size(x, 1))
 
 
    !##### ROUTINE CODE #######################################################

    ! check for sizes
    m = assert_eq(size(x, 1), size(il, 1), size(ir, 1), size(phi, 1), 'linint_Equi')
    m = m
         
    ! invert the grid to get point
    xinv = grid_Inv_Equi_m(min(max(x, left), right), left, right, n)
 
    ! get left and right gridpoint
    il = min(max(floor(xinv), 0), n-1)
    ir = il+1

    ! determine left and right gridpoint
    h  = (right-left)/n
    xl = h*dble(il)+left
    xr = h*dble(ir)+left

    ! get share on the left point
    phi = (xr-x)/(xr-xl)

end subroutine linint_Equi_m

    !##############################################################################
    ! SUBROUTINE grid_Cons_Equi
    !
    ! Constructs a whole equidistant grid on [left,right].
    !##############################################################################
subroutine grid_Cons_Equi(a, n, left, right)
     
     
    !##### INPUT/OUTPUT VARIABLES #############################################
    integer(int64), intent(in) :: n

    ! the array to fill
    real*8, intent(out) :: a(n)
 
    ! left and right interval point
    real*8, intent(in) :: left, right
 
 
    !##### OTHER VARIABLES ####################################################
 
    real*8 :: h
    integer :: j, n_new
 
 
    !##### ROUTINE CODE #######################################################
    
    ! get number of grid points
    n_new = size(a, 1)-1
 
    ! check for left <= right
    if(left >= right)call error('grid_Cons_Equi', &
        'left interval point greater than right point')
 
    ! calculate distance between grid points
    h = (right-left)/dble(n_new)
 
    ! calculate grid value
    a = h*(/(dble(j), j=0,n_new)/)+left
 
end subroutine grid_Cons_Equi
 

subroutine grid_Cons_Grow(a, n, left, right, growth)
    ! arguments
    integer(int64), intent(in) :: n
    ! left and right interval point
    real*8, intent(in) :: left, right, growth
     
    !integer(int64), intent(in) :: n
    real*8, intent(out) ::  a(n)

    ! the growth rate of the grid
    !real*8 :: growth

    !##### OTHER VARIABLES ####################################################
     
    real*8 :: h
    integer :: j, n_new

    !##### ROUTINE CODE #######################################################
        
        ! get number of grid points
    n_new = size(a, 1)-1

    ! check for left <= right
    if(left >= right)call error('grid_Cons_Grow', &
        'left interval point greater than right point')

    ! check for growth
    if(growth <= 0d0)call error('grid_Cons_Grow', &
        'growth rate must be greater than zero')

    ! calculate factor
    h = (right-left)/((1d0+growth)**n_new-1d0)
 
    ! calculate grid value
    a = h*((1d0+growth)**(/(dble(j), j=0,n_new)/)-1d0)+left
end subroutine grid_Cons_Grow


   
!##############################################################################
! subroutine linint_Grow_1
!
! Calculates linear interpolant on a growing grid.
!##############################################################################
subroutine linint_Grow_1(x, left, right, growth, n, il, ir, phi)
    
    !##### INPUT/OUTPUT VARIABLES #############################################
 
    ! point that shall be calculated
    real*8, intent(in) :: x
 
    ! left and right interval point
    real*8, intent(in) :: left, right, growth
 
    ! last grid point: 0,1,...,n
    !integer, intent(in) :: n
    integer(int64), intent(in) :: n

    ! left interpolation point
    !integer, intent(out) :: il
    integer(int64), intent(out) :: il(1)

    ! right interpolation point
    !integer, intent(out) :: ir
    integer(int64), intent(out) :: ir(1)

    ! interpolation fraction
    real*8, intent(out) :: phi(1)

    !##### OTHER VARIABLES ####################################################
 
    real*8 :: h, xinv, xl, xr
 
 
    !##### ROUTINE CODE #######################################################
 
    ! check for left <= right
    if(left >= right)call error('linint_Grow', &
        'left interval point greater than right point')
 
    ! invert the grid to get point
    xinv = grid_Inv_Grow_1(min(max(x, left), right), left, right, growth, n)
 
    ! get left and right gridpoint
    il = min(max(floor(xinv), 0), n-1)
    ir = il+1

    ! determine left and right gridpoint
    h = (right-left)/((1+growth)**n-1)     
    xl = h*((1+growth)**dble(il(1))-1d0)+left
    xr = h*((1+growth)**dble(ir(1))-1d0)+left

    ! get share on the left point
    phi = (xr-x)/(xr-xl)

end subroutine linint_Grow_1

    
    !##############################################################################
    ! FUNCTION grid_Inv_Grow_1
    !
    ! Calculates inverse of gridpoints of a growing grid.
    !##############################################################################
function grid_Inv_Grow_1(x, left, right, growth, n)
    
    
    !##### INPUT/OUTPUT VARIABLES #############################################
 
    ! point that shall be calculated
    real*8, intent(in) :: x
 
    ! left and right interval point
    real*8, intent(in) :: left, right

    ! growth rate
    real*8, intent(in) :: growth
 
    ! last grid point: 0,1,...,n
    !integer, intent(in) :: n             
    integer(int64), intent(in) :: n

    ! value of the inverse of the gridpoint x \in [left, right]
    real*8 :: grid_Inv_Grow_1
 
 
    !##### OTHER VARIABLES ####################################################
 
    real*8 :: h
 
 
    !##### ROUTINE CODE #######################################################
 
    ! check for left <= right
    if(left >= right)call error('grid_Inv_Grow', &
        'left interval point greater than right point')

    ! check for growth
    if(growth <= 0d0)call error('grid_Inv_Grow', &
        'growth rate must be greater than zero')
 
    ! calculate factor
    h = (right-left)/((1+growth)**n-1d0)
 
    ! calculate grid value
    grid_Inv_Grow_1 = log((x-left)/h+1d0)/log(1d0+growth)

end function grid_Inv_Grow_1


end module

module mod_julfort_wrapper

    use iso_c_binding

    use :: mod_julfort

    implicit none

contains


    subroutine grid_Cons_Equi_wrapper(a, n, left, right) bind(C, name="grid_Cons_Equi")
        ! arguments
        real(c_double), intent(in) :: left, right
        integer(c_long), intent(in) :: n
        real(c_double), intent(out) :: a(n)

        call grid_Cons_Equi(a, n, left, right)

    end subroutine grid_Cons_Equi_wrapper

    subroutine linint_Equi_wrapper(x, left, right, n, il, ir, phi) bind(C, name="linint_Equi")
 
        real(c_double), intent(in) :: x, left, right
        integer(c_long), intent(in) :: n
        integer(c_long), intent(out) :: il(1)
        integer(c_long), intent(out) :: ir(1)
        real(c_double), intent(out) :: phi(1)

        call linint_Equi_1(x, left, right, n, il, ir, phi)

    end subroutine linint_Equi_wrapper

    subroutine grid_Cons_Grow_wrapper(a, n, left, right, growth) bind(C, name="grid_Cons_Grow")
        ! arguments
        real(c_double), intent(in) :: left, right, growth
        integer(c_long), intent(in) :: n
        real(c_double), intent(out) :: a(n)

        call grid_Cons_Grow(a, n, left, right, growth)

    end subroutine grid_Cons_Grow_wrapper

    subroutine linint_Grow_wrapper(x, left, right, growth, n, il, ir, phi) bind(C, name="linint_Grow")
 
        real(c_double), intent(in) :: x, left, right, growth
        integer(c_long), intent(in) :: n
        integer(c_long), intent(out) :: il(1)
        integer(c_long), intent(out) :: ir(1)
        real(c_double), intent(out) :: phi(1)

        call linint_Grow_1(x, left, right, growth, n, il, ir, phi)

    end subroutine linint_Grow_wrapper

end module

