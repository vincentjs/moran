var tipuesearch = {"pages":[{"text":"peano A hodgepodge of mathematical routines written in modern Fortran. It provides functions for statistics, linear algebra, interpolation, search and sorts, root finding, fourier transforms, and differential equations. Some routines have been parallelized. This project's codename refers to Guiseppe Peano , a prolific late 19th century mathematician who contributed much to set theory and mathematical logic, as well as the fields of calculus and differential equations. In perhaps his most famous contribution, Peano introduces the definition of the natural numbers in terms of sets, resulting in a set of axioms that have become a cornerstone in modern set theory and mathematical induction. This project is licensed under the MIT License, and can be obtained at the project's github page, https://github.com/vincentjs/peano. Developer Info Vincent San Miguel","tags":"home","loc":"index.html","title":" peano "},{"text":"Procedures Procedure Location Procedure Type Description assert_same_bounds assert Subroutine Asserts that input vector /(/mathbf{x}/) has the same lower and upper bounds as input vector /(/mathbf{y}/) assert_same_rank assert Subroutine Asserts that input vector /(/mathbf{x}/) has the same rank as input vector /(/mathbf{y}/) assert_x_is_ge_y assert Interface assert_x_is_ge_y_int assert Subroutine [Integer] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) assert_x_is_ge_y_real assert Subroutine [Real] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) assert_x_is_gt_y assert Interface assert_x_is_gt_y_int assert Subroutine Asserts that input argument /(x/) is greater than input argument /(y/) assert_x_is_gt_y_real assert Subroutine Asserts that input argument /(x/) is greater than input argument /(y/) bl_growth airfoil_bl_growth Subroutine calc_bl_growth airfoil_bl_growth Subroutine Calculates boundary layer growth on an airfoil (modeled as an ellipse), beginning at a stagnation point. Uses Thwaite's method for the laminar flow region, Michel's method to correct transition, and Head's method for the turbulent flow region. ve airfoil_bl_growth Function VE of the ellipse,\n  v_e = \\left(1 + \\tau\\right) \\cdot \\sqrt{\\frac{1 - x_i&#94;2}{1-\\left(1-\\tau&#94;2\\right)x_i&#94;2}}  x airfoil_bl_growth Function X-coordinates of the ellipse,\n  x_i = -\\cos\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right)  y airfoil_bl_growth Function Y-coordinates of the ellipse,\n  y_i = \\sin\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right) ","tags":"list procedures","loc":"lists/procedures.html","title":"\nAll Procedures – peano\n"},{"text":"Source Files File Description airfoilBlGrowth.f90 assert.f90 constants.f90 error.f90 precision.f90","tags":"list files","loc":"lists/files.html","title":"\nAll Files – peano\n"},{"text":"Modules Module Source File Description airfoil_bl_growth airfoilBlGrowth.f90 Calculates the boundary-layer growth over an elliptical airfoil assuming a non-viscous fluid. For these geometries, a potential-flow closed-form solution can be obtained. We assume that the boundary layer growth begins at a stagnation point. assert assert.f90 Development tool used to ensure that predicates are as expected before proceeding. If an assertion fails, the error message will be logged and the program terminated. constants constants.f90 Provides constants used within this library. error error.f90 Development tool used to add and print error messages, as well as handle logic flow in the event an error is thrown. precision precision.f90 Provides kind attributes for setting machine-compiler-independent precision for real numbers. Here, we define \"single precision\" to mean 32-bit precision, \"double precision\" to mean 64-bit precision, and \"quadruple precision\" to mean 128-bit precision. This ensures that precision is preserved regardless of the compiler or computer architecture.","tags":"list modules","loc":"lists/modules.html","title":"\nAll Modules – peano\n"},{"text":"constants.f90 Source File Source File constants.f90 Modules constants All Source Files airfoilBlGrowth.f90 assert.f90 constants.f90 error.f90 precision.f90 module constants !! Provides constants used within this library. use precision implicit none real , parameter :: PI = 3.1415926535898 contains end module constants © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"sourcefile/constants.f90.html","title":"constants.f90 – peano"},{"text":"precision.f90 Source File Source File precision.f90 Modules precision All Source Files airfoilBlGrowth.f90 assert.f90 constants.f90 error.f90 precision.f90 module precision !! Provides kind attributes for setting machine-compiler-independent precision for real numbers. Here, we define \"single precision\" to mean 32-bit precision, \"double precision\" to mean 64-bit precision, and \"quadruple precision\" to mean 128-bit precision. This ensures that precision is preserved regardless of the compiler or computer architecture. use ISO_FORTRAN_ENV implicit none integer , parameter :: sp = REAL32 !selected_real_kind(6, 37) !! Single precision: 32-bit real integer , parameter :: dp = REAL64 !selected_real_kind(15, 307) !! Double precision: 64-bit real integer , parameter :: qp = REAL128 !selected_real_kind(33, 4931) !! Quadruple precision: 128-bit real integer , parameter :: kd = dp !! Internal kind used within this library. All real kinds in this library are defined with this parameter. contains end module precision © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"sourcefile/precision.f90.html","title":"precision.f90 – peano"},{"text":"assert.f90 Source File Source File assert.f90 Modules assert All Source Files airfoilBlGrowth.f90 assert.f90 constants.f90 error.f90 precision.f90 module assert !! Development tool used to ensure that predicates are as expected before proceeding. If an assertion fails, the error message will be logged and the program terminated. use error implicit none public interface assert_x_is_ge_y module procedure assert_x_is_ge_y_real , assert_x_is_ge_y_int end interface assert_x_is_ge_y interface assert_x_is_gt_y module procedure assert_x_is_gt_y_real , assert_x_is_gt_y_int end interface assert_x_is_gt_y contains subroutine assert_x_is_ge_y_real ( x , y ) !! [Real] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) real , intent ( in ) :: x , y character ( len = 256 ) :: err_s if ( x < y ) then err_s = \"Logical error - x must be greater than or equal to y.\" call add_error_message ( err_s ) call print_error_list_to_shell () call terminate_with_failure () end if end subroutine assert_x_is_ge_y_real subroutine assert_x_is_ge_y_int ( x , y ) !! [Integer] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) integer , intent ( in ) :: x , y character ( len = 256 ) :: err_s if ( x < y ) then err_s = \"Logical error - x must be greater than or equal to y.\" call add_error_message ( err_s ) call print_error_list_to_shell () call terminate_with_failure () end if end subroutine assert_x_is_ge_y_int subroutine assert_x_is_gt_y_real ( x , y ) !! Asserts that input argument /(x/) is greater than input argument /(y/) real , intent ( in ) :: x , y character ( len = 256 ) :: err_s if ( x <= y ) then err_s = \"Logical error - x must be greater than y.\" call add_error_message ( err_s ) call print_error_list_to_shell () call terminate_with_failure () end if end subroutine assert_x_is_gt_y_real subroutine assert_x_is_gt_y_int ( x , y ) !! Asserts that input argument /(x/) is greater than input argument /(y/) integer , intent ( in ) :: x , y character ( len = 256 ) :: err_s if ( x <= y ) then err_s = \"Logical error - x must be greater than y.\" call add_error_message ( err_s ) call print_error_list_to_shell () call terminate_with_failure () end if end subroutine assert_x_is_gt_y_int subroutine assert_same_bounds ( x , y ) !! Asserts that input vector /(/mathbf{x}/) has the same lower and upper bounds as input vector /(/mathbf{y}/) real , allocatable , intent ( in ) :: x (:), y (:) character ( len = 256 ) :: err_s if ( lbound ( x , 1 ) /= lbound ( y , 1 ) . or . ubound ( x , 1 ) /= ubound ( y , 1 )) then err_s = \"Array error - bounds must begin and end at the same position.\" call add_error_message ( err_s ) call print_error_list_to_shell () call terminate_with_failure () end if end subroutine assert_same_bounds subroutine assert_same_rank ( x , y ) !! Asserts that input vector /(/mathbf{x}/) has the same rank as input vector /(/mathbf{y}/) real , allocatable , intent ( in ) :: x (:), y (:) character ( len = 256 ) :: err_s if ( size ( x ) /= size ( y )) then err_s = \"Array error - ranks must be equivalent.\" call add_error_message ( err_s ) call print_error_list_to_shell () call terminate_with_failure () end if end subroutine assert_same_rank end module assert © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"sourcefile/assert.f90.html","title":"assert.f90 – peano"},{"text":"airfoilBlGrowth.f90 Source File Source File airfoilBlGrowth.f90 Modules airfoil_bl_growth All Source Files airfoilBlGrowth.f90 assert.f90 constants.f90 error.f90 precision.f90 module airfoil_bl_growth !! Calculates the boundary-layer growth over an elliptical airfoil assuming a non-viscous fluid. For these geometries, a potential-flow closed-form solution can be obtained. We assume that the boundary layer growth begins at a stagnation point. use precision use constants , only : PI implicit none real :: tau = 0.0 !! Thickness ratio integer :: nx = 0 !! Number of elements along the airfoil logical :: useMichelsCriterion , useFixedCriterion !! Criterion for handling boundary layer transition real :: xTrans = 0.0 !! Location of transition criterion (used if tc == 'F') real :: Re = 0.0 !! Reynolds number real :: xx ( 100 ) !! Domain coordinates of the boundary layer along the airfoil real :: yy ( 50 ) !! Range coordinates of the boundary layer along the airfoil real :: vgrad ( 100 ), theta ( 100 ) contains subroutine calc_bl_growth ( thicknessRatio , numElements , transitionCriterion , reynoldsNumber , transitionLocation ) !! Calculates boundary layer growth on an airfoil (modeled as an ellipse), beginning at a stagnation point. Uses Thwaite's method for the laminar flow region, Michel's method to correct transition, and Head's method for the turbulent flow region. real , intent ( in ) :: thicknessRatio !! Airfoil thickness ratio integer , intent ( in ) :: numElements !! Number of desired elements along the airfoil character , intent ( in ) :: transitionCriterion !! Boundary layer transition criterion, allowable: 'M' (Michels Criterion) or 'F' (Fixed Location) real , intent ( in ) :: reynoldsNumber !! Reynolds number (reference length), /(Re_L/) real , optional , intent ( in ) :: transitionLocation !! Boundary layer transition location along the airfoil (only used if transitionCriterion == 'F') character ( len = 256 ) :: err_s tau = thicknessRatio nx = numElements if ( transitionCriterion == 'F' ) then useFixedCriterion = . true . else if ( transitionCriterion == 'M' ) then useMichelsCriterion = . true . else err_s = \"Argument Error - Transition criterion must be either 'F' (fixed location) or 'M' (Michels Criterion).\" call add_error_message ( err_s ) call print_error_list_to_shell () call terminate_with_failure () end if re = reynoldsNumber xTrans = transitionLocation call bl_growth () end subroutine calc_bl_growth subroutine bl_growth () real :: dth2ve6 , dx , dy real :: H , L , fact , cf , lambda real :: x1 , x2 , x3 , v1 , v2 , v3 real :: Rex , Ret , RetMax integer :: i ! Calculate the domain coordinates along the airfoil, xx xx ( 1 ) = 0.0 do i = 2 , nx dx = x ( i ) - x ( i - 1 ) dy = y ( i ) - y ( i - 1 ) xx ( i ) = xx ( i - 1 ) + sqrt ( dx ** 2 + dy ** 2 ) end do ! Calculate the velocity gradient at each node v1 = ve ( 3 ); x1 = xx ( 3 ) v2 = ve ( 1 ); x2 = xx ( 1 ) ! ve(nx+1) = ve(nx-2) xx ( nx + 1 ) = xx ( nx - 2 ) do i = 1 , nx v3 = v1 ; x3 = x1 v1 = v2 ; v1 = x2 v2 = ve ( nx - 2 ) if ( i < nx ) v2 = ve ( i + 1 ) x2 = xx ( i + 1 ) fact = ( x3 - x1 ) / ( x2 - x1 ) vgrad ( i ) = ( ( v2 - v1 ) * fact - ( v3 - v1 ) / fact ) / ( x3 - x2 ) end do ! Laminar Flow Region theta ( 1 ) = sqrt ( 0.75 / Re / vgrad ( 1 )) i = 1 lambda = theta ( 1 ) ** 2 * vgrad ( 1 ) * Re do while ( lambda >= - 0.0842 ) lambda = theta ( i ) ** 2 * vgrad ( i ) * Re call thwats ( lambda , H , L ) cf = 2 * L / Re / theta ( i ) if ( i > 1 ) cf = cf / ve ( i ) print * , x ( i ), y ( i ), ve ( i ), vgrad ( i ), theta ( i ), h , cf i = i + 1 if ( i > nx ) stop dth2ve6 = 0.225 * ( ve ( i ) ** 5 + ve ( i - 1 ) ** 5 ) * ( xx ( i ) - xx ( i - 1 )) / Re theta ( i ) = sqrt ( (( theta ( i - 1 ) ** 2 ) * ( ve ( i - 1 ) ** 6 ) + dth2ve6 ) / ve ( i ) ** 6 ) if ( i == 2 ) theta ( 2 ) = theta ( 1 ) ! Test for transition if ( useFixedCriterion ) then Rex = Re * xx ( i ) * ve ( i ) Ret = Re * theta ( i ) * ve ( i ) RetMax = 1.174 * ( 1 + 22400 / Rex ) * Rex ** 0.46 end if end do end subroutine bl_growth elemental real function x ( i ) !! X-coordinates of the ellipse, !!  x_i = -\\cos\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right)  integer , intent ( in ) :: i x = - cos ( PI * ( i - 1 ) / real ( nx - 1 )) * tau end function x elemental real function y ( i ) !! Y-coordinates of the ellipse, !!  y_i = \\sin\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right)  integer , intent ( in ) :: i y = sin ( pi * ( i - 1 ) / real ( nx - 1 )) * tau end function y elemental real function ve ( i ) !! VE of the ellipse, !!  v_e = \\left(1 + \\tau\\right) \\cdot \\sqrt{\\frac{1 - x_i&#94;2}{1-\\left(1-\\tau&#94;2\\right)x_i&#94;2}}  integer , intent ( in ) :: i ve = ( 1 + tau ) * sqrt (( 1 - x ( i ) ** 2 ) / ( 1 - ( 1 - tau ** 2 ) * x ( i ) ** 2 )) end function ve end module airfoil_bl_growth © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"sourcefile/airfoilblgrowth.f90.html","title":"airfoilBlGrowth.f90 – peano"},{"text":"error.f90 Source File Source File error.f90 Modules error All Source Files airfoilBlGrowth.f90 assert.f90 constants.f90 error.f90 precision.f90 module error !! Development tool used to add and print error messages, as well as handle logic flow in the event an error is thrown. implicit none private character ( len = 1024 ) :: error_list ( 1 : 100 ) integer :: err_i = 0 contains subroutine add_error_message ( err_s ) character ( len =* ), intent ( in ) :: err_s err_i = err_i + 1 error_list ( err_i ) = \"Error: \" // trim ( err_s ) end subroutine add_error_message subroutine print_error_list_to_shell () integer :: i , l , u l = lbound ( error_list , 1 ) u = ubound ( error_list , 1 ) do i = l , u write ( 6 , * ) error_list ( i ) end do end subroutine print_error_list_to_shell subroutine terminate_with_failure () error stop 1 end subroutine terminate_with_failure end module error © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"sourcefile/error.f90.html","title":"error.f90 – peano"},{"text":"assert_x_is_ge_y_real Subroutine Source File assert.f90 assert assert_x_is_ge_y_real Variables err_s All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine assert_x_is_ge_y_real(x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description [Real] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) Variables Type Visibility Attributes Name Initial character(len=256), public :: err_s © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/assert_x_is_ge_y_real.html","title":"assert_x_is_ge_y_real – peano"},{"text":"assert_x_is_ge_y_int Subroutine Source File assert.f90 assert assert_x_is_ge_y_int Variables err_s All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine assert_x_is_ge_y_int(x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description [Integer] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) Variables Type Visibility Attributes Name Initial character(len=256), public :: err_s © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/assert_x_is_ge_y_int.html","title":"assert_x_is_ge_y_int – peano"},{"text":"assert_x_is_gt_y_real Subroutine Source File assert.f90 assert assert_x_is_gt_y_real Variables err_s All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine assert_x_is_gt_y_real(x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) Variables Type Visibility Attributes Name Initial character(len=256), public :: err_s © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/assert_x_is_gt_y_real.html","title":"assert_x_is_gt_y_real – peano"},{"text":"assert_x_is_gt_y_int Subroutine Source File assert.f90 assert assert_x_is_gt_y_int Variables err_s All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine assert_x_is_gt_y_int(x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) Variables Type Visibility Attributes Name Initial character(len=256), public :: err_s © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/assert_x_is_gt_y_int.html","title":"assert_x_is_gt_y_int – peano"},{"text":"assert_same_bounds Subroutine Source File assert.f90 assert assert_same_bounds Variables err_s All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine assert_same_bounds(x, y) Arguments Type Intent Optional Attributes Name real, intent(in), allocatable :: x (:) real, intent(in), allocatable :: y (:) Description Asserts that input vector /(/mathbf{x}/) has the same lower and upper bounds as input vector /(/mathbf{y}/) Variables Type Visibility Attributes Name Initial character(len=256), public :: err_s © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/assert_same_bounds.html","title":"assert_same_bounds – peano"},{"text":"assert_same_rank Subroutine Source File assert.f90 assert assert_same_rank Variables err_s All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine assert_same_rank(x, y) Arguments Type Intent Optional Attributes Name real, intent(in), allocatable :: x (:) real, intent(in), allocatable :: y (:) Description Asserts that input vector /(/mathbf{x}/) has the same rank as input vector /(/mathbf{y}/) Variables Type Visibility Attributes Name Initial character(len=256), public :: err_s © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/assert_same_rank.html","title":"assert_same_rank – peano"},{"text":"assert_x_is_ge_y Interface Source File assert.f90 assert assert_x_is_ge_y Module Procedures assert_x_is_ge_y_real assert_x_is_ge_y_int All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public interface assert_x_is_ge_y Module Procedures public  subroutine assert_x_is_ge_y_real (x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description [Real] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) public  subroutine assert_x_is_ge_y_int (x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description [Integer] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"interface/assert_x_is_ge_y.html","title":"assert_x_is_ge_y – peano"},{"text":"assert_x_is_gt_y Interface Source File assert.f90 assert assert_x_is_gt_y Module Procedures assert_x_is_gt_y_real assert_x_is_gt_y_int All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public interface assert_x_is_gt_y Module Procedures public  subroutine assert_x_is_gt_y_real (x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) public  subroutine assert_x_is_gt_y_int (x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"interface/assert_x_is_gt_y.html","title":"assert_x_is_gt_y – peano"},{"text":"x Function Source File airfoilBlGrowth.f90 airfoil_bl_growth x All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public elemental function x(i) Arguments Type Intent Optional Attributes Name integer, intent(in) :: i Return Value real Description X-coordinates of the ellipse,\n  x_i = -\\cos\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right)  © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/x.html","title":"x – peano"},{"text":"y Function Source File airfoilBlGrowth.f90 airfoil_bl_growth y All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public elemental function y(i) Arguments Type Intent Optional Attributes Name integer, intent(in) :: i Return Value real Description Y-coordinates of the ellipse,\n  y_i = \\sin\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right)  © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/y.html","title":"y – peano"},{"text":"ve Function Source File airfoilBlGrowth.f90 airfoil_bl_growth ve All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public elemental function ve(i) Arguments Type Intent Optional Attributes Name integer, intent(in) :: i Return Value real Description VE of the ellipse,\n  v_e = \\left(1 + \\tau\\right) \\cdot \\sqrt{\\frac{1 - x_i&#94;2}{1-\\left(1-\\tau&#94;2\\right)x_i&#94;2}}  © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/ve.html","title":"ve – peano"},{"text":"calc_bl_growth Subroutine Source File airfoilBlGrowth.f90 airfoil_bl_growth calc_bl_growth Variables err_s All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine calc_bl_growth(thicknessRatio, numElements, transitionCriterion, reynoldsNumber, transitionLocation) Arguments Type Intent Optional Attributes Name real, intent(in) :: thicknessRatio Airfoil thickness ratio integer, intent(in) :: numElements Number of desired elements along the airfoil character, intent(in) :: transitionCriterion Boundary layer transition criterion, allowable: 'M' (Michels Criterion) or 'F' (Fixed Location) real, intent(in) :: reynoldsNumber Reynolds number (reference length), /(Re_L/) real, intent(in), optional :: transitionLocation Boundary layer transition location along the airfoil (only used if transitionCriterion == 'F') Description Calculates boundary layer growth on an airfoil (modeled as an ellipse), beginning at a stagnation point. Uses Thwaite's method for the laminar flow region, Michel's method to correct transition, and Head's method for the turbulent flow region. Variables Type Visibility Attributes Name Initial character(len=256), public :: err_s © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/calc_bl_growth.html","title":"calc_bl_growth – peano"},{"text":"bl_growth Subroutine Source File airfoilBlGrowth.f90 airfoil_bl_growth bl_growth Variables dth2ve6 dx dy H L fact cf lambda x1 x2 x3 v1 v2 v3 Rex Ret RetMax i All Procedures assert_same_bounds assert_same_rank assert_x_is_ge_y assert_x_is_ge_y_int assert_x_is_ge_y_real assert_x_is_gt_y assert_x_is_gt_y_int assert_x_is_gt_y_real bl_growth calc_bl_growth ve x y public  subroutine bl_growth() Arguments None Variables Type Visibility Attributes Name Initial real, public :: dth2ve6 real, public :: dx real, public :: dy real, public :: H real, public :: L real, public :: fact real, public :: cf real, public :: lambda real, public :: x1 real, public :: x2 real, public :: x3 real, public :: v1 real, public :: v2 real, public :: v3 real, public :: Rex real, public :: Ret real, public :: RetMax integer, public :: i © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"proc/bl_growth.html","title":"bl_growth – peano"},{"text":"constants Module Source File constants.f90 constants Variables PI All Modules airfoil_bl_growth assert constants error precision Uses: precision Provides constants used within this library. Variables Type Visibility Attributes Name Initial real, public, parameter :: PI = 3.1415926535898 © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"module/constants.html","title":"constants – peano"},{"text":"precision Module Source File precision.f90 precision Variables sp dp qp kd All Modules airfoil_bl_growth assert constants error precision Uses: ISO_FORTRAN_ENV Provides kind attributes for setting machine-compiler-independent precision for real numbers. Here, we define \"single precision\" to mean 32-bit precision, \"double precision\" to mean 64-bit precision, and \"quadruple precision\" to mean 128-bit precision. This ensures that precision is preserved regardless of the compiler or computer architecture. Variables Type Visibility Attributes Name Initial integer, public, parameter :: sp = REAL32 Single precision: 32-bit real integer, public, parameter :: dp = REAL64 Double precision: 64-bit real integer, public, parameter :: qp = REAL128 Quadruple precision: 128-bit real integer, public, parameter :: kd = dp Internal kind used within this library. All real kinds in this library are defined with this parameter. © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"module/precision.html","title":"precision – peano"},{"text":"assert Module Source File assert.f90 assert Interfaces assert_x_is_ge_y assert_x_is_gt_y Subroutines assert_x_is_ge_y_real assert_x_is_ge_y_int assert_x_is_gt_y_real assert_x_is_gt_y_int assert_same_bounds assert_same_rank All Modules airfoil_bl_growth assert constants error precision Uses: error Development tool used to ensure that predicates are as expected before proceeding. If an assertion fails, the error message will be logged and the program terminated. Interfaces public interface assert_x_is_ge_y public  subroutine assert_x_is_ge_y_real (x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description [Real] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) public  subroutine assert_x_is_ge_y_int (x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description [Integer] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) public interface assert_x_is_gt_y public  subroutine assert_x_is_gt_y_real (x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) public  subroutine assert_x_is_gt_y_int (x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) Subroutines public  subroutine assert_x_is_ge_y_real (x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description [Real] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) public  subroutine assert_x_is_ge_y_int (x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description [Integer] Asserts that input argument /(x/) is greater than or equal to input argument /(y/) public  subroutine assert_x_is_gt_y_real (x, y) Arguments Type Intent Optional Attributes Name real, intent(in) :: x real, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) public  subroutine assert_x_is_gt_y_int (x, y) Arguments Type Intent Optional Attributes Name integer, intent(in) :: x integer, intent(in) :: y Description Asserts that input argument /(x/) is greater than input argument /(y/) public  subroutine assert_same_bounds (x, y) Arguments Type Intent Optional Attributes Name real, intent(in), allocatable :: x (:) real, intent(in), allocatable :: y (:) Description Asserts that input vector /(/mathbf{x}/) has the same lower and upper bounds as input vector /(/mathbf{y}/) public  subroutine assert_same_rank (x, y) Arguments Type Intent Optional Attributes Name real, intent(in), allocatable :: x (:) real, intent(in), allocatable :: y (:) Description Asserts that input vector /(/mathbf{x}/) has the same rank as input vector /(/mathbf{y}/) © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"module/assert.html","title":"assert – peano"},{"text":"airfoil_bl_growth Module Source File airfoilBlGrowth.f90 airfoil_bl_growth Variables tau nx useMichelsCriterion useFixedCriterion xTrans Re xx yy vgrad theta Functions x y ve Subroutines calc_bl_growth bl_growth All Modules airfoil_bl_growth assert constants error precision Uses: precision constants Calculates the boundary-layer growth over an elliptical airfoil assuming a non-viscous fluid. For these geometries, a potential-flow closed-form solution can be obtained. We assume that the boundary layer growth begins at a stagnation point. Variables Type Visibility Attributes Name Initial real, public :: tau = 0.0 Thickness ratio integer, public :: nx = 0 Number of elements along the airfoil logical, public :: useMichelsCriterion logical, public :: useFixedCriterion Criterion for handling boundary layer transition real, public :: xTrans = 0.0 Location of transition criterion (used if tc == 'F') real, public :: Re = 0.0 Reynolds number real, public :: xx (100) Domain coordinates of the boundary layer along the airfoil real, public :: yy (50) Range coordinates of the boundary layer along the airfoil real, public :: vgrad (100) real, public :: theta (100) Functions public elemental function x (i) Arguments Type Intent Optional Attributes Name integer, intent(in) :: i Return Value real Description X-coordinates of the ellipse,\n  x_i = -\\cos\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right)  public elemental function y (i) Arguments Type Intent Optional Attributes Name integer, intent(in) :: i Return Value real Description Y-coordinates of the ellipse,\n  y_i = \\sin\\left(\\frac{(i-1)\\pi}{nx-1}\\tau\\right)  public elemental function ve (i) Arguments Type Intent Optional Attributes Name integer, intent(in) :: i Return Value real Description VE of the ellipse,\n  v_e = \\left(1 + \\tau\\right) \\cdot \\sqrt{\\frac{1 - x_i&#94;2}{1-\\left(1-\\tau&#94;2\\right)x_i&#94;2}}  Subroutines public  subroutine calc_bl_growth (thicknessRatio, numElements, transitionCriterion, reynoldsNumber, transitionLocation) Arguments Type Intent Optional Attributes Name real, intent(in) :: thicknessRatio Airfoil thickness ratio integer, intent(in) :: numElements Number of desired elements along the airfoil character, intent(in) :: transitionCriterion Boundary layer transition criterion, allowable: 'M' (Michels Criterion) or 'F' (Fixed Location) real, intent(in) :: reynoldsNumber Reynolds number (reference length), /(Re_L/) real, intent(in), optional :: transitionLocation Boundary layer transition location along the airfoil (only used if transitionCriterion == 'F') Description Calculates boundary layer growth on an airfoil (modeled as an ellipse), beginning at a stagnation point. Uses Thwaite's method for the laminar flow region, Michel's method to correct transition, and Head's method for the turbulent flow region. public  subroutine bl_growth () Arguments None © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"module/airfoil_bl_growth.html","title":"airfoil_bl_growth – peano"},{"text":"error Module Source File error.f90 error All Modules airfoil_bl_growth assert constants error precision Development tool used to add and print error messages, as well as handle logic flow in the event an error is thrown. © 2015 peano was written by Vincent San Miguel. Documentation generated by FORD .","tags":"","loc":"module/error.html","title":"error – peano"}]}