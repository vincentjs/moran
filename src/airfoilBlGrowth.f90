module airfoil_bl_growth
  !! Calculates the boundary-layer growth over an elliptical airfoil assuming a non-viscous fluid. For these geometries, a potential-flow closed-form solution can be obtained. We assume that the boundary layer growth begins at a stagnation point. 

  use precision
  use constants, only : PI
  
  implicit none

  real :: tau = 0.0
  !! Thickness ratio
  integer :: nx = 0
  !! Number of elements along the airfoil
  logical :: useMichelsCriterion, useFixedCriterion
  !! Criterion for handling boundary layer transition
  real :: xTrans = 0.0
  !! Location of transition criterion (used if tc == 'F')
  real :: Re = 0.0
  !! Reynolds number
  
  real :: xx(100)
  !! Domain coordinates of the boundary layer along the airfoil
  real :: yy(50)
  !! Range coordinates of the boundary layer along the airfoil
  real :: vgrad(100), theta(100) 
  
contains

  subroutine calc_bl_growth(thicknessRatio, numElements, transitionCriterion, reynoldsNumber, transitionLocation)
    !! Calculates boundary layer growth on an airfoil (modeled as an ellipse), beginning at a stagnation point. Uses Thwaite's method for the laminar flow region, Michel's method to correct transition, and Head's method for the turbulent flow region.
    
    real, intent(in) :: thicknessRatio
    !! Airfoil thickness ratio
    integer, intent(in) :: numElements
    !! Number of desired elements along the airfoil
    character, intent(in) :: transitionCriterion
    !! Boundary layer transition criterion, allowable: 'M' (Michels Criterion) or 'F' (Fixed Location)
    real, intent(in) :: reynoldsNumber
    !! Reynolds number (reference length), /(Re_L/)
    real, optional, intent(in) :: transitionLocation
    !! Boundary layer transition location along the airfoil (only used if transitionCriterion == 'F')

    character(len=256) :: err_s
    
    tau = thicknessRatio
    nx = numElements
    if (transitionCriterion == 'F') then
       useFixedCriterion = .true.
    else if (transitionCriterion == 'M') then
       useMichelsCriterion = .true.
    else
       err_s = "Argument Error - Transition criterion must be either 'F' (fixed location) or 'M' (Michels Criterion)."
       call add_error_message(err_s)
       call print_error_list_to_shell()
       call terminate_with_failure()
    end if

    re = reynoldsNumber
    xTrans = transitionLocation
    
    call bl_growth()
    
  end subroutine calc_bl_growth

  subroutine bl_growth()

    real :: dth2ve6, dx, dy
    real :: H, L, fact, cf, lambda
    real :: x1, x2, x3, v1, v2, v3
    real :: Rex, Ret, RetMax
    real :: rtheta, H10FH, H0FH1, H

    integer :: i, itrans
    
    ! Calculate the domain coordinates along the airfoil, xx
    xx(1) = 0.0
    do i = 2, nx
       dx = x(i) - x(i-1)
       dy = y(i) - y(i-1)
       xx(i) = xx(i-1) + sqrt(dx**2 + dy**2)
    end do

    ! Calculate the velocity gradient at each node
    v1 = ve(3); x1 = xx(3)
    v2 = ve(1); x2 = xx(1)
    ! ve(nx+1) = ve(nx-2)
    xx(nx+1) = xx(nx-2)

    do i = 1, nx
       v3 = v1; x3 = x1
       v1 = v2; v1 = x2

       v2 = ve(nx-2)
       if (i < nx) v2 = ve(i+1)

       x2 = xx(i+1)
       fact = (x3 - x1) / (x2 - x1)

       vgrad(i) = ( (v2 - v1)*fact - (v3-v1)/fact ) / (x3 - x2)
    end do

    ! Boundary layer growth, laminar and turbulent
    theta(1) = sqrt(0.75 / Re / vgrad(1))

    if (useMichelsCriterion) ret = 0.0; retMax = 1E10
    
    i = 1

    ! Laminar Flow Region
    do while (lambda >= -0.0842 .and. &
       (useFixedCriterion .and. x(i) < xTrans) .or. &
       (useMichelsCriterion .and. ret < retMax))

       lambda = theta(i)**2 * vgrad(i) * Re

       if (lambda >= -0.0842) then
          call thwats(lambda, H, L)

          cf = 2 * L / Re / theta(i)
          if (i > 1) cf = cf / ve(i)

          print *, x(i), y(i), ve(i), vgrad(i), theta(i), h, cf

          i = i + 1

          if (i > nx) stop

          dth2ve6 = 0.225 * (ve(i)**5 + ve(i-1)**5) * (xx(i) - xx(i-1)) / Re
          theta(i) = sqrt( ((theta(i-1)**2)*(ve(i-1)**6) + dth2ve6) / ve(i)**6 )

          if (i == 2) theta(2) = theta(1)

          ! Test for transition
          if (useFixedCriterion) then
             Rex = Re*xx(i)*ve(i)
             Ret = Re*theta(i)*ve(i)
             RetMax = 1.174 * (1 + 22400/Rex) * Rex**0.46
          end if
       end if
    end do

    if (lambda < -0.0842) then
       print *, "Laminar Separation"
    else
       ! Turbulent Flow Region
       H = 1E10
       
       do while (H >= 1.0)
          itrans = i
          print *, "Input H at transition: "
          read H

          if (H < 1.0) stop

          i = itrans
          yy(2) = H10FH(H)
          yy(1) = theta(i-1)

          do while (H <= 2.4 .and. i <= nx)
             dx = xx(i) - xx(i-1)

             call runge2(i-1, i, dx, yy, 2)

             theta(i) = yy(1)
             H = H0FH1(yy(2))
             rtheta = Re * ve(i) * theta(i)
             cf = cfturb(rtheta, H)

             print *, x(i), y(i), ve(i), vgrad(i), theta(i), H, CF

             i = i + 1
          end do

          if (H > 2.4) then
             print *, "Turbulent separation."
          end if
       end do
    end if
    
  end subroutine bl_growth

  pure subroutine thwats(lambda, H, L)
    !! Thwaite's Correlation Formulas

    real, intent(in) :: lambda
    real, intent(out) :: H, L

    if (lambda < 0.0) then
       L = 0.22 + 1.402*lambda + 0.18*lambda / (0.107 + lambda)
       H = 2.088 + 0.0731 / (0.14 + lambda)
    else
       L = 0.22 + lambda*(1.57 - 1.8*lambda)
       H = 2.61 - lambda*(3.75 - 5.24*lambda)
    end if
    
  end subroutine thwats

  elemental real function H10FH(H)
    !! Head's Correlation Formula for H1(H)

    if (H > 1.6) then
       H10FH = 3.3 + 1.5501 * (H - 0.6778)**-3.064
    else
       H10FH = 3.3 + 0.8234 * (H - 1.1)**-1.287
    end if

  end function H10FH

  elemental real function H0FH1(H1)
    !! Inverse of H1(H)

    if (H1 < 3.3) then
       H0FH1 = 3.0
    else if (H1 >= 3.3 .or. H1 < 5.3) then
       H0FH1 = 0.6778 + 1.1536 * (H1 - 3.3)**-0.326
    else
       H0FH1 = 1.1 + 0.86 * (H1 - 3.3)**-0.777
    end if

  end function H0FH1

  pure real function cfturb(rtheta, H)
    !! Ludwieg-Tillman Skin Friction Formula

    cfturb = 0.246 * (10.0**(-0.678*H)) * (rtheta**-0.268)

  end function cfturb

  subroutine derivs(i)
    !! Set derivatives of vector Y
    H1 = yt(2)
    if (H1 > 3.0) then
       H = H0FH1(H1)
       rtheta = re * ve(i) * yt(1)
       yp(1) = -(H + 2.0) * yt(1) * vgrad(i) / ve(i) + 0.5 * cfturb(rtheta,H)
       yp(2) = -H1 * (vgrad(i) / ve(i) + yp(1)/yt(1)) + 0.0306 * ((H1-3)**-0.6169) / yt(1)
    end if

  end subroutine derivs
      
  elemental real function x(i)
    !! X-coordinates of the ellipse,
    !! $$ x_i = -\cos\left(\frac{(i-1)\pi}{nx-1}\tau\right) $$

    integer, intent(in) :: i

    x = -cos(PI * (i-1) / real(nx - 1)) * tau
  end function x

  elemental real function y(i)
    !! Y-coordinates of the ellipse,
    !! $$ y_i = \sin\left(\frac{(i-1)\pi}{nx-1}\tau\right) $$
    
    integer, intent(in) :: i

    y = sin(pi*(i-1) / real(nx-1)) * tau
  end function y
  
  elemental real function ve(i)
    !! VE of the ellipse,
    !! $$ v_e = \left(1 + \tau\right) \cdot \sqrt{\frac{1 - x_i^2}{1-\left(1-\tau^2\right)x_i^2}} $$
    integer, intent(in) :: i
    
    ve = (1 + tau) * sqrt((1 - x(i)**2) / (1 - (1-tau**2)*x(i)**2))
  end function ve
  
  
end module airfoil_bl_growth
