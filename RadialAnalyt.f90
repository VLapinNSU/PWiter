module radialAnalyt
implicit none 

    integer, parameter, public :: DP = SELECTED_REAL_KIND(15)                     ! Double precision real kind

    contains

! gives time, when radial fracture with given parameters reaches the given radius 
    function countTimeOfRadialR(qin,Ep,mu, R)
    real(DP), intent(IN) :: qin     ! inflow rate
    real(DP), intent(IN) :: Ep      ! elastic E/(1-nu^2)
    real(DP), intent(IN) :: mu      ! viscosity
    real(DP), intent(IN) :: R       ! radius
    real(DP) :: countTimeOfRadialR
        countTimeOfRadialR = ( (R/0.6944d0)**9 / (qin**3 * Ep) * (mu*12) )**(0.25d0)
    end function

! calcululates radius of the fracture using analytical formula for viscosity regime. Is taken from DoncovEV_2016_approximate
    function countROfRadialTime(qin, Ep, mu, t)
    real(DP), intent(IN) :: qin     ! inflow rate
    real(DP), intent(IN) :: Ep      ! elastic E/(1-nu^2)
    real(DP), intent(IN) :: mu      ! viscosity
    real(DP), intent(IN) :: t       ! time
    real(DP) :: countROfRadialTime
    countROfRadialTime =  0.6944d0 * ( qin**3 * Ep * t**4 / (mu*12) )**(1.d0/9.d0)
    end function
    
    function countWidthFromRadialAnalyt(mu, q_in, Ep, t, rho)
    ! calcululates width of the fracture using analytical formula for viscosity regime. Is taken from DoncovEV_2016_approximate
    real(DP), intent(IN) :: mu, q_in, Ep
    real(DP) :: t       ! time 
    real(DP) :: rho     ! ralative point radius R/Rfrac
    real(DP) :: countWidthFromRadialAnalyt
    if (rho>1.d0-0.01d0) then 
        countWidthFromRadialAnalyt = 0.d0
    else
        
        countWidthFromRadialAnalyt =  1.1901d0 * ( (mu*12.d0)**2 * q_in**3 * t / Ep**2 )**(1.d0/9.d0) &
                          * (1.d0 + rho)**(0.487d0) * (1.d0 - rho)**(2.d0/3.d0)
    endif
    end function countWidthFromRadialAnalyt


end module radialAnalyt
