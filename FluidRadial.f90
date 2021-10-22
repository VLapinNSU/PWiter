! make matrix A(w) where A(w)p = (w^n+1-w^n)/dt for radial flow in flat channel, and calculates pressure form width
! elements are constant, cells are filled by constant values, A - 3diagonal.
module fluidRadial
implicit none 

    real(8), parameter, public :: Pi = 3.1415926535897932384626433832795028841971d0  ! Pi number

    type TfluidParams
        real(8) :: qin, qout            ! inflow and outflow rate (should be set 2 of 4)
        real(8) :: pin, pout            ! inflow and outflow pressure (should be set 2 of 4)
        real(8) :: mu
        real(8) :: dt
        integer :: BoundaryTypes        ! which of p and q are set: 11 - qin, qout (incorrect problem without elastic)
                                        ! 12 - qin, poit, 21 - pin, qout, 22 - pin, pout
    end type

    contains     

!=================================================================================!
! Subprogram: test matrix calculating for radial flow in flat channel
!              A(w)p = (w^n+1-w^n)/dt
! The subroutine check matrix calculations with qin and pout conditions 
!=================================================================================! 
subroutine testFluidMatrixRadial()    
real(8), allocatable :: xx(:)                           ! mesh, cell boundaries [0;Rfrac], [0..NN]
real(8), allocatable :: press(:), pressEx(:)            ! pressure, width at cell centres [1..NN]
real(8), allocatable :: width(:), widthN(:)             ! width at current and previous steps [1..NN]
real(8), allocatable :: matrAp(:,:),matrApStorage(:,:)  ! fluid matrix [3,1..NN] with boundary conditions
real(8), allocatable :: RHS(:,:),RHSstorage(:,:)        ! Right hand side [1,1..NN] with boundary conditions
integer :: NN                                           ! mesh size
type(TfluidParams) :: fluidParams                       ! paranmeters
integer :: i
real(8) :: Err, MaxSol                                  ! Error 
real(8) :: rr, rout                                     ! tmp radius 
integer :: info                                         ! tmp for matrix SLAE solution

    NN = 2
    allocate(xx(0:NN))
    allocate(width(1:NN),widthN(1:NN),press(1:NN),pressEx(1:NN))
    allocate(matrAp(3,1:NN),matrApStorage(3,1:NN),RHS(1,1:NN),RHSstorage(1,1:NN))
    call setFluidParamsForTest(fluidParams)
    xx = (/ (i, i = 0, NN) /) * 10.d0/NN !+ 2.d0
    width = 0.01d0    ;     widthN = width
    
    call makeFluidMatrixAndRHSRadial(xx,NN,fluidParams,width,widthN, matrAp,RHS)
    matrApStorage = matrAp      ;   RHSstorage = RHS
    call dgtsv(NN,1,matrAp(1,2:NN),matrAp(2,1:NN),matrAp(3,1:NN-1),RHS(1,1:NN),NN,info )
    !check solution of dgtsv
    Err = 0.d0
    do i = 1, NN
        rr = matrApStorage(2,i)*RHS(1,i) - RHSstorage(1,i)
        if (i> 1) rr = rr + matrApStorage(1,i)*RHS(1,i-1)
        if (i<NN) rr = rr + matrApStorage(3,i)*RHS(1,i+1)
        Err = max(Err,rr)
    enddo
    Err = Err / maxval(abs(RHSstorage(1,1:NN)))     ! relative residual 
    press(1:NN) = RHS(1,1:NN)

    ! exact pressure for constant width and qin, pout conditions
    rout = 0.5d0*(xx(NN-1)+xx(NN))
    open(1, file = 'test.plt')
    do i = 1, NN
        rr = 0.5d0*(xx(i-1)+xx(i))
        pressEx(i) = fluidParams%pout + 6.d0*fluidParams%mu*fluidParams%qin / (pi*width(1)**3) * (log(rout)-log(rr))
        write(1,*) xx(i), press(i), pressEx(i)
    enddo
    Err = maxval(abs(pressEx-press)) 
    Err = Err / maxval(abs(pressEx))
    close(1)
    
end subroutine testFluidMatrixRadial    

! subroutine fills matrix and RHS and boundary conditions 
! A(w)p = (w^n+1-w^n)/dt
! equation for aproximation 
!  \frac{1}{r} \frac{\partial}{\partial r} \left( \frac{r width^3}{12 \mu} \frac{\partial p}{\partial r} \right) = 
!  \frac{\partial width}{\partial t} 
! boundary conditions for approximations
! p(N) = pout   ;   qin 
! TODO for Lapin: implement others types of boundary conditions
! TODO for Lapin: only 1-st order of convergence is observed - fix approximation
subroutine makeFluidMatrixAndRHSRadial(xx,NN,fluidParams,width,widthN, matrAp,RHS)
real(8), intent(IN ) :: xx(0:)                           ! mesh, cell boundaries [0;Rfrac], [0..NN]
real(8), intent(IN ) :: width(:), widthN(:)             ! width at current and previous steps [1..NN]
integer, intent(IN ) :: NN                              ! mesh size
type(TfluidParams), intent(IN ) :: fluidParams          ! paranmeters
real(8), intent(OUT) :: matrAp(:,:)                     ! fluid matrix [3,1..NN] with boundary conditions
real(8), intent(OUT) :: RHS(:,:)                        ! Right hand side [1,1..NN] with boundary conditions
integer :: i
real(8) :: rp, rm, wp, wm, wm3r, wp3r, hh                           ! r_{i+-1}, w_{i+-1}, hh=dx

    ! equation approximation
    hh = xx(1)-xx(0)
    matrAp = 0.d0   ;   rhs = 0.d0
    do i = 1, NN    
        !if (i> 1) then  ;   wm = 0.5d0*(width(i-1)+width(i))    ;   rm = 0.5d0*(xx(i-1)+xx(i))  ;   endif
        !if (i<NN) then  ;   wp = 0.5d0*(width(i+1)+width(i))    ;   rp = 0.5d0*(xx(i+1)+xx(i))  ;   endif
        !if (i> 1) matrAp(1,i) = 1.d0/xx(i) *  ( rm*wm**3 / (12.d0*fluidParams%mu) ) / hh**2
        !if (i<NN) matrAp(3,i) = 1.d0/xx(i) *  ( rp*wp**3 / (12.d0*fluidParams%mu) ) / hh**2
        if (i> 1)   wm3r = 0.5d0*(width(i-1)**3*xx(i-1) + width(i)**3*xx(i))
        if (i<NN)   wp3r = 0.5d0*(width(i+1)**3*xx(i+1) + width(i)**3*xx(i))
        if (i> 1) matrAp(1,i) = 1.d0/xx(i) *  ( wm3r / (12.d0*fluidParams%mu) ) / hh**2
        if (i<NN) matrAp(3,i) = 1.d0/xx(i) *  ( wp3r / (12.d0*fluidParams%mu) ) / hh**2
        matrAp(2,i) = - (matrAp(3,i)+matrAp(1,i))
        RHS(1,i) = (width(i)-widthN(i))/fluidParams%dt
    enddo
    ! inflow boundary
    RHS(1,1) = RHS(1,1) - fluidParams%qin / (2.d0*pi) / xx(1) / hh      ! аппроксимация потока в i+1/2 уже есть
    ! outflow boundary
    RHS(1,NN) = fluidParams%pout
    matrAp(2,NN) = 1.d0     ;       matrAp(1,NN) = 0.d0

end subroutine makeFluidMatrixAndRHSRadial





subroutine setFluidParamsForTest(fluidParams)    
type(TfluidParams) :: fluidParams                       ! paranmeters
    fluidParams%qin     = 0.01d0
    fluidParams%pout    = 0.d0
    fluidParams%mu      = 0.1d0
    fluidParams%dt      = 1.d0
    fluidParams%BoundaryTypes = 12      ! 11 - qin, qout, 12 - qin, poit, 21 - pin, qout, 22 - pin, pout
end subroutine setFluidParamsForTest








end module fluidRadial