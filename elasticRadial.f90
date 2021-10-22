! make matrix T where p = T*w for radial fracture, and calculates width from pressure
! elements are constant
 module elasticRadial
implicit none 
    
    contains     

!=================================================================================!
! Subprogram: test matrix calculating for elastic radial fracture
!             W_k = T_k_i * p_i 
! The subroutine calculates matrix matrWfromP where T_k_i = matrWfromP * 4*(1-\nu)/(\pi*\mu)
!=================================================================================! 
subroutine testElasticMatrixRadial()    
real(8), allocatable :: xx(:)                           ! mesh, cell boundaries [0;Rfrac], [0..NN]
real(8), allocatable :: press(:), width(:), widthEx(:)  ! pressure, width at cell centres [1..NN]
real(8), allocatable :: matrWfromP(:,:)                 ! elastic matrix [1..NN,1..NN]
integer :: NN                                           ! mesh size
real(8) :: Rfrac                                        ! fracture radius
integer :: i
real(8) :: ElasticCoef                                  ! 8/pi/E1, E1 = E/(1-\nu^2)
real(8) :: Err                                          ! Error
real(8) :: rr                                           ! tmp radius 

    NN = 2
    allocate(xx(0:NN))
    allocate(width(1:NN),press(1:NN),widthEx(1:NN))
    allocate(matrWfromP(1:NN,1:NN))
    Rfrac = 10.d0
    press = 1.d6
    xx = (/ (i, i = 0, NN) /) * Rfrac/NN
    ElasticCoef = 8.d0/3.14159265d0/20.d9*(1-0.25d0**2)
    
    call makeElasticMatrixRadial(xx,Rfrac,NN, matrWfromP)
    width(1:NN) = matmul(matrWfromP(1:NN,1:NN),press(1:NN))*ElasticCoef
    do i = 1, NN
        rr = 0.5d0*(xx(i-1)+xx(i))
        widthEx(i) = ElasticCoef * Rfrac * press(1) * sqrt(1-rr**2/Rfrac**2)
    enddo
    Err = maxval(abs(widthEx-width)) 
    Err = Err / maxval(abs(widthEx))
    
end subroutine testElasticMatrixRadial    
    
!=================================================================================!
! Subprogram: counts elements of auxiliary matrices of width-pressure connection 
!             W_k = T_k_i * p_i 
! The subroutine calculates matrix matrWfromP where T_k_i = matrWfromP * 4*(1-\nu)/(\pi*\mu)
!=================================================================================! 
subroutine makeElasticMatrixRadial(xx,Rfrac,NN, matrWfromP)
real(8), intent(IN ) :: xx(0:NN)            ! mesh [R_w;R_out], [0..NN]
integer, intent(IN ) :: NN                  ! mesh size
real(8), intent(IN ) :: Rfrac               ! frac radius

real(8), intent(OUT) :: matrWfromP(:,:)     ! matrix for w_k = T^s_k*p_s

real(8), allocatable :: R_c(:,:)                    ! Auxiliary matrices

real(8) :: ww(1:6), xi(1:6)  ! Weights and nodes in Gaussian formula
real(8) :: rr       ! 
integer :: k,i,j

    allocate(R_c(1:NN,1:NN))
    ! calculate nodes and weights for Gauss quadrature rule
    call Gauss_points(6, xi, ww )
    ! Calculate auxilary matrices R_a, R_b, R_c to combine them into connection matrix T_pw further
    do k = 1, NN
        ! for each point (k)
        do i = 1, NN
        do j = 1, NN
            rr = 0.5d0*(xx(k)+xx(k-1))
            R_c(i,j) = make_R(xx,i,j,1,rr,xi,ww,6)
        enddo
        enddo 
        ! Calculate T_pw(k,l,j,s)
        do j = 1, NN
            matrWfromP(k,j) = 0.d0
            do i = 1, NN
                matrWfromP(k,j) = matrWfromP(k,j) + R_c(i,j)
            enddo
        enddo 
        continue  ! matrWfromP(k,*) is calculated
    enddo
    deallocate(R_c)

end subroutine makeElasticMatrixRadial



!=================================================================================!
! Subprogram: counts elements of auxiliary matrices of width-pressure connection 
!             W_k = T_k_i * p_i 
! The subroutine is required to calculate integrals for obtaining elements
! of matrix T_k_i
!=================================================================================! 
    
function make_R(xx,i,j,ind,rr,xi,ww,Ngauss) result(Rez)
! External integral (see Slepyan 'Fracture Mechanics' in Russian, section: Radial fracture loaded by pressure)
! Integrate using Gaussian formula. We consider that nodes and weights (xi,ww) have been already calculated (call Gauss_points in advance)
! integral from x=elem(i)xb to x=elem(i)xe ( J_j^ind(x) / sqrt(x**2-rr**2)) dx
! integral from x=xx(i-1)   to x=xx(i)     ( J_j^ind(x) / sqrt(x**2-rr**2)) dx
! see J_j in procedure make_J
implicit none
! input
real(8), intent(IN) :: xx(0:)                        ! mesh [R_w;R_out], [0..NN]
integer, intent(IN) :: i,j                          ! Number of elements and number of inside integrals
integer, intent(IN) :: ind                          ! Power of function under integral
real(8), intent(IN) :: rr                           ! Parameter of function under integral
integer, intent(IN) :: Ngauss                       ! Quantity of nodes in Gaussian quadrature
real(8), intent(IN) :: ww(1:Ngauss), xi(1:Ngauss)   ! Weights and nodes in Gaussian formulas
!
! output
real(8) :: Rez                             ! Value of integral
! tmp
integer :: k                                        ! Index
double precision :: a,b                             ! Limits of integration 
double precision :: sum, sum_1, int_rr 
double precision :: xl, int_jj, ff                  ! Temporary variables

    if ( rr.lt.xx(i-1)-1.d-10 ) then 
        ! calculate using regular Gaussian formula
        a = xx(i-1)   ;   b = xx(i)   ;   sum = 0.d0
        do k = 1, Ngauss 
            xl = a + (b-a) * (xi(k)+1) / 2.d0
            int_Jj = make_J(xx(j-1),xx(j),xl,ind)
            ff = int_Jj / dsqrt(xl**2-rr**2)
            sum = sum + ff*ww(k)
        enddo 
        Rez = sum * 0.5d0*(b-a)
    elseif ( rr.le.xx(i)-1.d-10 ) then 
        ! it is necessary to extract singularity - 0 in denominator on boundary
        ! calculate part containing singularity
        int_Jj = make_J(xx(j-1),xx(j),rr,ind)
        if (dabs(int_Jj).ge.1.d-10) then 
            int_rr = int_Jj / dsqrt(2.d0*rr)
            sum_1 = int_rr * 2.d0 * dsqrt(xx(i)-rr)
        else
            int_rr = 0.d0
            sum_1 = 0.d0
        endif
        ! then calculate integral for funciton without singularity 
        a = rr    ;   b = xx(i)   ;   sum = 0.d0
        do k = 1, Ngauss 
            xl = a + (b-a) * (xi(k)+1) / 2.d0
            int_Jj = make_J(xx(j-1),xx(j),xl,ind)
            ff = ( int_Jj / dsqrt(xl+rr) - int_rr ) / dsqrt(xl-rr)
            sum = sum + ff*ww(k)
        enddo 
        Rez = sum * 0.5d0*(b-a) + sum_1
    elseif ( rr.gt.xx(i)-1.d-10 ) then 
        ! integral is equal to zero
        Rez = 0.d0
    else
        ! unexpected error
        write(*,*) 'make_R : unexpected choice ', i
        stop
   endif
 
end function make_R
    
    
!=================================================================================!
! Subprogram: analytically calculates internal 
!=================================================================================! 
function make_J(a,b_int,xl,ind) result(int_Jj)
! integral from ro = a to ro = b_int ( ro**ind / sqrt(xl**2-ro**2) d ro
! it is defined that if (xl-ro<0) then (ro**ind / sqrt(xl**2-ro**2)) = 0 
! it is the feature of solution (see Slepyan 'Fracture Mechanics' in Russian)
! input
real(8), intent(IN) :: a, b_int                         ! limits of integral 
real(8), intent(IN) :: xl                               ! value of parameter in denominator, see expression for integral above
integer, intent(IN) :: ind                              ! value with a power rî, see expression for integral above
! output 
real(8) :: int_Jj                                       ! Value of integral
! tmp 
real(8) :: b                               ! Left limit of integration, where we actually integrate
    if (xl.le.a) then
        ! value in radical expression is negative so integral is equal to zero
        int_Jj = 0.d0
        return
    elseif (xl.le.b_int) then 
        ! parameter is within the integration interval then
        ! integrate only to the left of the parameter
        b = xl
    elseif (xl.ge.b_int) then 
        ! parameter is out of the integration interval
        ! function is without singularity
        b = b_int
    else 
        ! unexpected error
        write(*,*) 'make_J : unexpected choice ', a, ' < ', xl, ' < ', b_int
        stop
    endif 
   if (ind.eq.1) then 
      ! integral of the first degree (for constant functions) 
      int_Jj = dsqrt(xl**2-a**2) - dsqrt(xl**2-b**2)
   elseif (ind.eq.2) then 
      ! integral of the second degree (for linear functions) 
      int_Jj = 0.5d0 * ( a*dsqrt(xl**2-a**2) - b*dsqrt(xl**2-b**2) ) &
               + 0.5d0*xl**2 * ( dacos(a/xl) - dacos(b/xl) )
   elseif (ind.eq.3) then 
      ! integral of the third degree (for parabolas)
      int_Jj = ( a**2*dsqrt(xl**2-a**2) - b**2*dsqrt(xl**2-b**2) ) / 3.d0 &
               + 2.d0 * xl**2 * ( dsqrt(xl**2-a**2) - dsqrt(xl**2-b**2) ) / 3.d0
   else 
      ! unexpected error
      int_Jj = 0.d0
      write(*,*) 'make_J : unexpected index = ', ind
      stop
   endif
end function make_J


!****************************************************************************
! Soubroutine to compute the Gausse quadrature rule
!****************************************************************************
subroutine Gauss_points( order, x, w )
! The weight function is w(x) = 1.0. 
! The integral to approximate: Integral (-1 <= x <= 1) F(x) dx
! The quadrature rule: Sum (1 <= i <= ORDER) weight(i)*F(xtab(i))
!****************************************************************************
!    Input, integer ORDER, the order of the rule.
!    ORDER must be greater than 0.
!    At current version order should be equal to 6 (it is enougth)
!    Output, double precision X(ORDER), the abscissas.
!    Output, double precision W(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
implicit none
integer, intent(IN) :: order
real(8) :: x(1:order), w(1:order)
integer :: i 
   if (order.eq.6) then 
      x(1) = 0.9324695142d0
      x(2) = 0.6612093864d0
      x(3) = 0.2386191861d0
      w(1) = 0.1713244924d0
      w(2) = 0.3607615730d0
      w(3) = 0.4679139346d0
   endif 
   do i =  1, Order / 2
      w(Order+1-i) = w(i)
      x(Order+1-i) = x(i)
      x(i) =-x(i)
   enddo
end subroutine Gauss_points


end module
    