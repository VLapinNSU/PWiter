PROGRAM demo
use functionsFG                     ! функции удобно выновить в модули, тогда файл меньше
use elasticRadial
use fluidRadial
implicit none                       ! это запрещает использование неописанных переменных 
integer, parameter :: Nmax = 4    ! позволяет менять размеры сразу у всех массивов
! real(8) - двойная точность. real - одинарная. обычно используют двойную
real(8), dimension(Nmax,Nmax) :: F, G, A
real(8), dimension(2*Nmax, 2*Nmax) :: B
real(8), dimension(1:2*Nmax, 1:2*Nmax) :: C
real(8), dimension(Nmax) :: P, W, Wn, P0, W0, P1, W1
real(8), dimension(2*Nmax) :: relax  
real(8) :: eps, pout, hh
integer :: i
integer :: NN           ! size of system of nonlinear equations
integer :: NN05         ! size of subvectors w, and p
! параметры для новых функций
real(8) :: ElasticCoef, Rfrac       ! elastic and mesh parameters    !TODO: make common subroutine for all parameters
real(8) :: matrWfromP(Nmax,Nmax), matrPfromW(Nmax,Nmax)          ! matrix for w_k = T^s_k*p_s
type(TfluidParams) :: fluidParams   ! parameters for fluid
real(8) :: matrAp(3,Nmax)           ! fluid matrix [3,1..NN] with boundary conditions
real(8) :: RHS(1,1:Nmax)              ! Right hand side [1,1..NN] with boundary conditions
real(8) :: XX(1:Nmax)               ! mesh, cell boundaries [0;Rfrac], [0..NN]

NN05 = Nmax/2 ! только для новых функций
! subroutines to check quality of fluid and elastic matrixes
!call testElasticMatrixRadial()     
!call testFluidMatrixRadial()    
! set parameters (!!!!! p[MPa], w[mm], mu[MPa*s])
ElasticCoef = 8.d0/3.14159265d0/20.d0*(1-0.25d0**2)
Rfrac = 10.d0   
call setFluidParamsForTest(fluidParams)         
fluidParams%pout = 0.01d0   ;   fluidParams%mu = 1.d0/1.d6
xx = (/ (i, i = 0, NN05) /) * Rfrac/NN05    ! build mesh
Wn = 0.d0
! initial step, to make fluid matrix 
call makeElasticMatrixRadial(xx,Rfrac,NN05, matrWfromP)    
matrPfromW = invertMatrix(matrWfromP(1:NN05,1:NN05), NN05)  
pout = 0.01d0
P = pout
W(1:NN05) = matmul(matrWfromP(1:NN05,1:NN05),P(1:NN05))*ElasticCoef
call makeFluidMatrixAndRHSRadial(xx,NN05,fluidParams,W,Wn,matrAp,RHS)
!TODO: новые функции считаются как 
!Ffunc(1) = W/dt - matrAp*P - fluidParams%qin / (2.d0*pi) / xx(1) / hh
!Ffunc(остальные) = W/dt - matrAp*P 
! матрица matrAp записана как трёхдиагональная, надо конвертировать 
! Ax=b => matrAp(1,i)*x(i-1)+matrAp(2,i)*x(i)+matrAp(3,i)*x(i+1) = b(i)
! Gfunc = W - matrWfromP * P или Gfunc = P - matrPfromW * W
! для обобщения будет удобнее с matrPfromW, но можно начать с любого варианта
! matrPfromW, matrWfromP считаются один раз и записаны как обычные матрицы

eps = 0.001d0
hh = xx(2)-xx(1) ! не нашла, где взять значение
NN = Nmax
do i = 1, NN  ! начальное приближение
    P0(i) = 0.0d0
    W0(i) = 0.0d0
    relax(i) = 0.01d0
    relax(i+NN) = 0.01d0
end do  

A = ConvertMatrix(matrAp, NN05) ! сначала конвертирую матрицу matrAp
B = dFdxForNewton(A, matrPfromW, fluidParams%dt, NN05) !Считаю матрицу производных для метода Ньютона
call PrintMatrix(B,2*NN05,2*NN05) !Вывожу ее
C = invertMatrix(B, 2*NN05) !Ищу обратную к ней
call PrintMatrix(C,2*NN05, 2*NN05) !Но он почему - то считает обратной единичную матрицу, то есть что-то не так с функцией InvertMatrix
call PrintMultMatrix(B,C,2*NN05) ! Вот их произведение, видно, что матрица получилась не единичной
PRINT*,"Newton's method:"
call NewtonMethod(P0, W0, NN, fluidParams%dt, fluidParams%qin, pi, hh, eps, matrPfromW, xx, A)
PRINT*,"Relaxation method :"
call RelaxMethod(P0, W0, NN, fluidParams%dt, fluidParams%qin, pi, hh, eps, relax, matrPfromW, xx, A)

END PROGRAM demo
