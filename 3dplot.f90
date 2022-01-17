PROGRAM demo
use functionsFG                     ! функции удобно выносить в модули, тогда файл меньше
use elasticRadial
implicit none                       ! это запрещает использование неописанных переменных 
integer, parameter :: Nmax = 4 ! позволяет менять размеры сразу у всех массивов
real(8), dimension(Nmax/2) :: P, W, Wn, P0, W0, P1, W1
real(8), dimension(Nmax) :: relax  
real(8) :: eps, pout, hh, lambda
integer :: i
integer :: NN           ! size of system of nonlinear equations
integer :: NN05         ! size of subvectors w, and p
! параметры для новых функций
real(8) :: ElasticCoef, Rfrac       ! elastic and mesh parameters    !TODO: make common subroutine for all parameters
real(8) :: matrWfromP(Nmax/2,Nmax/2), matrPfromW(Nmax/2,Nmax/2)          ! matrix for w_k = T^s_k*p_s
type(TfluidParams) :: fluidParams   ! parameters for fluid
real(8) :: matrAp(3,Nmax/2), matrAp2(3,Nmax/2), matrAp3(3,Nmax/2)       ! fluid matrix [3,1..NN] with boundary conditions
real(8) :: RHS(1,1:Nmax)              ! Right hand side [1,1..NN] with boundary conditions
real(8) :: Xcentr(1:Nmax/2), Xbound(0:Nmax/2)               ! mesh, cell centers [1..NN] and cell boundaries [0;Rfrac], [0..NN]
real    :: timeNewton, timeRelax, timePrev, t(2), etime, timeLevenbergMarkvardt ! Нужно ли здесь сделать t(3)?

NN = Nmax       ! размер всей системы уравнений
NN05 = Nmax/2   ! размер каждого вектора w, p

! subroutines to check quality of fluid and elastic matrixes
!call testElasticMatrixRadial()     
!call testFluidMatrixRadial()    
! set parameters (!!!!! p[MPa], w[mm], mu[MPa*s])
ElasticCoef = 8.d0/3.14159265d0/20.d0*(1-0.25d0**2)
call setFluidParamsForTest(fluidParams)         
Rfrac = 4.d0   
fluidParams%pout = 0.01d0
fluidParams%mu = 1.d0/1.d6  ! mu[MPa*s]
fluidParams%dt = 100.d0
Xbound(0:NN05) = (/ (i, i = 0, NN05) /) * Rfrac/NN05
!Xcentr(1:NN05) = (/ (i, i = 1, NN05) /) * Rfrac/NN05 + (Rfrac/NN05)/2
do i = 1, NN05  ;   Xcentr(i) = 0.5d0*(Xbound(i-1)+Xbound(i))    ;   enddo
Wn = 0.d0
! initial step, to make fluid matrix 
call makeElasticMatrixRadial(Xbound,Rfrac,NN05, matrWfromP)    
matrPfromW = invertMatrix(matrWfromP(1:NN05,1:NN05), NN05)  
!call PrintMatrix(matrPfromW,NN05,NN05)
pout = 0.01d0
P = pout
W(1:NN05) = matmul(matrWfromP(1:NN05,1:NN05),P(1:NN05))*ElasticCoef
call makeFluidMatrixAndRHSRadial(Xcentr,NN05,fluidParams,W,Wn,matrAp,RHS)
!TODO: новые функции считаются как
!Ffunc(1) = W/dt - matrAp*P - fluidParams%qin / (2.d0*pi) / Xcentr(1) / hh
!Ffunc(остальные) = W/dt - matrAp*P 
!Ffunc(N) = matrAp*P - p_out = 0
! матрица matrAp записана как трёхдиагональная, надо конвертировать 
! Ax=b => matrAp(1,i)*x(i-1)+matrAp(2,i)*x(i)+matrAp(3,i)*x(i+1) = b(i)
! Gfunc = W - matrWfromP * P или Gfunc = P - matrPfromW * W
! для обобщения будет удобнее с matrPfromW, но можно начать с любого варианта
! matrPfromW, matrWfromP считаются один раз и записаны как обычные матрицы

eps = 1.d-8
hh = Xcentr(2)-Xcentr(1)
do i = 1, Nmax/2
    print *, Xcentr(i)
end do
print *, hh
!call PrintMatrix(matrPfromW, NN05, NN05)
P0(1:NN05) = 1.d-3  ! начальное приближение
W0(1:NN05) = 1.d-2
PRINT*,"Newton's method:"

matrAp2(1:3,1:NN05) = matrAp(1:3,1:NN05)  ! matrAp после метода Ньютона возвращается видоизмененной, поэтому запоминаю ее первоначальный вид в matrAp2
matrAp3(1:3,1:NN05) = matrAp(1:3,1:NN05)

timePrev = etime( t )
call NewtonMethod(P0(1:NN05), W0(1:NN05), NN, fluidParams%dt, fluidParams%qin, pi, hh, eps, matrPfromW(1:NN05,1:NN05), Xcentr(1:NN05), matrAp(1:3,1:NN05), fluidParams)
timeNewton = etime( t ) - timePrev
call Grafik1D(Xcentr,P0,NN05,'NewtonP.plt')
call Grafik1D(Xcentr,W0,NN05,'NewtonW.plt')

PRINT*,"Relaxation method :"
P0(1:NN05) = 1.d-3  ! начальное приближение
W0(1:NN05) = 1.d-2
relax(1:NN) = 0.01d0
timePrev = etime( t )
call RelaxMethod(P0(1:NN05), W0(1:NN05), NN, fluidParams%dt, fluidParams%qin, pi, hh, eps, relax(1:NN), matrPfromW(1:NN05,1:NN05), Xcentr(1:NN05), matrAp2(1:3,1:NN05), fluidParams)
timeRelax = etime( t ) - timePrev

call Grafik1D(Xcentr,P0,NN05,'RelaxP.plt')
call Grafik1D(Xcentr,W0,NN05,'RelaxW.plt')

PRINT*,"LevenbergMarkvardts method :"
P0(1:NN05) = 1.d-3  ! начальное приближение
W0(1:NN05) = 1.d-2
lambda = 0.00001d0 ! этот параметр нужно брать как можно меньшим, тогда точность лучше и сходимость быстрее
timePrev = etime( t )
call MethodLevenbergMarkvardt(P0(1:NN05), W0(1:NN05), NN, fluidParams%dt, fluidParams%qin, pi, hh, eps, lambda, matrPfromW(1:NN05,1:NN05), Xcentr(1:NN05), matrAp3(1:3,1:NN05), fluidParams)
timeLevenbergMarkvardt = etime( t ) - timePrev

call Grafik1D(Xcentr,P0,NN05,'LevenbergMarkvardtP.plt')
call Grafik1D(Xcentr,W0,NN05,'LevenbergMarkvardtW.plt')

! Кажется, я нашла причину, почему метод релаксации и метод Ньютона выдавали разные результаты. Система имеет несколько решений и если взять начальные приближения не достаточно близкими к некоторому корню системы, то методы выдают разные корни. Я взяла начальное приближение близкое к тому корню, что выдал метод релаксации и метод Ньютона сошелся к тому же корню, что и метод релаксации.
! Я нашла онлайн решение системы для N=2 (т.е. система из 4 уравнений) и у нее есть два решения. Думаю, аналогично и для больших N.
    
write(*,'(A,F8.3,A)') 'time Newton = ', timeNewton, ' s'
write(*,'(A,F8.3,A)') 'time Relax  = ', timeRelax , ' s'
write(*,'(A,F8.3,A)') 'time LevenbergMarkvardt  = ', timeLevenbergMarkvardt , ' s'
END PROGRAM demo
