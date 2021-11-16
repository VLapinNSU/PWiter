PROGRAM demo
use functionsFG                     ! функции удобно выновить в модули, тогда файл меньше
use elasticRadial
use fluidRadial
implicit none                       ! это запрещает использование неописанных переменных 
integer, parameter :: Nmax = 4    ! позволяет менять размеры сразу у всех массивов
! real(8) - двойная точность. real - одинарная. обычно используют двойную
real, dimension(Nmax,Nmax) :: F 
real, dimension(Nmax,Nmax) :: G  
real, dimension(Nmax) :: P, W, Wn
real, dimension(2,2) :: T, relax  
real, dimension(2) :: P0, W0, P1, W1
real :: eps, qin, pout
integer :: i, mu, E, h, dt  
integer :: NN          ! используемый размер массива не обязан совпадать с самим размером. Позволяет менять его без перекомпиляции программы

! параметры для новых функций
real(8) :: ElasticCoef, Rfrac       ! elastic and mesh parameters    !TODO: make common subroutine for all parameters
real(8) :: matrWfromP(Nmax,Nmax), matrPfromW(Nmax,Nmax)          ! matrix for w_k = T^s_k*p_s
type(TfluidParams) :: fluidParams   ! parameters for fluid
real(8) :: matrAp(3,Nmax)           ! fluid matrix [3,1..NN] with boundary conditions
real(8) :: RHS(1,1:Nmax)              ! Right hand side [1,1..NN] with boundary conditions
real(8) :: XX(1:Nmax)               ! mesh, cell boundaries [0;Rfrac], [0..NN]

NN = Nmax/2 ! только для новых функций
! subroutines to check quality of fluid and elastic matrixes
!call testElasticMatrixRadial()     
!call testFluidMatrixRadial()    
! set parameters (!!!!! p[MPa], w[mm], mu[MPa*s])
ElasticCoef = 8.d0/3.14159265d0/20.d0*(1-0.25d0**2)
Rfrac = 10.d0   
call setFluidParamsForTest(fluidParams)         
fluidParams%pout = 0.01d0   ;   fluidParams%mu = 1.d0/1.d6
xx = (/ (i, i = 0, NN) /) * Rfrac/NN    ! build mesh
Wn = 0.d0
! initial step, to make fluid matrix 
call makeElasticMatrixRadial(xx,Rfrac,NN, matrWfromP)    
!matrPfromW = invertMatrix(matrWfromP(1:NN,1:NN), NN)   ! У меня выдает ошибку, говорит, что что-то не так с типом matrWfromP
pout = 0.01d0
P = pout
W(1:NN) = matmul(matrWfromP(1:NN,1:NN),P(1:NN))*ElasticCoef
!call makeFluidMatrixAndRHSRadial(xx,NN,fluidParams,W,Wn,matrAp,RHS)  ! У меня еще вот эта функция не работает (говорит, что тип W и Wn не тот, что в функции)
!TODO: новые функции считаются как 
!Ffunc(1) = W/dt - matrAp*P - fluidParams%qin / (2.d0*pi) / xx(1) / hh
!Ffunc(остальные) = W/dt - matrAp*P 
! матрица matrAp записана как трёхдиагональная, надо конвертировать 
! Ax=b => matrAp(1,i)*x(i-1)+matrAp(2,i)*x(i)+matrAp(3,i)*x(i+1) = b(i)
! Gfunc = W - matrWfromP * P или Gfunc = P - matrPfromW * W
! для обобщения будет удобнее с matrPfromW, но можно начать с любого варианта
! matrPfromW, matrWfromP считаются один раз и записаны как обычные матрицы
    
NN = Nmax
call initParam(eps, mu, E, qin, h, pout, dt, T, NN) !задание параметров

!Метод Ньютона:
call NewtonInit(P0, W0, NN) !начальное приближение
PRINT*,"Newton's method:"
call NewtonStart(P0, W0, NN, dt, qin, mu, h, T, E, eps)

!Метод релаксации:
relax(1,1) = 0.001
relax(1,2) = -0.001
relax(2,1) = -0.001
relax(2,2) = 0.001
call RelaxInit(P1, W1, NN) !начальное приближение
PRINT*,"Relaxation method :"
call RelaxStart(P1, W1, NN, dt, qin, mu, h, T, E, eps, relax)

END PROGRAM demo
