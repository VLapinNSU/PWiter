PROGRAM demo
use functionsFG                     ! функции удобно выновить в модули, тогда файл меньше
use elasticRadial
use fluidRadial
implicit none                       ! это запрещает использование неописанных переменных 
integer, parameter :: Nmax = 4    ! позволяет менять размеры сразу у всех массивов
! real(8) - двойная точность. real - одинарная. обычно используют двойную
real(8), dimension(1:Nmax,1:Nmax) :: Ainvert, D
real(8), dimension(Nmax,Nmax) :: F,A ! если не используете описание переменных в комментариях, то можно писать через запятую 
real(8), dimension(Nmax,Nmax) :: G  
real(8), dimension(Nmax) :: P, W     ! как здесь
real(8), dimension(2,2) :: T, relax  
real(8), dimension(2) :: P0, W0, Pk, Wk, Pn, P1, Wn, W1, Fn, Gn, Fr, Gr
integer, dimension(Nmax) :: IPIV
real(8) :: eps, differenceP, differenceW, qin, pout    ! но лучше комментарии какие-то делать
integer :: i, j, iter, l, mu, E, h, dt, INFO, k          ! надо все переменные описывать (мы же запретили implicit none)
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
    matrPfromW = invertMatrix(matrWfromP(1:NN,1:NN), NN)            !TODO: у меня обращение матрицы работает. Не могу воспроизвести ошибку
    pout = 0.01d0
    P = pout
    W(1:NN) = matmul(matrWfromP(1:NN,1:NN),P(1:NN))*ElasticCoef
    call makeFluidMatrixAndRHSRadial(xx,NN,fluidParams,W,Wn, matrAp,RHS)
    !TODO: новые функции считаются как 
    !Ffunc(1) = W/dt - matrAp*P - fluidParams%qin / (2.d0*pi) / xx(1) / hh
    !Ffunc(остальные) = W/dt - matrAp*P 
    ! матрица matrAp записана как трёхдиагональная, надо конвертировать 
    ! Ax=b => matrAp(1,i)*x(i-1)+matrAp(2,i)*x(i)+matrAp(3,i)*x(i+1) = b(i)
    ! Gfunc = W - matrWfromP * P или Gfunc = P - matrPfromW * W
    ! для обобщения будет удобнее с matrPfromW, но можно начать с любого варианта
    ! matrPfromW, matrWfromP считаются один раз и записаны как обычные матрицы
    
!TODO: задание параметров в отдельную subroutine initParam
NN = Nmax
eps = 0.001
mu=1
E=20
qin=0.05
h=1
pout = 0.01
dt = 1
T(1,1)= 2/3*(2*h*(1-1/16)**0.5/E)
T(1,2)= 1/3*(2*h*(1-1/16)**0.5/E)
T(2,1)= 1/3*(2*h*(1-9/16)**0.5/E)
T(2,2)= 2/3*(2*h*(1-9/16)**0.5/E)
P = (/ (i, i = 0, NN-1) /) * 1.0/100.0        ! в фортране можно немного короче 
W = (/ (i, i = 0, NN-1) /) * 1.0/100.0

! TODO: этот текст что делает? до --------------
! рисование графика лучше убрать в отдельную процедуру (тоже положил)
! ну и уж точно не использовать явный вид функций в основной программе
do i = 1, NN
  do j = 1, NN
        if(ABS(W(i)) >= eps) then
           !Ffunc(P, W, dt, qin, mu, h) = 0.0
        else 
          !F(P, W, dt, qin, mu, h) = 0.0
        end if
        !Gfunc(P, W, T, h, E) = 0.0
  end do
end do 
!call Grafik2D(P,W,F,NN,'F_grafik.dat')       
!call Grafik2D(P,W,G,NN,'G_grafik.dat')    
!TODO: -------------- до этого места?

!TODO: Задание начального приближения вынести в отдельную subroutine NewtonInit
!Метод Ньютона:
P0(1) = 0.3 ! начальное приближение P0, W0
P0(2) = 0.2
W0(1) = 0.02
W0(2) = 0.02
differenceP = 0.0
differenceW = 0.0     
A = dFdx(P0, W0, NN)
Ainvert = invertMatrix(A, NN)
do i = 1, NN
          do j = 1, NN
              PRINT*, A(i,j)
          end do
          PRINT*, '  '
end do 
!TODO: печать матрицы - в отдельную subroutine 
do i = 1, NN
          do j = 1, NN
              PRINT*, Ainvert(i,j)
          end do
          PRINT*, '  '
end do ! произведение матриц
!TODO: произведение матриц - в отдельную subroutine multMatrAB
do i = 1, NN
    do j = 1, NN
    D(i, j) = 0
        do k = 1, NN
        D(i, j) = D(i, j) + A(i, k) * Ainvert(k, j)
        end do
        PRINT*, D(i,j)
    end do
    PRINT*, '  '
end do

!TODO: сделал здесь форматирование отступов, постарайтесь так по всему тексту сделать
!TODO: Метод Ньютона в отдельную subroutine NewtonStart
PRINT*,"Newton's method:"
!do while(max(differenceP, differenceW) > eps)             !Нужно выбрать норму, в которой посчитаем Pk-P0 и Wk-W0 (взяла пока максимум)
    do k= 1,10       
        PRINT*, P0(1), P0(2), W0(1), W0(2)
        Fn = Ffunc(P0, W0, dt, qin, mu, h)
        Gn = Gfunc(P0, W0, T, h, E)
        A = dFdx(P0, W0, NN)
        Ainvert = invertMatrix(A, NN)
        !TODO: печать матрицы - в отдельную subroutine, не надо снова писать код
        do i = 1, NN
            do j = 1, NN
                PRINT*, Ainvert(i,j)
            end do
        end do 
        do i = 1, 2
            Pk(i) = P0(i) - Ainvert(i,1)*Fn(1)-Ainvert(i,2)*Fn(2)-Ainvert(i,3)*Gn(1)-Ainvert(i,4)*Gn(2)
            Wk(i) = W0(i) - Ainvert(i+1,1)*Fn(1)-Ainvert(i+1,2)*Fn(2)-Ainvert(i+1,3)*Gn(1)-Ainvert(i+1,4)*Gn(2)
        end do
        differenceP = ABS(Pk(1)-P0(1))
        differenceW = ABS(Wk(1)-W0(1))
        do i = 1, 2
            P0(i) = Pk(i)
            W0(i) = Wk(i)
        end do
    end do
PRINT*, P0(1), P0(2), W0(1), W0(2)

!TODO: Задание начального приближения вынести в отдельную subroutine RelaxInit
!Метод релаксации:
    relax(1,1) = 0.01
    relax(1,2) = -0.01
    relax(2,1) = 0.01
    relax(2,2) = 0.01
    do i = 1, 2 ! начальное приближение P1, W1
        P1(i) = 0.25
        W1(i) = 0.01
    end do
    differenceP = 1.0
    differenceW = 1.0
    l= 0.0 ! счетчик итераций    !TODO: счетчик итераций не integer???

!TODO: метод релаксации в отдельную subroutine RelaxStart
    PRINT*,"Relaxation method :"
    do while(max(differenceP, differenceW) >= eps)
        l = l+1
        !PRINT*, P1(1), P1(2), W1(1), W1(2)
        Fr = Ffunc(P1, W1, dt, qin, mu, h)
        Gr = Gfunc(P1, W1, T, h, E)
        Pn(1) = P1(1) - relax(1,1)*Fr(1)
        Pn(2) = P1(2) - relax(1,2)*Fr(2)
        Wn(1) = W1(1) - relax(2,1)*Gr(1)
        Wn(2) = W1(2) - relax(2,2)*Gr(2)
        differenceP = max(ABS(P1(1)-Pn(1)), ABS(P1(2)-Pn(2)))
        differenceW = max(ABS(W1(1)-Wn(1)), ABS(W1(2)-Wn(2)))
        do i = 1, 2
          P1(i) = Pn(i)
          W1(i) = Wn(i)
             end do
end do
PRINT*, Pn(1), Pn(2), Wn(1), Wn(2)
PRINT*, l
PRINT*, 'Hello, world' 
END PROGRAM demo