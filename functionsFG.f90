module functionsFG
implicit none 

contains 

function Ffunc(P, W, dt, qin, mu, h) result(func)
implicit none 
real ::  qin
integer:: mu, h, dt 
real, dimension(2) :: P, W, func
    func(1) = W(1)/dt - qin/h - 12*mu / (0.5*(W(1)+W(2))**3) * (P(1)-P(2))/h**2 ! не должно быть деления на ноль
    func(2) = W(2)/dt + 0 + 12*mu / (0.5*(W(1)+W(2))**3) * (P(1)-P(2))/h**2
end function
    
function Gfunc(P, W, T, h, E) result(func)
implicit none 
integer:: h, E 
real, dimension(2,2) :: T
real, dimension(2) :: P, W, func
    func(1) = T(1,1)*P(1)+ T(1,2)*P(2) - W(1)
    func(2) = T(2,1)*P(1)+ T(2,2)*P(2) - W(2)
end function

function Ffunction(P, W, dt, hh, pi, qin, xx, matrAp, NN) result(func)
implicit none 
integer:: i, j, NN
real(8) :: hh, dt, qin, pi
real, dimension(NN) :: P, W, func, xx
real(8) :: matrAp(NN,NN)   
do i = 1, NN
    do j = 1, NN
        func(i) = W(i)/dt - matrAp(i,j)*P(j)
    end do
    if(i==1) then
        func(1) = func(1) - qin / (2.d0*pi) / xx(1) / hh
    end if
end do
end function
    
function Gfunction(P, W, matrPfromW, NN) result(func)
implicit none 
real, dimension(NN) :: P, W
real(8) :: matrPfromW(NN,NN)    
integer:: i, j, NN
real, dimension(NN) ::  func
do i = 1, NN
    func(i) = P(i)
    do j = 1, NN
        func(i) = func(i) - matrPfromW(i,j) * W(j)
    end do
end do
end function

subroutine RelaxMethod(P0, W0, N, dt, qin, hh, eps, relax, matrPfromW, xx, matrAp)
integer:: N, i
real(8) :: hh, dt, qin, pi, eps, differenceP, differenceW
real, dimension(2*N) :: relax !(1:N) - для F, (N+1:2N) - для G
real, dimension(N) :: P0, W0, F, G, Pn, Wn, xx
real(8) :: matrPfromW(N,N), matrAp(N, N)      
differenceP = 1.0
differenceW = 1.0 
do while(max(differenceP, differenceW) > eps)   
        F = Ffunction(P0, W0, dt, hh, pi, qin, xx, matrAp, N)
        G = Gfunction(P0, W0, matrPfromW, N)
        do i = 1, N
            Pn(i) = P0(i) - relax(i)*F(i)
            Wn(i) = W0(i) - relax(i+N)*G(i)
        end do
        differenceP = 0.0
        differenceW = 0.0
        do i = 1, N
            differenceP = max(differenceP, ABS(P0(i)-Pn(i)))
            differenceW = max(differenceW, ABS(W0(i)-Wn(i)))
        end do
        do i = 1, N
          P0(i) = Pn(i)
          W0(i) = Wn(i)
        end do
end do
end subroutine

subroutine initParam(eps, mu, E, qin, h, pout, dt, T, N)
real :: eps, qin, pout 
real, dimension(N/2, N/2) :: T
integer :: mu, E, h, dt, N
eps = 0.00001
mu=1
E=20
qin=0.05
h=1
pout = 0.01
dt = 1
T(1,1)= 2.0*(2.0*h*(1.0-1.0/16.0)**0.5/E)/3.0
T(1,2)= (2.0*h*(1.0-1.0/16.0)**0.5/E)/3.0
T(2,1)= (2.0*h*(1.0-9.0/16.0)**0.5/E)/3.0
T(2,2)= 2.0*(2.0*h*(1.0-9.0/16.0)**0.5/E)/3.0
end subroutine

function dFdx(Pk, Wk, NN) result(Aout) 
real, intent (IN) :: Pk(:), Wk(:)
integer, intent (IN) :: NN
real, dimension(NN,NN) :: Aout
integer:: i, j
do i = 1, NN
  do j = 1, NN
      if((i==1).and.(j==1)) then 
          Aout(i,j) = -24/(Wk(1)+Wk(2))**3
      else if((i==2).and.(j==2)) then 
          Aout(i,j) = Aout(1,1)
      else if ((i==1).and.(j==2)) then 
          Aout(i,j) = 24/(Wk(1)+Wk(2))**3
      else if ((i==2).and.(j==1)) then
          Aout(i,j) = Aout(1,2)
      else if ((i==1).and.(j==3)) then 
          Aout(i,j)=1 + 72*(Pk(1)-Pk(2))/(Wk(1)+Wk(2))**4
      else if ((i==2).and.(j==4)) then 
          Aout(i,j)=1 - 72*(Pk(1)-Pk(2))/(Wk(1)+Wk(2))**4
      else if ((i==1).and.(j==4)) then 
          Aout(i,j) = 72*(Pk(1)-Pk(2))/(Wk(1)+Wk(2))**4
      else if ((i==2).and.(j==3)) then 
          Aout(i,j) = - 72*(Pk(1)-Pk(2))/(Wk(1)+Wk(2))**4
      else if ((i==3).and.(j==1)) then 
          Aout(i,j) = 15**0.5/60.0
      else if ((i==3).and.(j==2)) then 
          Aout(i,j) = 15**0.5/120.0
      else if ((i==4).and.(j==1)) then 
          Aout(i,j) = 7**0.5/120.0
      else if ((i==4).and.(j==2)) then 
          Aout(i,j) = 7**0.5/60.0   
      else if ((i==3).and.(j==3)) then 
          Aout(i,j) = -1.0
      else if ((i==4).and.(j==4)) then
          Aout(i,j) = -1.0
      else 
          Aout(i,j) = 0.0
      end if
  end do
end do
end function

subroutine NewtonInit(P, W, N)
real, dimension(N/2) :: P, W
integer, intent (IN) :: N
integer :: i
do i = 1, N/2
    P(i)= 0.2
    W(i)= 0.01
end do
end subroutine

subroutine RelaxInit(P, W, N)
real, dimension(N/2) :: P, W
integer, intent (IN) :: N
integer :: i
do i = 1, N/2
    P(i)= 0.2
    W(i)= 0.01
end do
end subroutine

subroutine PrintMatrix(A,N)
real, intent (IN) :: A(:,:)      
integer, intent (IN) :: N
integer :: i, j
Print *, 'Matrix: '
do i = 1, N
    do j = 1,N
        Print *, A(i, j)
    end do
    Print *, '   '
end do    
end subroutine

subroutine PrintMultMatrix(A,B,N)
real, dimension(N,N) :: A
real(8), dimension(1:N,1:N) :: B, D
integer, intent (IN) :: N
integer :: i, j, k
Print *, 'Multiplication result : '
do i = 1, N
    do j = 1, N
        D(i, j) = 0
        do k = 1, N
        D(i, j) = D(i, j) + A(i, k) * B(k, j)
        end do
        PRINT*, D(i,j)
    end do
    Print*, ' '
end do   
end subroutine
    
function invertMatrix(Ain, NN) result(Aout) 
real, intent (IN) :: Ain(:,:)      
integer, intent (IN) :: NN
real(8), dimension(1:NN,1:NN) :: A, Aout
integer:: INFO, i, j
integer, dimension(NN) :: IPIV
do i = 1, NN ! единичная матрица
  do j = 1, NN
      if((i>j).or.(j>i)) then
      Aout(i,j)=0.0
      else if(i==j) then
          Aout(i,j)= 1.0
      end if
      A(i,j)= Ain(i,j)
      end do
end do
CALL DGESV(NN, NN, A, NN, IPIV, Aout, NN, INFO)
end function 

function ConvertMatrix(Ain, NN) result(A) ! По заданным элементам, расположенным на трех диагоналях, строим матрицу
real, intent (IN) :: Ain(:,:)      
integer, intent (IN) :: NN
real(8), dimension(1:NN,1:NN) :: A
integer:: i, j
do i = 1, NN 
  do j = 1, NN
      A(i, j)= 0.0
      A(i, i-1) = Ain(1,i)
      A(i, i) = Ain(2,i)
      A(i, i+1) = Ain(3,i)   
  end do
end do
A(1, 1) = Ain(2,1)
A(1, 2) = Ain(3,1)
A(NN, NN-1) = Ain(1,NN)
A(NN, NN) = Ain(2,NN)
end function 

subroutine NewtonStart(P0, W0, N, dt, qin, mu, h, T, E, eps)
real ::  qin, eps, differenceP, differenceW
integer:: mu, h, dt, E, N, i
real, dimension(N/2, N/2) :: T
real, dimension(N/2) :: P0, W0, Fn, Gn, Pk, Wk
real(8), dimension(1:N, 1:N) :: Ainvert
real, dimension(N, N) :: A
differenceP = 1.0
differenceW = 1.0 
do while(max(differenceP, differenceW) > eps)    
        PRINT*, P0(1), P0(2), W0(1), W0(2)
        Fn = Ffunc(P0, W0, dt, qin, mu, h)
        Gn = Gfunc(P0, W0, T, h, E)
        A = dFdx(P0, W0, N)
        Ainvert = invertMatrix(A, N)
        do i = 1, 2
          Pk(i) = P0(i) - Ainvert(i,1)*Fn(1)-Ainvert(i,2)*Fn(2)-Ainvert(i,3)*Gn(1)-Ainvert(i,4)*Gn(2)
          Wk(i) = W0(i) - Ainvert(i+2,1)*Fn(1)-Ainvert(i+2,2)*Fn(2)-Ainvert(i+2,3)*Gn(1)-Ainvert(i+2,4)*Gn(2)
        end do
        differenceP = ABS(Pk(1)-P0(1))
        differenceW = ABS(Wk(1)-W0(1))
        do i = 1, 2
          P0(i) = Pk(i)
          W0(i) = Wk(i)
        end do
end do
PRINT*, P0(1), P0(2), W0(1), W0(2)
end subroutine

subroutine RelaxStart(P1, W1, N, dt, qin, mu, h, T, E, eps, relax)
real ::  qin, eps, differenceP, differenceW
integer:: mu, h, dt, E, N, k, i
real, dimension(N/2, N/2) :: T, relax
real, dimension(N/2) :: P1, W1, Fr, Gr, Pn, Wn
real(8), dimension(1:N, 1:N) :: Ainvert
real, dimension(N, N) :: A
differenceP = 1.0
differenceW = 1.0 
!do while(max(differenceP, differenceW) >= eps)
do k = 1,70
        PRINT*, P1(1), P1(2), W1(1), W1(2)
        Fr = Ffunc(P1, W1, dt, qin, mu, h)
        Gr = Gfunc(P1, W1, T, h, E)
        Pn(1) = P1(1) - relax(1,1)*Gr(1)
        Pn(2) = P1(2) - relax(1,2)*Gr(2)
        Wn(1) = W1(1) - relax(2,1)*Fr(2)
        Wn(2) = W1(2) - relax(2,2)*Fr(1)
        differenceP = max(ABS(P1(1)-Pn(1)), ABS(P1(2)-Pn(2)))
        differenceW = max(ABS(W1(1)-Wn(1)), ABS(W1(2)-Wn(2)))
        do i = 1, 2
          P1(i) = Pn(i)
          W1(i) = Wn(i)
        end do
end do
PRINT*, Pn(1), Pn(2), Wn(1), Wn(2)
end subroutine

! вывод в файл графика двумерной функции на равномерной ортогональной сетке
subroutine Grafik2D(X,Y,Func,NN,filename)       
implicit none 
real, intent (IN) :: Func(:,:)                ! процедуре не обязательно знать максимальный размер, но используемый размер нужен
real, intent (IN) :: X(:), Y(:)      
integer, intent (IN) :: NN
character(*), intent (IN) :: filename 
integer :: i,j 
open(unit=48,file=filename)
!write(48, '(A)') 'Variables = X, Y, F'                      ! не знаю, как оформить аналог для gnuplot 
!write(48, '(A,I4,A,I4)') 'zone I = ', NN, ' J = ', NN       ! но как gnuplot знает как точки соединять в сетку?
do i = 1, NN
    do j = 1, NN 
      if (Func(i,j) < 0) then
        write (48,'(F8.6,F9.6,E12.4)') X(j), Y(i), Func(i,j)
      else
        write (48,'(F8.6,F9.6,E11.4)') X(j), Y(i), Func(i,j)
      end if
          
        !write (48,'(3(E12.4))') X(j), Y(i), Func(i,j)
        ! в фортране есть форматная запись, что позволяет регулировать "количество пробелов"
        ! например здесь "E14.4" - значит, что числа выводятся в экспоненциальном формате "   -0.9984E+00"
        ! всего число занимает 14 символов и имеет 4 знака после точки, получчается 3 пробела 
        ! если поставить меньше, чем 14, то и пробелов будет меньше. 
        ! если поствить "F" вместо "E", то будт запись с фиксированной точкой, посмотрите, что получится. 
        ! 3 - в "3(E14.4)" значит,что таких чисел будет 3 в строе
    end do
    write (48,*)
end do
close(48)
end subroutine 

end module 
 