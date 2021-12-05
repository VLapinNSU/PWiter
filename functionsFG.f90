﻿module functionsFG
implicit none 

contains 

function Ffunction(P, W, dt, hh, pi, qin, xx, matrAp, NN) result(func)
implicit none 
integer:: i, j, NN
real(8) :: hh, dt, qin, pi
real(8), dimension(NN) :: P, W, func, xx
real(8) :: matrAp(NN,NN)   
do i = 1, NN
    func(i) = W(i)/dt
    do j = 1, NN
        func(i) = func(i) - matrAp(i,j)*P(j)
    end do
    if(i==1) then
        func(1) = func(1) - qin / (2.d0*pi) / xx(1) / hh
        !print*, qin / (2.d0*pi) / xx(1) / hh
    end if
end do
end function
    
function Gfunction(P, W, matrPfromW, NN) result(func)
implicit none 
real(8), dimension(NN) :: P, W
real(8) :: matrPfromW(NN,NN)    
integer:: i, j, NN
real(8), dimension(NN) ::  func
do i = 1, NN
    func(i) = P(i)
    do j = 1, NN
        func(i) = func(i) - matrPfromW(i,j) * W(j)
    end do
end do
end function

subroutine RelaxMethod(P0, W0, N, dt, qin, pi, hh, eps, relax, matrPfromW, xx, matrAp)
integer:: N, i
real(8) :: hh, dt, qin, pi, eps, differenceP, differenceW
real(8), dimension(N) :: relax !(1:N/2) - для F, (N/2+1:N) - для G
real(8), dimension(N/2) :: P0, W0, F, G, Pn, Wn, xx
real(8) :: matrPfromW(N/2,N/2), matrAp(3, N/2), matrApconvert(N/2, N/2)        
differenceP = 1.0
differenceW = 1.0 
matrApconvert(1: N/2,1: N/2) = ConvertMatrix(matrAp, N/2)
do while(max(differenceP, differenceW) > eps)  
        do i = 1, N/2
            Print*, P0(i)
        end do
        do i = 1, N/2
            Print*, W0(i)
        end do
        Print*, '    '
        F = Ffunction(P0, W0, dt, hh, pi, qin, xx, matrApconvert, N/2)
        G = Gfunction(P0, W0, matrPfromW, N/2)
        do i = 1, N/2
            Pn(i) = P0(i) - relax(i)*G(i)
            Wn(i) = W0(i) - relax(i+N/2)*F(i)
        end do
        differenceP = 0.0
        differenceW = 0.0
        do i = 1, N/2
            differenceP = max(differenceP, ABS(P0(i)-Pn(i)))
            differenceW = max(differenceW, ABS(W0(i)-Wn(i)))
        end do
        do i = 1, N/2
          P0(i) = Pn(i)
          W0(i) = Wn(i)
        end do
end do   
do i = 1, N/2
            Print*, P0(i)
        end do
        do i = 1, N/2
            Print*, W0(i)
        end do
        Print*, '    '
end subroutine

function dFdxForNewton(matrAp, matrPfromW, dt, N) result(Aout) 
integer, intent (IN) :: N
real(8) :: dt
real(8) :: matrPfromW(N,N), matrAp(N, N) 
real(8), dimension(2*N,2*N) :: Aout
integer:: i, j
Aout = 0.d0
do i = 1, 2*N
  do j = 1, 2*N
      if (i <= N) then
          if (j <= N) then
                Aout(i,j) = - matrAp(i,j)
          else
              Aout(i,j) = 1/dt
          end if    
      else
          if (j <= N) then
              Aout(i,j) = 1 
          else
                Aout(i,j) = - matrPfromW(i-N,j-N)
          end if  
      end if    
  end do
end do
end function

subroutine NewtonMethod(P0, W0, N, dt, qin, pi, hh, eps, matrPfromW, xx, matrAp)
integer:: N, i, j, Iter
real(8) :: hh, dt, pi, qin, eps, differenceP, differenceW
real(8), dimension(N/2) :: P0, W0, F, G, Pn, Wn, xx
real(8) :: matrPfromW(N/2,N/2), matrAp(3, N/2), matrApconvert(N/2, N/2)  
real(8), dimension(1:N, 1:N) :: Ainvert
real(8), dimension(N, N) :: A
differenceP = 1.0
differenceW = 1.0 
matrApconvert(1: N/2, 1: N/2) = ConvertMatrix(matrAp, N/2)  !конвертирую матрицу matrAp
!do while(max(differenceP, differenceW) > eps)   
do Iter = 1, 10
    do i = 1, N/2
        Print*, P0(i), ' P0 '
    end do
    do i = 1, N/2
        Print*, W0(i), ' W0 '
    end do
    Print*, '    '
    F = Ffunction(P0, W0, dt, hh, pi, qin, xx, matrApconvert, N/2)
    G = Gfunction(P0, W0, matrPfromW, N/2)
    A = dFdxForNewton(matrApconvert, matrPfromW, dt, N/2)  ! Матрица получается не зависящей от P0 и W0, т.е. постоянной
    Ainvert(1:N,1:N) = invertMatrix(A, N) ! Тоже постоянная матрица
    !call PrintMatrix(A,N,N)
    !call PrintMatrix(Ainvert,N,N)
    !call PrintMultMatrix(A,Ainvert,N)
    do i = 1, N/2
        Pn(i) = P0(i)
        Wn(i) = W0(i)
        do j = 1, N/2
            Pn(i) = Pn(i) - Ainvert(i,j)*F(j) - Ainvert(i,j+N/2)*G(j)
            Wn(i) = Wn(i) - Ainvert(i+N/2,j)*F(j) - Ainvert(i+N/2,j+N/2)*G(j)
        end do
    end do
    differenceP = 0.0
    differenceW = 0.0
    do i = 1, N/2
       differenceP = max(differenceP, ABS(P0(i)-Pn(i)))
       differenceW = max(differenceW, ABS(W0(i)-Wn(i)))
    end do
    do i = 1, N/2
        P0(i) = Pn(i)
        W0(i) = Wn(i)
    end do
end do
do i = 1, N/2
    Print*, P0(i), ' P0 '
end do
do i = 1, N/2
    Print*, W0(i), ' W0 '
end do
Print*, '    '
end subroutine
    
function invertMatrix(Ain, NN) result(Aout) 
real(8), intent (IN) :: Ain(:,:)      
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
real(8), intent (IN) :: Ain(3,NN)      
integer, intent (IN) :: NN
real(8), dimension(NN,NN) :: A
integer:: i, j
do i = 1, NN 
  do j = 1, NN
      A(i, j) = 0.0
      if ((i>1).and.(i<NN)) then
        A(i, i-1) = Ain(1,i)
        A(i, i) = Ain(2,i)
        A(i, i+1) = Ain(3,i) 
      end if
  end do
end do
A(1, 1) = Ain(2,1)
A(1, 2) = Ain(3,1)
A(NN, NN-1) = Ain(1,NN)
A(NN, NN) = Ain(2,NN)
end function 

subroutine PrintMatrix(A,N,M)
real(8), intent (IN) :: A(:,:)      
integer, intent (IN) :: N, M
integer :: i, j
Print *, 'Matrix: '
do i = 1, N
    do j = 1,M
        Print *, A(i, j)
    end do
    Print *, '   '
end do    
end subroutine

subroutine PrintMultMatrix(A,B,N)
real(8), dimension(N,N) :: A
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

function Ffunc(P, W, dt, qin, mu, h) result(func)
implicit none 
real(8) ::  qin
integer:: mu, h, dt 
real(8), dimension(2) :: P, W, func
    func(1) = W(1)/dt - qin/h - 12*mu / (0.5*(W(1)+W(2))**3) * (P(1)-P(2))/h**2 ! не должно быть деления на ноль
    func(2) = W(2)/dt + 0 + 12*mu / (0.5*(W(1)+W(2))**3) * (P(1)-P(2))/h**2
end function
    
function Gfunc(P, W, T, h, E) result(func)
implicit none 
integer:: h, E 
real(8), dimension(2,2) :: T
real(8), dimension(2) :: P, W, func
    func(1) = T(1,1)*P(1)+ T(1,2)*P(2) - W(1)
    func(2) = T(2,1)*P(1)+ T(2,2)*P(2) - W(2)
end function

subroutine initParam(eps, mu, E, qin, h, pout, dt, T, N)
real(8) :: eps, qin, pout 
real(8), dimension(N/2, N/2) :: T
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
real(8), intent (IN) :: Pk(:), Wk(:)
integer, intent (IN) :: NN
real(8), dimension(NN,NN) :: Aout
integer:: i, j
Aout = 0.d0
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
real(8), dimension(N/2) :: P, W
integer, intent (IN) :: N
integer :: i
do i = 1, N/2
    P(i)= 0.2
    W(i)= 0.01
end do
end subroutine

subroutine RelaxInit(P, W, N)
real(8), dimension(N/2) :: P, W
integer, intent (IN) :: N
integer :: i
do i = 1, N/2
    P(i)= 0.2
    W(i)= 0.01
end do
end subroutine

subroutine NewtonStart(P0, W0, N, dt, qin, mu, h, T, E, eps)
real(8) ::  qin, eps, differenceP, differenceW
integer:: mu, h, dt, E, N, i
real(8), dimension(N/2, N/2) :: T
real(8), dimension(N/2) :: P0, W0, Fn, Gn, Pk, Wk
real(8), dimension(1:N, 1:N) :: Ainvert
real(8), dimension(N, N) :: A
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
real(8) ::  qin, eps, differenceP, differenceW
integer:: mu, h, dt, E, N, k, i
real(8), dimension(N/2, N/2) :: T, relax
real(8), dimension(N/2) :: P1, W1, Fr, Gr, Pn, Wn
real(8), dimension(1:N, 1:N) :: Ainvert
real(8), dimension(N, N) :: A
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
real(8), intent (IN) :: Func(:,:)                ! процедуре не обязательно знать максимальный размер, но используемый размер нужен
real(8), intent (IN) :: X(:), Y(:)      
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