module functionsFG
use fluidRadial
implicit none 

contains 

function Ffunction(P, W, WprevTimeStep, dt, hh, pi, qin, pout, xx, matrAp, NN) result(func)
implicit none 
integer:: i, j, NN
real(8) :: hh, dt, qin, pi, pout
real(8), dimension(NN) :: P, W, func, xx, WprevTimeStep
real(8) :: matrAp(NN,NN)   

    do i = 1, NN
        func(i) = W(i)/dt - WprevTimeStep(i) / dt
        if (i==NN) func(i) = -pout       ! the last equation is P-p_out=0, no dW/dt, matrAp(NN,NN)=-1
        do j = 1, NN
            func(i) = func(i) - matrAp(i,j)*P(j)
        end do
    end do
    func(1) = func(1) - qin / (2.d0*pi) / xx(1) / hh
    !print *, qin / (2.d0*pi) / xx(1) / hh
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

subroutine RelaxMethod(P0, W0, N, dt, qin, pout, pi, hh, eps, relax, matrPfromW, xx, matrAp, fluidParams)
type(TfluidParams),intent(IN) :: fluidParams   ! parameters for fluid
integer:: N, i, k
real(8) :: hh, dt, qin, pout, pi, eps, differenceP, differenceW
real(8), dimension(N) :: relax !(1:N/2) - для F, (N/2+1:N) - для G
real(8), dimension(N/2) :: P0, W0, F, G, Pn, Wn, xx
real(8) :: matrPfromW(N/2,N/2), matrAp(3, N/2), matrApconvert(N/2, N/2)        
integer :: Iter
real(8) :: WprevTimeStep(1:N/2)     ! width distribution at previous time step
real(8) :: RHS(1,1:N/2)               ! right hand side of fluid equation. is not used here
!real(8) :: xBound(0:N/2)            ! coordinates of cell boundaries (not centr)

!xBound(1:N/2) = xx(1:N/2)+0.5d0*(xx(2)-xx(1))   ;   xBound(0) = xx(1)-0.5d0*(xx(2)-xx(1))     ! temporary calculate boundaries here (do it at mash forming)
differenceP = 1.0
differenceW = 1.0 
WprevTimeStep = 0.d0
open(10,file = 'RelaxHist.plt')
write(10,'(A)') 'Variables = Iter, difP, difW'
Iter = 0
do while(max(differenceP, differenceW) > eps)  
    Iter = Iter + 1
    !do i = 1, N/2
    !    Print*, P0(i)
    !end do
    !do i = 1, N/2
    !    Print*, W0(i)
    !end do
    !Print*, '    '
    call makeFluidMatrixAndRHSRadial(xx,N/2,fluidParams,W0,WprevTimeStep,matrAp,RHS)
    matrApconvert(1: N/2, 1: N/2) = ConvertMatrix(matrAp, N/2)  !конвертирую матрицу matrAp
    !call PrintMatrix(matrApconvert, N/2, N/2)
    F = Ffunction(P0, W0, WprevTimeStep, dt, hh, pi, qin,pout, xx, matrApconvert, N/2)
    G = Gfunction(P0, W0, matrPfromW, N/2)
    do i = 1, N/2
        Pn(i) = P0(i) - relax(i)*G(i)
        Wn(i) = W0(i) - relax(i+N/2)*F(i)
    end do
    differenceP = maxval( ABS(P0(1:N/2)-Pn(1:N/2)) )    
    differenceW = maxval( ABS(W0(1:N/2)-Wn(1:N/2)) )
    P0(1:N/2) = Pn(1:N/2)
    W0(1:N/2) = Wn(1:N/2)
    write(10,'(I6, 2(E16.6))') Iter, differenceP, differenceW
end do   
close(10)
do i = 1, N/2
    Print*, P0(i)
end do
do i = 1, N/2
    Print*, W0(i)
end do
Print*, Iter, '= Relaxs iterations'
end subroutine

function dFdxForNewton(matrAp, matrPfromW, dt, N, xx, fluidParams, P, W, hh) result(Aout) 
type(TfluidParams),intent(IN) :: fluidParams   ! parameters for fluid
integer, intent (IN) :: N
real(8) :: dt, hh
real(8) :: matrPfromW(N,N), matrAp(N, N) 
real(8), dimension(2*N,2*N) :: Aout
real(8) :: xx(1:N)          ! in makeFluidMatrixAndRHSRadial cell centres are used
real(8), dimension(N) :: P, W
integer:: i, j
Aout = 0.d0
!(dF/dp)
do i = 1, N
    do j = 1, N
        Aout(i,j) = - matrAp(i,j)
    enddo
enddo
!Aout(N,N) = 1.d0                            ! last F equation is p-p_out = 0 => dF/dp =1
!(dF/dw)
Aout(1, N+1) = ( (0.5d0*3.d0*xx(1  )*W(1  )**2) / (12.d0*fluidParams%mu*hh**2*xx(1)) ) * (P(1)-P(2)) 
Aout(1, N+2) = ( (0.5d0*3.d0*xx(2  )*W(2  )**2) / (12.d0*fluidParams%mu*hh**2*xx(1)) ) * (P(1)-P(2))
Aout(1, N+1) = Aout(1,N+1) + 1.d0/dt
do i = 2, N-1
    Aout(i, N+i-1) = ( (0.5d0*3.d0*xx(i-1)*W(i-1)**2) / (12.d0*fluidParams%mu*hh**2*xx(i)) ) * (P(i)-P(i-1))
    Aout(i, N+i+1) = ( (0.5d0*3.d0*xx(i+1)*W(i+1)**2) / (12.d0*fluidParams%mu*hh**2*xx(i)) ) * (P(i)-P(i+1))
    Aout(i, N+i  ) = ( (0.5d0*3.d0*xx(i  )*W(i  )**2) / (12.d0*fluidParams%mu*hh**2*xx(i)) ) * (2.d0*P(i)-P(i-1)-P(i+1)) 
    Aout(i, N+i  ) = Aout(i,N+i) + 1.d0/dt
enddo
!Aout(N,N+N) = Aout(N,N+N) + 1/dt            ! для N не надо, последнее уравнение p - p_out = 0
!(dG/dp)
do j = 1, N
    Aout(j+N,j) = 1.d0
enddo
!(dG/dw)
do j = 1, N
do i = 1, N
    Aout(i+N,j+N) = - matrPfromW(i,j)
enddo
enddo
end function

subroutine NewtonMethod(P0, W0, N, dt, qin, pout, pi, hh, eps, matrPfromW, xx, matrAp, fluidParams,WprevTimeStep)
type(TfluidParams),intent(IN) :: fluidParams   ! parameters for fluid
integer:: N, i, j, Iter
real(8) :: hh, dt, pi, qin, pout, eps, differenceP, differenceW
real(8), dimension(N/2) :: P0, W0, F, G, Pn, Wn, xx
real(8) :: matrPfromW(N/2,N/2), matrAp(3, N/2), matrApconvert(N/2, N/2)  
real(8), dimension(1:N, 1:N) :: Ainvert
real(8), dimension(N, N) :: A
real(8), intent(IN) :: WprevTimeStep(1:N/2)     ! width distribution at previous time step
real(8) :: RHS(1,1:N/2)               ! right hand side of fluid equation. is not used here
!real(8) :: xBound(0:N/2)            ! coordinates of cell boundaries (not centr)

!xBound(1:N/2) = xx(1:N/2)+0.5d0*(xx(2)-xx(1))   ;   xBound(0) = xx(1)-0.5d0*(xx(2)-xx(1))     ! temporary calculate boundaries here (do it at mash forming)
differenceP = 1.0
differenceW = 1.0 
!WprevTimeStep = 0.d0
open(10,file = 'NewtonHist.plt')
!write(10,'(A)') 'Variables = Iter, difP, difW'
Iter = 0
do while(max(differenceP, differenceW) > eps)  
    Iter = Iter + 1
    !do i = 1, N/2
    !    Print*, P0(i), ' P0 '
    !end do
    !do i = 1, N/2
    !    Print*, W0(i), ' W0 '
    !end do
    !Print*, '    '
    call makeFluidMatrixAndRHSRadial(xx,N/2,fluidParams,W0,WprevTimeStep,matrAp,RHS)
    matrApconvert(1: N/2, 1: N/2) = ConvertMatrix(matrAp, N/2)  !конвертирую матрицу matrAp
    F = Ffunction(P0, W0,WprevTimeStep, dt, hh, pi, qin, pout, xx, matrApconvert, N/2)
    G = Gfunction(P0, W0, matrPfromW, N/2)
    A = dFdxForNewton(matrApconvert, matrPfromW, dt, N/2, xx, fluidParams, P0, W0, hh)
    Ainvert(1:N,1:N) = invertMatrix(A, N)
    !call PrintMatrix(matrPfromW,N/2,N/2)
    !call PrintMatrix(A,N,N)
    !call PrintMatrix(Ainvert,N,N)
    !call PrintMultMatrix(A,Ainvert,N)
    do i = 1, N/2
        Pn(i) = P0(i)
        Wn(i) = W0(i)
        do j = 1, N/2
            Pn(i) = Pn(i) - (Ainvert(i,j)*F(j) + Ainvert(i,j+N/2)*G(j))
            Wn(i) = Wn(i) - (Ainvert(i+N/2,j)*F(j) + Ainvert(i+N/2,j+N/2)*G(j))
        end do
    end do
    differenceP = maxval( ABS(P0(1:N/2)-Pn(1:N/2)) )    
    differenceW = maxval( ABS(W0(1:N/2)-Wn(1:N/2)) )
    P0(1:N/2) = Pn(1:N/2)
    W0(1:N/2) = Wn(1:N/2)
    write(10,'(I6, 2(E16.6))') Iter, differenceP, differenceW
end do
close(10)
do i = 1, N/2
    Print*, P0(i), ' P0 '
end do
do i = 1, N/2
    Print*, W0(i), ' W0 '
end do
Print*, Iter, '= Newtons iterations'
end subroutine

subroutine MethodLevenbergMarkvardt(P0, W0, N, dt, qin, pout, pi, hh, eps, lambda, matrPfromW, xx, matrAp, fluidParams)
type(TfluidParams),intent(IN) :: fluidParams   ! parameters for fluid
integer:: N, i, k, Iter, nu, r, l
real(8) :: hh, dt, pi, qin,pout, eps, differenceP, differenceW, lambda, E, En
real(8), dimension(N/2) :: P0, W0, F, G, Pn, Wn, xx, Fn, Gn
real(8) :: matrPfromW(N/2,N/2), matrAp(3, N/2), matrApconvert(N/2, N/2)  
real(8), dimension(1:N, 1:N) :: Ainvert
real(8), dimension(N, N) :: J, JT, A, B ! JT = J^T, A = (J^T)*J+lambda*I, B = ((J^T)*J+lambda*I)^(-1)*(J^T)
real(8) :: WprevTimeStep(1:N/2)     ! width distribution at previous time step
real(8) :: RHS(1,1:N/2)               ! right hand side of fluid equation. is not used here
differenceP = 1.d0
differenceW = 1.d0
E = 1.d0
En = 2.d0
WprevTimeStep = 0.d0
open(10,file = 'LevenbergMarkvardtHist.plt')
write(10,'(A)') 'Variables = Iter, difP, difW'
Iter = 0
nu = 10
r = -1
do while ((E > eps).and.(Iter < 1000)) 
!do while ((max(differenceP, differenceW) > eps).and.(Iter < 10000000))  
    Iter = Iter + 1
    call makeFluidMatrixAndRHSRadial(xx,N/2,fluidParams,W0,WprevTimeStep,matrAp,RHS)
    matrApconvert(1: N/2, 1: N/2) = ConvertMatrix(matrAp, N/2)  !конвертирую матрицу matrAp
    F = Ffunction(P0, W0,WprevTimeStep, dt, hh, pi, qin,pout, xx, matrApconvert, N/2)
    G = Gfunction(P0, W0, matrPfromW, N/2)
    J(1:N, 1:N) = dFdxForNewton(matrApconvert, matrPfromW, dt, N/2, xx, fluidParams, P0, W0, hh)
    JT(1:N, 1:N) = transpose(J(1:N, 1:N)) ! JT = J^T, транспонирую матрицу J
    A(1:N, 1:N) = matmul(JT(1:N, 1:N), J(1:N, 1:N))   
    do i = 1, N
        !A(i,i) = A(i,i) + lambda        ! A = (J^T)*J+lambda*I 
        A(i,i) = A(i,i)*(1 + lambda)     ! A = (J^T)*J+lambda*Diag[(J^T)*J] - так лучше, учитываются скачки
    end do    
    Ainvert(1:N,1:N) = invertMatrix(A(1:N, 1:N), N)
    B(1:N,1:N) = matmul(Ainvert(1:N,1:N), JT(1:N,1:N))  ! B = ((J^T)*J+lambda*I)^(-1)*(J^T)
    do i = 1, N/2 
        Pn(i) = P0(i)
        Wn(i) = W0(i)
        do k = 1, N/2
            Pn(i) = Pn(i) - (B(i,k)*F(k) + B(i,k+N/2)*G(k))
            Wn(i) = Wn(i) - (B(i+N/2,k)*F(k) + B(i+N/2,k+N/2)*G(k))
            !Pn(i) = Pn(i) - (1/lambda)*2*(JT(i,k)*F(k) + JT(i,k+N/2)*G(k)) ! метод градиентого спуска
            !Wn(i) = Wn(i) - (1/lambda)*2*(JT(i+N/2,k)*F(k) + JT(i+N/2,k+N/2)*G(k))
        end do
    end do
    call makeFluidMatrixAndRHSRadial(xx,N/2,fluidParams,Wn,WprevTimeStep,matrAp,RHS)
    matrApconvert(1: N/2, 1: N/2) = ConvertMatrix(matrAp, N/2)  !конвертирую матрицу matrAp
    Fn = Ffunction(Pn, Wn, WprevTimeStep, dt, hh, pi, qin,pout, xx, matrApconvert, N/2)
    Gn = Gfunction(Pn, Wn, matrPfromW, N/2)
    En=0.d0
    E=0.d0
    do i = 1, N/2 
        En = En + Fn(i)*Fn(i) + Gn(i)*Gn(i) ! берем такую фукцию ошибки
        E = E + F(i)*F(i) + G(i)*G(i)
    end do
    !print*, 'En = ', En, 'E = ',  E
    if (En < E) then
        P0(1:N/2) = Pn(1:N/2)
        W0(1:N/2) = Wn(1:N/2)  
        r = -1
        lambda = lambda*(1.d0/nu)!lambda_{k} = lambda_{k-1}*nu^r - меняем параметр
        print*, 'here'
    else  
        r = r + 1
        l = r
        do while (l > 0)
            lambda = lambda*nu !lambda_{k} = lambda_{k-1}*nu^r - меняем параметр
            l = l - 1
        end do
    end if
    print *, 'l = ', lambda
    differenceP = maxval( ABS(P0(1:N/2)-Pn(1:N/2)) )    
    differenceW = maxval( ABS(W0(1:N/2)-Wn(1:N/2)) )
    write(10,'(I6, 2(E16.6))') Iter, differenceP, differenceW
end do
close(10)
do i = 1, N/2
    Print*, P0(i), ' P0 '
end do
do i = 1, N/2
    Print*, W0(i), ' W0 '
end do
Print*, Iter, '= LevenbergMarkvardts iterations'
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


! вывод в файл график одномерной функции 
subroutine Grafik1D(X,Func,NN,filename)       
implicit none 
real(8), intent (IN) :: Func(:)                ! процедуре не обязательно знать максимальный размер, но используемый размер нужен
real(8), intent (IN) :: X(:)
integer, intent (IN) :: NN
character(*), intent (IN) :: filename 
integer :: i
open(unit=48,file=filename)
write(48, '(A)') 'Variables = X, F'                      ! не знаю, как оформить аналог для gnuplot 
do i = 1, NN
    write (48,'(2(E14.6))') X(i), Func(i)
end do
close(48)
end subroutine 
end module 