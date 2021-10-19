PROGRAM demo
use functionsFG                     ! ������� ������ �������� � ������, ����� ���� ������
implicit none                       ! ��� ��������� ������������� ����������� ���������� 
integer, parameter :: Nmax = 4    ! ��������� ������ ������� ����� � ���� ��������
! real(8) - ������� ��������. real - ���������. ������ ���������� �������
real(8), dimension(1:Nmax,1:Nmax) :: Ainvert, D
real, dimension(Nmax,Nmax) :: F,A ! ���� �� ����������� �������� ���������� � ������������, �� ����� ������ ����� ������� 
real, dimension(Nmax,Nmax) :: G  
real, dimension(Nmax) :: P, W     ! ��� �����
real, dimension(2,2) :: T, relax  
real, dimension(2) :: P0, W0, Pk, Wk, Pn, P1, Wn, W1, Fn, Gn, Fr, Gr
integer, dimension(Nmax) :: IPIV
real :: eps, differenceP, differenceW, qin, pout    ! �� ����� ����������� �����-�� ������
integer :: i, j, iter, l, mu, E, h, dt, INFO, k          ! ���� ��� ���������� ��������� (�� �� ��������� implicit none)
integer :: NN          ! ������������ ������ ������� �� ������ ��������� � ����� ��������. ��������� ������ ��� ��� �������������� ���������
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
P = (/ (i, i = 0, NN-1) /) * 1.0/100.0        ! � �������� ����� ������� ������ 
W = (/ (i, i = 0, NN-1) /) * 1.0/100.0

! ��������� ������� ����� ������ � ��������� ��������� (���� �������)
! �� � �� ����� �� ������������ ����� ��� ������� � �������� ���������
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

!����� �������:
P0(1) = 0.3 ! ��������� ����������� P0, W0
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
do i = 1, NN
          do j = 1, NN
              PRINT*, Ainvert(i,j)
          end do
          PRINT*, '  '
end do ! ������������ ������
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

PRINT*,"Newton's method:"
!do while(max(differenceP, differenceW) > eps)             !����� ������� �����, � ������� ��������� Pk-P0 � Wk-W0 (����� ���� ��������)
 do k= 1,10       
        PRINT*, P0(1), P0(2), W0(1), W0(2)
        Fn = Ffunc(P0, W0, dt, qin, mu, h)
        Gn = Gfunc(P0, W0, T, h, E)
        A = dFdx(P0, W0, NN)
        Ainvert = invertMatrix(A, NN)
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
!����� ����������:
relax(1,1) = 0.01
relax(1,2) = -0.01
relax(2,1) = 0.01
relax(2,2) = 0.01
do i = 1, 2 ! ��������� ����������� P1, W1
P1(i) = 0.25
W1(i) = 0.01
end do
differenceP = 1.0
differenceW = 1.0
l= 0.0 ! ������� ��������
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