module functionsFG
implicit none 

contains 

function Ffunc(P, W, dt, qin, mu, h) result(func)
implicit none 
real ::  qin
integer:: mu, h, dt 
real, dimension(2) :: P, W, func
    func(1) = W(1)/dt - qin/h - 12*mu / (0.5*(W(1)+W(2))**3) * (P(1)-P(2))/h**2 ! �� ������ ���� ������� �� ����
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
    
function invertMatrix(Ain, NN) result(Aout) 
real, intent (IN) :: Ain(:,:)      
integer, intent (IN) :: NN
real(8), dimension(1:NN,1:NN) :: A, Aout
integer:: INFO, i, j
integer, dimension(NN) :: IPIV
do i = 1, NN ! ��������� �������
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

! ����� � ���� ������� ��������� ������� �� ����������� ������������� �����
subroutine Grafik2D(X,Y,Func,NN,filename)       
implicit none 
real, intent (IN) :: Func(:,:)                ! ��������� �� ����������� ����� ������������ ������, �� ������������ ������ �����
real, intent (IN) :: X(:), Y(:)      
integer, intent (IN) :: NN
character(*), intent (IN) :: filename 
integer :: i,j 
open(unit=48,file=filename)
!write(48, '(A)') 'Variables = X, Y, F'                      ! �� ����, ��� �������� ������ ��� gnuplot 
!write(48, '(A,I4,A,I4)') 'zone I = ', NN, ' J = ', NN       ! �� ��� gnuplot ����� ��� ����� ��������� � �����?
do i = 1, NN
    do j = 1, NN 
      if (Func(i,j) < 0) then
        write (48,'(F8.6,F9.6,E12.4)') X(j), Y(i), Func(i,j)
      else
        write (48,'(F8.6,F9.6,E11.4)') X(j), Y(i), Func(i,j)
      end if
          
        !write (48,'(3(E12.4))') X(j), Y(i), Func(i,j)
        ! � �������� ���� ��������� ������, ��� ��������� ������������ "���������� ��������"
        ! �������� ����� "E14.4" - ������, ��� ����� ��������� � ���������������� ������� "   -0.9984E+00"
        ! ����� ����� �������� 14 �������� � ����� 4 ����� ����� �����, ����������� 3 ������� 
        ! ���� ��������� ������, ��� 14, �� � �������� ����� ������. 
        ! ���� �������� "F" ������ "E", �� ���� ������ � ������������� ������, ����������, ��� ���������. 
        ! 3 - � "3(E14.4)" ������,��� ����� ����� ����� 3 � �����
    end do
    write (48,*)
end do
close(48)
end subroutine 

end module 