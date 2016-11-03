! сервисные функции

module servicefunc

use names

implicit none

contains

pure function vect_mult2d (a,b) result (c)

  real (kind=precis),dimension(1:2), Intent(In) :: a,b
  real (kind=precis) :: c

  c = a(1)*b(2)-a(2)*b(1)

end function vect_mult2d


pure function vect_mult (a,b) result (c)

  real (kind=precis),dimension(1:3), Intent(In) :: a,b
  real (kind=precis),dimension(1:3) :: c

  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)

end function vect_mult

pure function TriangleVol(a,b) result (S)
  real (kind=precis), dimension(1:2), Intent(In):: a,b
  real (kind=precis) :: S

  S = 0.5*abs(a(1)*b(2)-b(1)*a(2))
end function TriangleVol

pure function TriangleVol3(a,b,c) result (S)
  real (kind=precis), dimension(1:3), Intent(In):: a,b,c
  real (kind=precis), dimension(1:3) :: d,e,f
  real (kind=precis) :: S

  d = b-a
  e = c-a
  f(1) = d(2)*e(3)-e(2)*d(3)
  f(2) = d(3)*e(1)-e(3)*d(1)
  f(3) = d(1)*e(2)-e(1)*d(2)
  S = 0.5*norm2_3(f)
end function TriangleVol3

pure function distance(a,b) result (c)
  real (kind=precis), dimension(1:2), Intent(In):: a,b
  real (kind=precis) :: c

  c = sqrt(sum((b-a)**2))
end function distance

pure function distance3(a,b) result (c)
  real (kind=precis), dimension(1:3), Intent(In):: a,b
  real (kind=precis) :: c

  c = sqrt(sum((b-a)**2))
end function distance3

function inverseM3(A,stat) result(AI)
  real(kind=precis) :: A(3,3)
  real(kind=precis) :: AI(3,3)
  real(kind=precis) :: znam
  integer :: ii
  logical :: stat
  
  znam = -A(1,3)*A(2,2)*A(3,1) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) + A(1,1)*A(2,2)*A(3,3)
  
  AI(1,:) = (/ -A(2,3)*A(3,2)+A(2,2)*A(3,3),  A(1,3)*A(3,2)-A(1,2)*A(3,3), -A(1,3)*A(2,2)+A(1,2)*A(2,3) /)
  AI(2,:) = (/  A(2,3)*A(3,1)-A(2,1)*A(3,3), -A(1,3)*A(3,1)+A(1,1)*A(3,3),  A(1,3)*A(2,1)-A(1,1)*A(2,3) /)
  AI(3,:) = (/ -A(2,2)*A(3,1)+A(2,1)*A(3,2),  A(1,2)*A(3,1)-A(1,1)*A(3,2), -A(1,2)*A(2,1)+A(1,1)*A(2,2) /)
  AI = AI/znam
  stat = .true.
  
  if (abs(znam)<1D-10) then
    stat = .false.
!     print *, "bida-bida"
!     print *, znam
!     do ii=1,3
!       print *, A(ii,:)
!     end do
!     read*
  end if
  
end function inverseM3



pure function signum(numb) result(sig)

  integer, Intent(In) :: numb
  integer :: sig

  sig = 1
  if (numb<0) sig = -1

end function signum

pure function signumr(numb) result(sig)

  real(kind=precis), Intent(In) :: numb
  real(kind=precis) :: sig

  sig = 1.
  if (numb<0._precis) sig = -1.

end function signumr

pure function norm2(a) result (c)
  real (kind=precis), dimension(1:2), Intent(In):: a
  real (kind=precis) :: c

  c = sqrt(sum(a**2))
end function norm2

pure function norm2_3(a) result (c)
  real (kind=precis), dimension(1:3), Intent(In):: a
  real (kind=precis) :: c

  c = sqrt(sum(a**2))
end function norm2_3

pure function pintr (a,t1,t2,t3) result (pin)
  real (kind=precis), dimension(1:2), Intent(In):: a,t1,t2,t3
  logical :: pin

  pin = (signumr(vect_mult2d(a-t1,t2-t1))<1E-14 .and. signumr(vect_mult2d(a-t2,t3-t2))<1E-14 .and. signumr(vect_mult2d(a-t3,t1-t3))<1E-14)
end function


  subroutine T_rot(n, T, Ti)
    implicit none

    real (kind=precis), dimension(1:2), intent(in):: n ! n = (cos(theta), sin(theta)) -- вектор нормали

    real (kind=precis), dimension(1:5,1:5), intent(out):: T ! матрица преобразования
    real (kind=precis), dimension(1:5,1:5), intent(out):: Ti ! обратная к ней

        T = 0.0

        T(1,1)=1.0
        T(2,2)=n(1)
        T(2,3)=n(2)
        T(3,2)=-n(2)
        T(3,3)=n(1)
        T(4,4)=1.0
        T(5,5)=1.0

        Ti = transpose(T)

  end subroutine T_rot


  FUNCTION i2c(n) RESULT(str)
    IMPLICIT NONE

    INTEGER n
    INTEGER, PARAMETER:: r = 8 ! количество разрядов -- соотвествует 10000000

    CHARACTER(LEN=r):: str


    INTEGER, DIMENSION(0:r-1):: digit = 0
    INTEGER i
    INTEGER:: foo

    IF( (n>10**r-1).OR.(n<0)) THEN
       str = ''
    ELSE

       digit(r-1) = n/10**(r-1)
       foo = n - digit(r-1)*10**(r-1)

       DO i=r-2,0,-1
          digit(i) = foo/10**i
          foo = foo - digit(i)*10**i
       ENDDO

       str = ''

       str = ACHAR(IACHAR('0') + digit(7)) // &
             ACHAR(IACHAR('0') + digit(6)) // &
             ACHAR(IACHAR('0') + digit(5)) // &
             ACHAR(IACHAR('0') + digit(4)) // &
             ACHAR(IACHAR('0') + digit(3)) // &
             ACHAR(IACHAR('0') + digit(2)) // &
             ACHAR(IACHAR('0') + digit(1)) // &
             ACHAR(IACHAR('0') + digit(0))
    ENDIF

  END FUNCTION i2c


end module servicefunc