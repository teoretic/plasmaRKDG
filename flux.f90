! ===========================================================
! РЕШЕНИЕ ОДНОМЕРНОЙ ЗАДАЧИ О РАСПАДЕ РАЗРЫВА и вычисление потоков
!
! - расчет физических и численных потоков для решения 1D задачи
!   на границах ячеек
! - расчет собственных чисел и векторов
! ===========================================================

! ==========================================================
! вектор консервативных переменных: U = (rho, rho*u, rho*v, rho*w, E)
! ==========================================================

! ==========================================================
! Используется консервативая форма записи уравнений вида
! dU/dt + dF/dr + dG/dz + alpha/r*H = 0
! Toro, с. 28, раздел 1.6.3, (1.102)-(1.104)
! ==========================================================


module flux

  use names

  implicit none

  integer, parameter:: alpha = 0 ! 0 -- (x,y) геометрия, 1 -- (r,z) геометрия

contains

  ! ======================================================
  ! вычисление физического потока, собственных чисел и векторов
  ! входные переменные: UU
  ! выходные переменные: F, D, R, L
  ! UU -- вектор консервативных переменных
  ! F -- вектор физического потока
  ! D -- диагональная матрица собственных чисел
  ! R -- матрица правых собственных векторов
  ! L -- матрица левых собственных веторов, L = R^{-1}
  ! Toro, с. 103, раздел 3.2.1

  pure subroutine FL(UU,F,D,R,L)
    implicit none

    ! входные переменные
    real (kind=precis), dimension(1:5), intent(in):: UU

    ! выходные переменные
    real (kind=precis), dimension(1:5), intent(out):: F
    real (kind=precis), dimension(1:5,1:5), intent(out):: D, R, L

    real (kind=precis):: rho, u, v, w, E, H2, a, p, k_e, mult
    real (kind=precis):: courant, verif

    ! вычисляем все простые переменные

    rho = UU(1)

    u = UU(2)/rho
    v = UU(3)/rho
    w = UU(4)/rho
    ! полная энергия (внутренняя плюс кинетическая)
    E = UU(5)
    k_e = u*u+v*v+w*w
    p = (E-0.5*rho*k_e)*(gmm-1.0)         ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!
    H2 = (E+p)/rho                                                                ! энтальпия

    verif = p*gmm/rho
    if (verif.LT.0.0) then
      verif = 0.0
    end if

    ! ошибка вычисления SQRT из-за деления на нуль в предыдущих формулах  !!!!
    a = SQRT(verif)       ! скорость звука


    F = (/ UU(2), UU(2)*u+p, UU(2)*v, UU(2)*w, u*(E+p) /)

    D = 0.0
    D(1,1) = u-a
    D(2,2) = u
    D(3,3) = u
    D(4,4) = u
    D(5,5) = u+a


    courant = D(5,5)*tau/0.2

    R(1,1) =1D0;    R(1,2) =0D0; R(1,3) =0D0; R(1,4) =1D0;     R(1,5) =1D0;
    R(2,1) =u-a;    R(2,2) =0D0; R(2,3) =0D0; R(2,4) =u;       R(2,5) =u+a;
    R(3,1) =v;      R(3,2) =1D0; R(3,3) =0D0; R(3,4) =v;       R(3,5) =v;
    R(4,1) =w;      R(4,2) =0D0; R(4,3) =1D0; R(4,4) =w;       R(4,5) =w;
    R(5,1) =H2-a*u; R(5,2) =v;   R(5,3) =w;   R(5,4) =0.5*k_e; R(5,5) =H2+a*u;

    mult = (gmm-1.0)/(2.0*a*a)

    L(1,1) =H2*mult+1.0/(a*2.0)*(u-a);  L(1,2) =-(u*mult+1.0/(a*2.0));   L(1,3) =-v*mult;    L(1,4) = -w*mult;       L(1,5) =1D0*mult;
    L(2,1) =-v;                         L(2,2) =0D0;                L(2,3) =1.0;        L(2,4) =0D0;            L(2,5) =0D0;
    L(3,1) =-w;                         L(3,2) =0D0;                L(3,3) =0D0;        L(3,4) = 1.0;      L(3,5) =0D0;
    L(4,1) =-2.0*H2*mult+2.0;           L(4,2) = 2.0*u*mult;        L(4,3) =2.0*v*mult; L(4,4) =2.0*w*mult;     L(4,5) =-2D0*mult;
    L(5,1) =H2*mult-1.0/(a*2.0)*(u+a);  L(5,2) =-u*mult+1.0/(a*2.0);L(5,3) =-v*mult;    L(5,4) =-w*mult;        L(5,5) =1D0*mult;

  end subroutine FL

  function FL_matrix(UU) result(Fl)
    real (kind=precis) :: UU(5), Fl(5,3)
    
    real (kind=precis) :: u, v, w, p, rho, E, k_e
    
    rho = UU(1)
    u = UU(2)/rho
    v = UU(3)/rho
    w = UU(4)/rho
    E = UU(5)
    k_e = u*u+v*v+w*w
    p = (E-0.5*rho*k_e)*(gmm-1.0)    ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!
    
    Fl(:,1) = (/ rho*u, rho*u*u+p, rho*u*v  , rho*u*w  , (E+p)*u /)
    Fl(:,2) = (/ rho*v, rho*u*v  , rho*v*v+p, rho*v*w  , (E+p)*v /)
    Fl(:,3) = (/ rho*w, rho*u*w  , rho*w*v  , rho*w*w+p, (E+p)*w /)
    
  end function
  
  
  pure subroutine FL_hllc(UU,F,u,v,w,a,p,rho,E)
    implicit none

    ! входные переменные
    real (kind=precis), dimension(1:5), intent(in):: UU

    ! выходные переменные
    real (kind=precis), dimension(1:5), intent(out):: F
    real (kind=precis), intent(out):: u, v, w, a, p, rho, E

    real (kind=precis):: H2, k_e, mult
    real (kind=precis):: courant

    ! вычисляем все простые переменные

    rho = UU(1)

    u = UU(2)/rho
    v = UU(3)/rho
    w = UU(4)/rho
    E = UU(5)
    k_e = u*u+v*v+w*w
    p = (E-0.5*rho*k_e)*(gmm-1.0)    ! давление  ИЗМЕНЕНИЕ !!!!!!!!!!!!!!!
    H2 = (E+p)/rho  ! энтальпия

    if (p*gmm/rho<0.0) then
      p = 0.0
    end if

    a = SQRT(p*gmm/rho)                                           ! скорость звука

    F(1)=UU(2)
    F(2)=UU(2)*u+p
    F(3)=UU(2)*v
    F(4)=UU(2)*w
    F(5)=u*(E+p)

  end subroutine FL_hllc

  ! ======================================================
  pure function flux_num_hllc(Ul,Ur) result(fr)
    implicit none

    real (kind=precis), dimension(1:5):: FR ! численный поток
    real (kind=precis), dimension(1:5), intent(in):: Ul, Ur
    real (kind=precis), dimension(1:5)::U_lr, Ur_et, Ul_et
    real (kind=precis), dimension(1:5):: F_l, F_r, F_lr
    real (kind=precis), dimension(1:5,1:5)::  R_lr, L_lr, D_lr

    real (kind=precis):: S_l, S_r, S_et, q_l, q_r

    real (kind=precis):: u_l,a_l,p_l,rho_l,E_l,v_l,w_l
    real (kind=precis):: u_r,a_r,p_r,rho_r,E_r,v_r,w_r
    real (kind=precis):: a_lr,p_lr,rho_lr,E_lr,v_lr,ulr, w_lr
    real (kind=precis):: u_et, p_et
    real (kind=precis):: numb1, numb2

    call FL_hllc(Ul,F_l,u_l,v_l,w_l,a_l,p_l,rho_l,E_l)
    call FL_hllc(Ur,F_r,u_r,v_r,w_r,a_r,p_r,rho_r,E_r)

    rho_lr = 0.5*(rho_l+rho_r)
    a_lr = 0.5*(a_l+a_r)
    p_et = 0.5*(p_l+p_r)+0.5*(u_l-u_r)*(rho_lr*a_lr)
    u_et = 0.5*(u_l+u_r)+0.5*(p_l-p_r)/(rho_lr*a_lr)

    if (p_et .LE. p_l) then
      q_l = 1.0
    else
      q_l = sqrt(1.0+(gmm+1.0)/(2.0*gmm)*(p_et/p_l-1.0))
    end if

    if (p_et .LE. p_r) then
      q_r = 1.0
    else
      q_r = sqrt(1.0+(gmm+1.0)/(2.0*gmm)*(p_et/p_r-1.0))
    end if

    S_l = u_l-a_l*q_l
    S_r = u_r+a_r*q_r
    S_et = u_et
    !S_et = (p_r-p_l+rho_l*u_l*(S_l-u_l)-rho_r*u_r*(S_r-u_r))/(rho_l*(S_l-u_l)-rho_r*(S_r-u_r))

    numb1 = rho_r*(S_r-u_r)/(S_r-S_et)
    numb2 = rho_l*(S_l-u_l)/(S_l-S_et)

    Ur_et(1) =numb1*1D0
    Ur_et(2) =numb1*S_et
    Ur_et(3) =numb1*v_r
    Ur_et(4) =numb1*w_r
    Ur_et(5) =numb1*(E_r/rho_r+(S_et-u_r)*(S_et+p_r/(rho_r*(S_r-u_r))))
    Ul_et(1) =numb2*1D0
    Ul_et(2) =numb2*S_et
    Ul_et(3) =numb2*v_l
    Ul_et(4) =numb2*w_l
    Ul_et(5) =numb2*(E_l/rho_l+(S_et-u_l)*(S_et+p_l/(rho_l*(S_l-u_l))))

    if (S_l .ge. 0D0) then
      FR = F_l
    end if

    if (S_l .le. 0D0 .and. S_et .ge. 0D0) then
      FR = F_l + S_l*(Ul_et-Ul)
    end if

    if (S_et .le. 0D0 .and. S_r .ge. 0D0) then
      FR = F_r + S_r*(Ur_et-Ur)
    end if

    if (S_r .le. 0D0) then
      FR= F_r
    end if

  end function flux_num_hllc

end module flux