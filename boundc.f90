! ==========================================================
! Содержит процедyру и данные, необходимые для задания
! граничных и начальных условий
! ==========================================================

module boundc
  use names
  use servicefunc
  use basisRKDG
!  use flux

  implicit none

contains

subroutine update_bc(solArr)
  implicit none

  
  real (kind=precis), Intent(InOut) :: solArr(:,:,:,:)
  real (kind=precis) :: V(1:3) = 0.0, normal(1:3) = 0.0,rho, pressure,vvel
  integer:: elem,i,surf,point
  real (kind=precis), dimension(1:3,1:3) :: T, Ti
!   
  real (kind=precis) :: cosz, cosy, sinz, siny, xyproek, soundspeed, dist
! 
!   real (kind=precis), dimension(1:3) :: surf_cent
!     ! выходные переменные
!     real (kind=precis), dimension(1:5):: F
!     real (kind=precis) :: uu, vv, ww, aa, pp, rho, EE, ion, pressure, temper
!     logical :: stat
!   
!   
  do elem=1,nElems_b
    i = bc_elems(elem)
    surf = bc_surf(elem)
!                 ! U(:,nElems + elem) = (rho, rho*u, rho*v, rho*w, E)
    select case (bc_type(elem))
      case (1)
!                !       втекание
    rho = 88.6853D0
    pressure = 5529D3
    vvel = 1480D0

        do point=1,trianGaussN
          solArr(1,point,2,surf) = rho
          solArr(2,point,2,surf) = rho*vvel
          solArr(3,point,2,surf) = 0D0
          solArr(4,point,2,surf) = 0D0
          solArr(5,point,2,surf) = pressure/(gmm-1D0)+0.5*rho*(vvel**2)
        end do

      case (2)
!                          ! стенка
        do point=1,trianGaussN
          V = solArr(2:4,point,1,surf)
          normal(:) = surf_normal(:,surf)

          xyproek = sqrt(Sum(normal(1:2)**2))
          if (xyproek > 1e-7) then
            cosz = normal(1)/xyproek
            sinz = normal(2)/xyproek
          else
            cosz = 1D0
            sinz = 0D0
          endif

          cosy = xyproek
          siny = -normal(3)

          !поворот

          T(1,1:3) = (/ cosy*cosz,  -sinz, cosz*siny /)
          T(2,1:3) = (/ sinz*cosy,   cosz, sinz*siny /)
          T(3,1:3) = (/     -siny,    0D0,      cosy /)

          Ti = Transpose(T)
          V = matmul(Ti,V)

          V(1) = -V(1)

          ! обратный поворот
          V = matmul(T,V)
          solArr(1,point,2,surf) = solArr(1,point,1,surf)
          solArr(2:4,point,2,surf) = V
          solArr(5,point,2,surf) = solArr(5,point,1,surf)
        end do

      case (3)
         ! вытекание
        do point=1,trianGaussN
          solArr(:,point,2,surf) = solArr(:,point,1,surf)
        end do
        
!       case (4)
!          ! исторические условия
!         U_til(:,nElems + elem) = U(:,nElems + elem)
! 
! 
      case default
          print *, 'ERROR: update_bc: incorrect BC type', elem, bc_type(elem)
          read *
    end select
  end do
end subroutine update_bc

end module boundc