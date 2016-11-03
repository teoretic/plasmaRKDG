! =========================================================
! В модуле содержится информация, относящаяяся к временной сетке
!   - шаг по времени
!   - момент времени начала расчета
!   - момент времени окончания расчета
!   - номер временного слоя
!   и т.д.
! =========================================================

module time
  use names
  use servicefunc
  use solution_export
  use basisRKDG

  implicit none

contains

  subroutine CFLRestore
    if (mod(ntime-CFLChangeStep+1,50) == 0) then
      CFL = min(CFLInit,CFL/CFLFactor)
    end if
  end subroutine

  ! ==============================================
  ! Проверка условия окончания расчета по времени
  function is_finished_time() result (tm_flag)
    implicit none

    logical:: tm_flag

    tm_flag = .FALSE.

    if (t_current > t_final) tm_flag = .TRUE.

  end function is_finished_time

  subroutine time_step_compute

      tau = CFL*hmin/SpecRadius_gas()
  end subroutine time_step_compute


  function SpecRadius_gas() result (SpecRadius1)
    real (kind=precis) :: SpecRadius1
    integer :: elem, nEdge, point
    real (kind=precis) :: rho, v1, v2, v3, E, p, k_e, asound
    real (kind=precis), dimension (1:5) :: UU
    real (kind=precis), dimension (1:3,1:3) :: Tb
    real (kind=precis), dimension (1:2) :: nT
    integer :: bn
    real (kind=precis), dimension(1:3):: cellCenterLoc=(/0.25D0,0.25D0,0.25D0/)

    SpecRadius1 = 0.0
    do elem=1,nElems
      do point = 1,4
        UU = 0D0
        do bn=1,nBF
          UU = UU + U(:,elem,bn)*basis(bn,baseT(:,point))
        end do
        UU = UU/sqrt(elems_vol(elem))
        
        
        rho = UU(1)
        v1 = UU(2)/rho
        v2 = UU(3)/rho
        v3 = UU(4)/rho
        E = UU(5)
        k_e = v1*v1+v2*v2+v3*v3 !CHECK!
        p = (E-0.5*rho*k_e)*(gmm-1.0)
  !       print *, p
        if (p*gmm/rho<0.0) then
          print *, '***************'
          print *, 'Sound speed < 0'
          print *, '***************'
          print *, 'ntime = ',ntime
          print *, 'tau = ',tau
          print *, 'nElem = ', elem
          print *, 'Elems_center = ',elems_center(:,elem)
          print *, '***************'
          print *, 'rho = ', rho
          print *, 'u = ', v1
          print *, 'v = ', v2
          print *, 'w = ', v3
          print *, 'E = ', E
          print *, trouble(elem)
          print *, getUPresLoc(elem,cellCenterLoc,U)
          print *, '***************'
          print *, 'p = ', p

          print *, '***************'
          print *, 'Solution print'
          call save_vtk_solution_binary_big !_big
          print *, '***************'

          read *
        end if
        if (SpecRadius1<sqrt(k_e)+sqrt(p*gmm/rho)) then
          SpecRadius1 = sqrt(k_e)+sqrt(p*gmm/rho)
          nmsr = elem
        end if
      end do
    end do
  end function

end module time