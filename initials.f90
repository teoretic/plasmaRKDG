! initial conditions

module initials

  use names
  use servicefunc
  use solution_export
  use solvers
  use time
  use basisRKDG

implicit none

contains

  ! ============================================
  ! процедура задает начальные условия
  ! вызывается один раз в начале программы перед
  ! циклом пересчета по времени
  ! -- только для внутренних ячеек !!!

  ! см. пример на странице 128 (Chapter 7) "Part 2_3"
  subroutine init_ic
!     call blast
    call sphereOut
  end subroutine init_ic

  subroutine implosion
    integer:: elem
    real (kind=precis), dimension(1:3):: center=(/0D0,0D0,0D0/)
    real (kind=precis) :: temper, innerEn, pressure, ion, rho, bufSum, globals(3), angle,rad
    integer :: i, bn, point, info
    integer :: ipiv(4)
    
    real (kind=precis) :: pkoefM(4,4), bfRS(4), koefs(4)
    
    real (kind=precis) :: rhoM(4), rhouM(4), rhovM(4), rhowM(4), presM(4), enerM(4)


    gmm = 1.4d0 !5D0/3D0

    call allocate_solution_arrays

    print *, '*****************'
    print *, ''

    tau_out = 0.001   ! output timesteps
    t_start = 0.0    ! начальные момент времени
    t_final = 1.5    ! конечный момент времени
    t_current = t_start  ! текущий момент времени
    ntime = 1  ! Номер текущего временного слоя
    CFL = 0.3
    tau = 1D-5

    rho = 1.0D0

    do elem = 1,nElems
      do point = 1,4
        globals = loc2glob(elem,baseT(:,point))
        if (globals(1)+globals(2)> 0.15) then
          pressure = 1D0
          rho = 1D0
        else
          pressure = 0.14D0
          rho = 0.125D0          
        end if
        rhoM(point) = rho
        rhouM(point) = 0d0
        rhovM(point) = 0d0 !-rho*globals(3)
        rhowM(point) = 0d0 !rho*globals(2)
        enerM(point) = pressure/(gmm-1d0) + 0.5*(rhouM(point)*rhouM(point)+rhovM(point)*rhovM(point)+rhowM(point)*rhowM(point))/rhoM(point)
      end do
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv,  rhoM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhouM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhovM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhowM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, enerM, 4, info)
      U(1,elem,:) = rhoM
      U(2,elem,:) = rhouM
      U(3,elem,:) = rhovM
      U(4,elem,:) = rhowM
      U(5,elem,:) = enerM
    end do

    do elem = 1,nElems_b
!       bc_type(elem) = 3
!       if( abs(dot_product(surf_normal(:,bc_surf(elem)),center))>0.5) then
        bc_type(elem) = 2
!       end if
    end do

  end subroutine implosion
  
  subroutine blast
    integer:: elem
    real (kind=precis), dimension(1:3):: center=(/0.15D0,0.15D0,0D0/)
    real (kind=precis) :: temper, innerEn, pressure, ion, rho, bufSum, globals(3), angle,rad
    integer :: i, bn, point, info
    integer :: ipiv(4)
    
    real (kind=precis) :: pkoefM(4,4), bfRS(4), koefs(4)
    
    real (kind=precis) :: rhoM(4), rhouM(4), rhovM(4), rhowM(4), presM(4), enerM(4)


    gmm = 5D0/3D0

    call allocate_solution_arrays

    print *, '*****************'
    print *, ''

    tau_out = 0.001   ! output timesteps
    t_start = 0.0    ! начальные момент времени
    t_final = 1.5    ! конечный момент времени
    t_current = t_start  ! текущий момент времени
    ntime = 1  ! Номер текущего временного слоя
    CFL = 0.5
    tau = 1D-5

    rho = 1.0D0

    do elem = 1,nElems
      do point = 1,4
        globals = loc2glob(elem,baseT(:,point))
        rho = 1D0
        if (distance3(globals,center)< 0.03) then
          pressure = 10D0
        else
          pressure = 0.3D0
        end if
        rhoM(point) = rho
        rhouM(point) = 0d0
        rhovM(point) = 0d0 !-rho*globals(3)
        rhowM(point) = 0d0 !rho*globals(2)
        enerM(point) = pressure/(gmm-1d0) + 0.5*(rhouM(point)*rhouM(point)+rhovM(point)*rhovM(point)+rhowM(point)*rhowM(point))/rhoM(point)
      end do
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv,  rhoM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhouM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhovM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhowM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, enerM, 4, info)
      U(1,elem,:) = rhoM
      U(2,elem,:) = rhouM
      U(3,elem,:) = rhovM
      U(4,elem,:) = rhowM
      U(5,elem,:) = enerM
    end do

    do elem = 1,nElems_b
!       bc_type(elem) = 3
!       if( abs(dot_product(surf_normal(:,bc_surf(elem)),center))>0.5) then
        bc_type(elem) = 2
!       end if
    end do

  end subroutine blast  
  

  subroutine sphereOut
    integer:: elem
    real (kind=precis), dimension(1:3):: center=(/0D0,0D0,0D0/), ex=(/1D0,0D0,0D0/), symm=(/0d0,0d0,-1d0/)
    real (kind=precis) :: temper, innerEn, pressure, ion, rho, bufSum, globals(3), angle,rad
    integer :: i, bn, point, info
    integer :: ipiv(4)
    
    real (kind=precis) :: pkoefM(4,4), bfRS(4), koefs(4)
    
    real (kind=precis) :: rhoM(4), rhouM(4), rhovM(4), rhowM(4), presM(4), enerM(4)


    gmm = 7D0/5D0

    call allocate_solution_arrays

    print *, '*****************'
    print *, ''

    tau_out = 0.00002   ! output timesteps
    t_start = 0.0    ! начальные момент времени
    t_final = 0.02    ! конечный момент времени
    t_current = t_start  ! текущий момент времени
    ntime = 1  ! Номер текущего временного слоя
    CFL = 0.5
    tau = 1D-8

    rho = 88.6853D0
    pressure = 5529D3

    do elem = 1,nElems
      do point = 1,4
        globals = loc2glob(elem,baseT(:,point))
        rhoM(point) = rho
        rhouM(point) = 0d0
        rhovM(point) = 0d0 !-rho*globals(3)
        rhowM(point) = 0d0 !rho*globals(2)
        enerM(point) = pressure/(gmm-1d0) + 0.5*(rhouM(point)*rhouM(point)+rhovM(point)*rhovM(point)+rhowM(point)*rhowM(point))/rhoM(point)
      end do
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv,  rhoM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhouM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhovM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, rhowM, 4, info)
      pkoefM = bfkoefM/sqrt(elems_vol(elem))
      call dgesv(4, 1, pkoefM, 4, ipiv, enerM, 4, info)
      U(1,elem,:) = rhoM
      U(2,elem,:) = rhouM
      U(3,elem,:) = rhovM
      U(4,elem,:) = rhowM
      U(5,elem,:) = enerM
    end do

    do elem = 1,nElems_b
!       bc_type(elem) = 3
!       if( abs(dot_product(surf_normal(:,bc_surf(elem)),center))>0.5) then
      rad = distance3(elems_center(:,bc_elems(elem)),center)
      if (rad<0.2 .or. dot_product(surf_normal(:,bc_surf(elem)),symm)>0.5) then
        bc_type(elem) = 2
      else
        if(dot_product(surf_normal(:,bc_surf(elem)),ex)>0.5) then
          bc_type(elem) = 3
        else
          bc_type(elem) = 1
        end if
      end if
      
!       end if
    end do

  end subroutine  
  
  
  
end module initials
