
! ========================================
! Описатели сетки и работа с ними.
! Работа с геометрией.
! ========================================

! =====================================================================
! 1. Узлы КЭ обходятся против часовой стрелки
! 2. Ребра КЭ обходятся против часовой стрелки -- см. комментарии ниже
! 3. Внешняя часть расчетной области соотвествует КЭ с номером 0 (ноль) --
!    всюду, где это может потребоваться
! 4. Файлы с данными для этой программы -- в форматном (текстовом) виде --
!    это нужно для переносимости программы между системами и компиляторами
! 5. Для (r,z) геометрии ось r направлена горизонтально вправо,
!    ось z -- вертикально вверх
!    Для (x,y) направление осей обычное.
!    x = r, y = z
! ======================================================================

! =================================================================================
! Раскомментированны массивы, которые _уже_ используются в главной программе
! на данный момент.
! Для расчета восполнения могут понадобиться и другие.
! =================================================================================

module mesh3D

  use names

  use servicefunc

  implicit none

 contains

  ! ================================
  ! Инициализация всех описателей сетки
  subroutine init_mesh()

    ! used:
    !     nodes, elems, elems_elems, elems_vol, bc_elems, bc_nodes,
    !     bc_type, n, h, is_bc_elem

    implicit none

    integer :: node, elem, surf

    print *
    print *, 'init_mesh()'
    print *, '==========='

    open(101, file='./data/mesh/nodes.dat')
    open(102, file='./data/mesh/elems.dat')
    open(103, file='./data/mesh/surfaces.dat')
    open(104, file='./data/mesh/bc_elems.dat')

    read(101,*) nNodes
    print *, '  > nNodes   =', nNodes
    read(102,*) nElems
    print *, '  > nElems   =', nElems
    read(103,*) nSurfaces
    print *, '  > nSurfaces =', nSurfaces
    read(104,*) nElems_b
    print *, '  > nElems_b =', nElems_b

    print *, '== nodes =='
    allocate(nodes(1:3,1:nNodes), stat = err)
    print *, "Nodes allocate stat = ", err

    print *, '== elems =='
    allocate(elems(1:4,1:nElems), stat = err)
    print *, "elems allocate stat = ", err
    allocate(elems_elems(1:4,1:nElems), stat = err)
    print *, "elems_elems allocate stat = ", err
    allocate(elems_vol(1:nElems+nElems_b), stat = err)
    print *, "elems_vol allocate stat = ", err
    allocate(elems_center(1:3,1:nElems), stat = err)
    print *, "elems_center allocate stat = ", err
    allocate(elems_subdom(1:nElems), stat = err)
    print *, "elems_subdom allocate stat = ", err

    print *, '== surfaces =='
    allocate(surfaces(1:3,nSurfaces), stat = err)
    print *, "surfaces allocate stat = ", err
    allocate(surf_elems(1:2,nSurfaces), stat = err)
    print *, "surf_elems allocate stat = ", err
    allocate(surf_elemsLocNum(1:2,nSurfaces), stat = err)
    print *, "surf_elemsLocNum allocate stat = ", err
    allocate(elems_surfs(1:4,nElems+nElems_b), stat = err)
    print *, "elems_surfs allocate stat = ", err
    allocate(elems_surfs_num(1:4,nElems+nElems_b), stat = err)
    print *, "elems_surfs_num allocate stat = ", err
    allocate(elems_surfs_normSign(1:4,nElems+nElems_b), stat = err)
    print *, "elems_surfs_normSign allocate stat = ", err
    allocate(surf_area(nSurfaces), stat = err)
    print *, "surf_area allocate stat = ", err
    allocate(surf_normal(1:3,nSurfaces), stat = err)
    print *, "surf_normal allocate stat = ", err
    allocate(surf_center(1:3,nSurfaces), stat = err)
    print *, "surf_center allocate stat = ", err
    

    print *, '== boundary =='
    allocate(bc_type(nElems_b), stat = err)
    print *, "bc_type allocate stat = ", err
    allocate(bc_elems(nElems_b), stat = err)
    print *, "bc_elems allocate stat = ", err
    allocate(bc_nodes(1:3,nElems_b), stat = err)
    print *, "bc_nodes allocate stat = ", err
    allocate(bc_surf(nElems_b), stat = err)
    print *, "bc_surf allocate stat = ", err
    allocate(bc_elem_surf(nElems_b), stat = err)
    print *, "bc_elem_surf allocate stat = ", err
    allocate(bc_surf_center(1:3,1:nElems_b), stat = err)
    print *, "bc_surf_center allocate stat = ", err
    
    do node = 1,nNodes
      read (101,*) nodes(:,node)
    enddo

    close(101)

    do elem = 1,nElems
      read (102,*) elems(:,elem), elems_elems(:,elem), elems_vol(elem), elems_center(:,elem), &
                   elems_subdom(elem), elems_surfs(1:4,elem)
    enddo

    close(102)

    do surf = 1,nSurfaces
      read (103,*) surfaces(1:3,surf), surf_elems(1:2,surf), surf_elemsLocNum(1:2, surf), &
                   surf_normal(1:3,surf), surf_area(surf)

      surf_center(:,surf) = 0D0
      do node=1,3
        surf_center(:,surf) = surf_center(:,surf) + nodes(:, surfaces(node,surf))
      end do
      surf_center(:,surf) = surf_center(:,surf)/3D0
!      elems_surfs_normSign(surf_elemsLocNum(1, surf),surf_elems(1,surf)) = 1
!      elems_surfs_normSign(surf_elemsLocNum(2, surf),surf_elems(2,surf)) = -1
    enddo



    close(103)

    do elem = 1,nElems_b
      read (104,*) bc_nodes(:,elem), bc_elems(elem), bc_type(elem), bc_surf(elem), bc_elem_surf(elem)

      elems_vol(nElems+elem) = elems_vol(bc_elems(elem))
      
      surf_elems(2,bc_surf(elem)) = nElems+elem
      surf_elemsLocNum(2,bc_surf(elem)) = 1
      bc_surf_center(:,elem) = 0D0
      do node = 1,3
        bc_surf_center(:,elem) = bc_surf_center(:,elem) + nodes(:,bc_nodes(node,elem))/3D0
      end do
    enddo

    do surf = 1,nSurfaces
      elems_surfs_num(surf_elemsLocNum(1, surf),surf_elems(1,surf)) = 1
      elems_surfs_num(surf_elemsLocNum(2, surf),surf_elems(2,surf)) = 2
      elems_surfs_normSign(surf_elemsLocNum(1, surf),surf_elems(1,surf)) = 1
      elems_surfs_normSign(surf_elemsLocNum(2, surf),surf_elems(2,surf)) = -1
    enddo

    close(104)

    print *
    print *, 'init_mesh() DONE'
    print *, '================'

  end subroutine init_mesh


  subroutine hmin_compute
    integer :: elem, surf
    real (kind=precis) :: cand

    print *
    print *, 'hmin_compute()'
    print *, '==========='

    hmin = 1.0

    do elem = 1,nElems
      do surf = 1,4
        cand = 3.0_precis*elems_vol(elem)/surf_area(elems_surfs(surf,elem))
        if (cand < hmin) then
          hmin = cand
        endif
      end do
    end do

    print *, 'hmin = ', hmin
    print *
    print *, 'hmin_compute() DONE'
    print *, '================'
  end subroutine hmin_compute

end module mesh3D
