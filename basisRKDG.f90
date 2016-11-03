module basisRKDG

  use names
  use servicefunc

  
  contains

  subroutine init_quadrature
  
    tetraGaussN = 4
    allocate(tetraGaussW(tetraGaussN), stat = err)
    allocate(tetraGaussP(3, tetraGaussN), stat = err)
    tetraGaussW(:) = 1D0/24D0
    tetraGaussP(:,1) = (/ 5D0 + 3D0*sq5, 5D0 -     sq5, 5D0 -     sq5 /)
    tetraGaussP(:,2) = (/ 5D0 -     sq5, 5D0 + 3D0*sq5, 5D0 -     sq5 /)
    tetraGaussP(:,3) = (/ 5D0 -     sq5, 5D0 -     sq5, 5D0 + 3D0*sq5 /)
    tetraGaussP(:,4) = (/ 5D0 -     sq5, 5D0 -     sq5, 5D0 -     sq5 /)
    tetraGaussP(:,:) = tetraGaussP(:,:)/20D0
    
    trianGaussN = 3
    allocate(trianGaussW(trianGaussN), stat = err)
    allocate(trianGaussP(2, trianGaussN), stat = err)
    trianGaussW(:) = 1D0/6D0
    trianGaussP(:,1) = (/ 1D0, 1D0 /)
    trianGaussP(:,2) = (/ 4D0, 1D0 /)
    trianGaussP(:,3) = (/ 1D0, 4D0 /)
    trianGaussP(:,:) = trianGaussP(:,:)/6D0
    
  end subroutine
    
  subroutine init_basis
  
    integer :: bn, point
!     if (order == 1) then
    nBF = 4
!     end if
    baseT(:,1) = (/ 0d0, 0d0, 0d0 /)
    baseT(:,2) = (/ 1d0, 0d0, 0d0 /)
    baseT(:,3) = (/ 0d0, 1d0, 0d0 /)
    baseT(:,4) = (/ 0d0, 0d0, 1d0 /)
    
    do point=1,4
      do bn=1,nBF
        bfkoefM(point,bn) = basis(bn,baseT(:,point))
      end do
    end do

    gradBasisLoc(:,1) = 0D0
    gradBasisLoc(:,2) = (/ 0D0,    0D0,      4D0*sq5/sq3 /)
    gradBasisLoc(:,3) = (/ 0D0,    sq3*sq10, sq10/sq3 /)
    gradBasisLoc(:,4) = (/ 2*sq10, sq10,     sq10 /)
  end subroutine
  
  function basis (k, p) result (res)
    integer :: k ! номер базисной функции
    real(kind=precis), dimension(3) :: p ! локальные координаты точки
    real(kind=precis) :: res
  
    if (k==1) res = 1D0
    if (k==2) res = sq5*(4*p(3)-1)/sq3
    if (k==3) res = sq10*(p(3)+3*p(2)-1)/sq3
    if (k==4) res = sq10*(p(3)+p(2)+2*p(1)-1)
  end function
  
  function gradBasis(elem,bn) result (res)
    integer :: elem ! номер элемента
    integer :: bn ! номер базисной функции
    real(kind=precis), dimension(3) :: res
    real(kind=precis) :: trM(3,3)
    integer :: node
    logical :: stat

    res = gradBasisLoc(:,bn)/sqrt(elems_vol(elem))
    do node=1,3
      trM(:,node) = nodes(:,elems(node+1,elem)) - nodes(:,elems(1,elem))
    end do
    res = matmul(Transpose(inverseM3(trM, stat)),res)
  end function
  
  function glob2loc(elem, glob) result (loc)
    integer :: elem
    real(kind=precis) :: glob(3)
    real(kind=precis) :: loc(3)
    real(kind=precis) :: trM(3,3)
    integer :: node
    logical :: stat
    
    
    do node=1,3
      trM(:,node) = nodes(:,elems(node+1,elem)) - nodes(:,elems(1,elem))
    end do
    
    loc = matmul(inverseM3(trM, stat), glob-nodes(:,elems(1,elem)) )
    
  end function
  
  function loc2glob(elem,loc) result(glob)
    integer :: elem
    real(kind=precis) :: glob(3)
    real(kind=precis) :: loc(3)
    real(kind=precis) :: trM(3,3)
    integer :: i

    do i=1,3
      trM(:,i) = nodes(:,elems(i+1,elem)) - nodes(:,elems(1,elem))
    end do        
    
    glob = matmul(trM, loc) + nodes(:,elems(1,elem))
  end function

  function locTri2glob(surf, loc) result(glob)
    integer:: surf
    real(kind=precis) :: glob(3)
    real(kind=precis) :: loc(2)
    real(kind=precis) :: trM(3,2)    
  
    do i=1,2
      trM(:,i) = nodes(:,surfaces(i+1,surf)) - nodes(:,surfaces(1,surf))
    end do

    glob = matmul(trM, loc) + nodes(:,surfaces(1,surf))    
  end function
  
  function getUCompLoc(elem,comp,locals,arr) result (UComp)
    integer :: elem, comp
    real(kind=precis) :: locals(3), arr(:,:,:)
    real(kind=precis) :: UComp
    integer :: bn

    UComp = 0D0 
    do bn = 1, nBF
      UComp = UComp + arr(comp,elem,bn) * basis(bn,locals)
    end do      
    UComp = UComp/sqrt(elems_vol(elem))
  end function
    
  function getUPresLoc(elem,locals,arr) result (Pres)
    integer :: elem, comp
    real(kind=precis) :: locals(3), arr(:,:,:)
    real(kind=precis) :: UComp(5), Vel(3), Pres, Rho, KinE
    integer :: bn
    
    UComp = 0D0 
    do bn = 1, nBF
      UComp = UComp + arr(:,elem,bn) * basis(bn,locals)
    end do      
    UComp = UComp/sqrt(elems_vol(elem))
    Rho = UComp(1)
    Vel = UComp(2:4)/Rho
    KinE = 0.5d0*Rho*Sum(Vel**2)
    Pres = (UComp(5)-KinE)*(gmm-1.0)
  end function      
    
  function getVelLoc(elem,locals,arr) result (Vel)
    integer :: elem, comp
    real(kind=precis) :: locals(3), arr(:,:,:)
    real(kind=precis) :: UComp(5), Vel(3)
    integer :: bn
    
    UComp = 0D0 
    do bn = 1, nBF
      UComp = UComp + arr(:,elem,bn) * basis(bn,locals)
    end do      
    UComp = UComp/sqrt(elems_vol(elem))
    Vel = UComp(2:4)/UComp(1)
  end function  
  
  function getULoc(elem,locals,arr) result (UComp)
    integer :: elem, comp
    real(kind=precis) :: locals(3), arr(:,:,:)
    real(kind=precis) :: UComp(5)
    integer :: bn
    
    UComp = 0D0 
    do bn = 1, nBF
      UComp = UComp + arr(:,elem,bn) * basis(bn,locals)
    end do      
    UComp = UComp/sqrt(elems_vol(elem))
  end function  
  
  function getUComp(elem,comp,globals,arr) result (UComp)
    integer :: elem, comp
    real(kind=precis) :: globals(3), locals(3), arr(:,:,:)
    real(kind=precis) :: UComp
    integer :: bn
    
    locals = glob2loc(elem,globals)
    UComp = 0D0 
    do bn = 1, nBF
      UComp = UComp + arr(comp,elem,bn) * basis(bn,locals)
    end do      
    UComp = UComp/sqrt(elems_vol(elem))
  end function

  function getU(elem,globals,arr) result (UComp)
    integer :: elem, comp
    real(kind=precis) :: globals(3), locals(3), arr(:,:,:)
    real(kind=precis) :: UComp(5)
    integer :: bn
    
    locals = glob2loc(elem,globals)
    UComp = 0D0 
    do bn = 1, nBF
      UComp = UComp + arr(:,elem,bn) * basis(bn,locals)
    end do      
    UComp = UComp/sqrt(elems_vol(elem))
  end function

  function getVel(elem,globals,arr) result (Vel)
    integer :: elem, comp
    real(kind=precis) :: globals(3), locals(3), arr(:,:,:)
    real(kind=precis) :: UComp(5), Vel(3)
    integer :: bn
    
    locals = glob2loc(elem,globals)
    UComp = 0D0 
    do bn = 1, nBF
      UComp = UComp + arr(:,elem,bn) * basis(bn,locals)
    end do      
    UComp = UComp/sqrt(elems_vol(elem))
    Vel = UComp(2:4)/UComp(1)
  end function

  function getMedianOnCellComp(elem,comp,arr) result (UComp)
    integer :: elem, comp
    real(kind=precis) :: arr(:,:,:)
    real(kind=precis) :: UComp
    real (kind=precis), dimension(1:3) :: local_center = (/0.25,0.25, 0.25/)

    UComp = arr(comp,elem,1)*basis(1,local_center)/sqrt(elems_vol(elem))
  end function

  function getGradComp(elem,comp,arr) result (grad)
    integer :: elem, comp
    real(kind=precis) :: arr(:,:,:)
    real(kind=precis) :: grad(3)
    integer :: bn
    
    grad = 0D0
    
    do bn=1,nBF
      grad = grad + arr(comp,elem,bn)*gradBasis(elem,bn)
    end do
    
  end function  
  
end module
