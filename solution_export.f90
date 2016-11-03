! ==========================================================================
! Модуль вывода данных
! ==========================================================================


module solution_export
  use names
  use servicefunc
  use basisRKDG

  implicit none

contains


  ! функция возвращает имя файла с решением на временном слое t-1 в формате VTK
  SUBROUTINE filename_se_sol_vtk(time, name)
    IMPLICIT NONE

    INTEGER, INTENT(in):: time
    CHARACTER(LEN=*), INTENT(out):: name

    name = ''
    name = root // '/' // fname_sol // &
         delim // 'T' // i2c(time) // delim // 'vtk' // extentionVTK

  END SUBROUTINE filename_se_sol_vtk

function FixTheLine_CharArray(s) result (sout)
  implicit none
    character s*(*)
    character stemp*(Len(Trim(ADJUSTL(s)))+2)
    character,allocatable:: sout(:)
    integer :: i
    Character LF,CR
    allocate(sout(1:(Len(Trim(ADJUSTL(s)))+2)))

        LF=Char(10)
        CR=Char(13)
    stemp=Trim(ADJUSTL(s))
    stemp(Len(Trim(ADJUSTL(s)))+1:Len(Trim(ADJUSTL(s)))+1)=CR
    stemp(Len(Trim(ADJUSTL(s)))+2:Len(Trim(ADJUSTL(s)))+2)=LF
    do i=1,Len(Trim(ADJUSTL(s)))+2
    sout(i)=stemp(i:i)
    enddo
end function
! 
! 
! 
! subroutine save_vtk_solution
!   implicit none
! 
!   integer :: node, elem
!   real (kind=precis):: rho, ux, v, w, Ek, p, k_e, sounds, Mach
! 
!   call filename_se_sol_vtk(ntime,current_sol_name)
! 
!   open (unit=1,file=current_sol_name)
! 
!   write(1,"(A)") '# vtk DataFile Version 3.0'
!   write(1,"(A,f6.3)") 'Plasma3D output, t= ',t_current
!   write(1,"(A)") 'ASCII'
!   write(1,*) ''
!   write (1,"(A)") 'DATASET UNSTRUCTURED_GRID'
!   write (1,"(A,I7,A)") 'POINTS ', nNodes, ' float'
! 
!   do node=1,nNodes
!     write (1,*) nodes(:,node)
!   end do
! 
!   write (1,"(A,I7,A,I7)") 'CELLS ', nElems, ' ', nElems+4*nElems
! 
!   do elem=1,nElems
!     write (1,*) 4, elems(1,elem)-1, elems(2,elem)-1, elems(3,elem)-1, elems(4,elem)-1
!   end do
! 
!   write (1,"(A,I7)") 'CELL_TYPES ', nElems
!   do elem=1,nElems
!     write (1,*) 10
!   end do
! 
!   write (1,"(A,I7)") 'CELL_DATA ', nElems
! 
!   write (1,"(A)") 'SCALARS density float'
!   write (1,"(A)") 'LOOKUP_TABLE default'
!   do elem = 1,nElems
!     write(1,*)  U(1,elem)
!   end do
! 
!   write (1,"(A)") 'VECTORS velocity float'
!   do elem=1,nElems
!     write(1,*) U(2,elem)/U(1,elem), U(3,elem)/U(1,elem), U(4,elem)/U(1,elem)
!   end do
! 
!   write (1,"(A)") 'SCALARS energy float'
!   write (1,"(A)") 'LOOKUP_TABLE default'
!   do elem=1,nElems
!     write(1,*)  U(5,elem)
!   end do
! 
!   write (1,"(A)") 'SCALARS pressure float'
!   write (1,"(A)") 'LOOKUP_TABLE default'
!   do elem=1,nElems
!     rho = U(1,elem)
!     ux = U(2,elem)/rho
!     v = U(3,elem)/rho
!     w = U(4,elem)/rho
!     Ek = U(5,elem)
!     k_e = ux*ux+v*v+w*w
!     p = (Ek-0.5*rho*k_e)*(gmm-1.0)
! 
!     write(1,*) p
!   end do
! 
!   close (unit=1,status='Keep')
!   write (*,*) 'File ',current_sol_name,' has been created'
!   write (*,*) 't=',t_current
!   write (*,*) '************************'
! end subroutine save_vtk_solution
! 
! 
! 
! 
! 
subroutine save_vtk_solution_binary
  use ifport

  implicit none

  integer :: node, elem, bn
  real (kind=precis):: rho, ux, v, w, Ek, p, k_e, sounds, Mach, buf, UU(5)
  character out_string*200
  real (kind=precis), dimension(1:3):: cellCenterLoc=(/0.25D0,0.25D0,0.25D0/), buf_v

  call filename_se_sol_vtk(ntime,current_sol_name)

  open (unit=1,file=current_sol_name,form='BINARY')

  write (out_string,"(A)") '# vtk DataFile Version 3.0'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A,f10.7)") 'Plasma3D output, t= ',t_current
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'BINARY'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A,I7,A)") 'POINTS ', nNodes, ' double'
  write (1) FixTheLine_CharArray(out_string)

  do node=1,nNodes
    write (1) nodes(:,node)
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A,I7,A,I7)") 'CELLS ', nElems, ' ', nElems+4*nElems
  write (1) FixTheLine_CharArray(out_string)

  do elem=1,nElems
    write (1) 4, elems(1,elem)-1, elems(2,elem)-1, elems(3,elem)-1, elems(4,elem)-1
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A,I7)") 'CELL_TYPES ', nElems
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    write (1) 10
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A,I7)") 'CELL_DATA ', nElems
  write (1) FixTheLine_CharArray(out_string)

  write (out_string,"(A)") 'SCALARS density double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem = 1,nElems
    write (1) getUCompLoc(elem,1,cellCenterLoc,U) 
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A)") 'VECTORS velocity double'
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    write (1) getVelLoc(elem,cellCenterLoc,U)
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A)") 'SCALARS energy double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    write (1)  getUCompLoc(elem,5,cellCenterLoc,U) 
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A)") 'SCALARS pressure double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    
    UU = getULoc(elem,cellCenterLoc,U)
    rho = UU(1)
    Ek = UU(5)
    ux = UU(2)/rho
    v  = UU(3)/rho
    w  = UU(4)/rho
    k_e = ux*ux+v*v+w*w
    p = (Ek-0.5*rho*k_e)*(gmm-1.0)

    write (1) p
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A)") 'SCALARS trInd double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem = 1,nElems
      write (1) trInd(elem)
  end do
  write (1) FixTheLine_CharArray(' ')  
  
  write (out_string,"(A)") 'SCALARS KXRCF double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem = 1,nElems
    if (trouble(elem)) then
      write (1) 1D0
    else
      write (1) 0D0
    end if
  end do
  write (1) FixTheLine_CharArray(' ')      
  
  
  close (unit=1,status='Keep')
  write (*,*) 'File ',current_sol_name,' has been created'
  write (*,*) 't=',t_current
  write (*,*) '************************'
end subroutine save_vtk_solution_binary




subroutine save_vtk_solution_binary_big
  use ifport

  implicit none

  integer :: node, elem, bn
  real (kind=precis):: rho, ux, v, w, Ek, p, k_e, sounds, Mach, buf
  character out_string*200
  real (kind=precis), dimension(1:3):: cellCenterLoc=(/0.25D0,0.25D0,0.25D0/), buf_v
  real (kind=precis) :: vertices(3,4)

  vertices(:,1) = (/ 0D0, 0D0, 0D0 /)
  vertices(:,2) = (/ 1D0, 0D0, 0D0 /)
  vertices(:,3) = (/ 0D0, 1D0, 0D0 /)
  vertices(:,4) = (/ 0D0, 0D0, 1D0 /)
  
  call filename_se_sol_vtk(ntime,current_sol_name)

  open (unit=1,file=current_sol_name,form='BINARY')

  write (out_string,"(A)") '# vtk DataFile Version 3.0'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A,f6.3)") 'Plasma3D output, t= ',t_current
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'BINARY'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'DATASET UNSTRUCTURED_GRID'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A,I7,A)") 'POINTS ', nElems*4, ' double'
  write (1) FixTheLine_CharArray(out_string)

  do elem=1,nElems
    do node=1,4
      write (1) nodes(:,elems(node,elem))
    end do
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A,I7,A,I7)") 'CELLS ', nElems, ' ', nElems+4*nElems
  write (1) FixTheLine_CharArray(out_string)

  node = 0
  do elem=1,nElems
    write (1) 4, node, node+1, node+2, node+3
    node = node + 4
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A,I7)") 'CELL_TYPES ', nElems
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    write (1) 10
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A,I7)") 'POINT_DATA ', nElems*4
  write (1) FixTheLine_CharArray(out_string)

  write (out_string,"(A)") 'SCALARS density double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem = 1,nElems
    do node = 1,4
      rho = 0D0
      do bn=1,nBF
        rho = rho + U(1,elem,bn) * basis(bn,vertices(:,node))/sqrt(elems_vol(elem))
      end do
      write (1) rho
    end do
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A)") 'VECTORS velocity double'
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    do node = 1,4
      buf_v = 0D0
      rho = 0D0
      do bn=1,nBF
        buf_v = buf_v + U(2:4,elem,bn) * basis(bn,vertices(:,node))/sqrt(elems_vol(elem))
        rho = rho + U(1,elem,bn) * basis(bn,vertices(:,node))/sqrt(elems_vol(elem))
      end do
      buf_v = buf_v/rho
      write (1) buf_v
    end do
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A)") 'SCALARS energy double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    do node = 1,4    
      buf = 0D0
      do bn=1,nBF
        buf = buf + U(5,elem,bn) * basis(bn,vertices(:,node))/sqrt(elems_vol(elem))
      end do
      write (1)  buf
    end do
  end do
  write (1) FixTheLine_CharArray(' ')

  write (out_string,"(A)") 'SCALARS pressure double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem=1,nElems
    do node = 1,4
      rho = 0D0
      buf_v = 0D0  
      Ek = 0D0
      do bn=1,nBF
        rho = rho + U(1,elem,bn) * basis(bn,vertices(:,node))/sqrt(elems_vol(elem))
        buf_v = buf_v + U(2:4,elem,bn) * basis(bn,vertices(:,node))/sqrt(elems_vol(elem))
        Ek = Ek + U(5,elem,bn) * basis(bn,vertices(:,node))/sqrt(elems_vol(elem))
      end do
!       buf_v = buf_v/rho
        
      ux = buf_v(1)/rho
      v  = buf_v(2)/rho
      w  = buf_v(3)/rho
      k_e = ux*ux+v*v+w*w
      p = (Ek-0.5*rho*k_e)*(gmm-1.0)

      write (1) p
    end do
  end do
  write (1) FixTheLine_CharArray(' ')

  
  write (out_string,"(A,I7)") 'CELL_DATA ', nElems
  write (1) FixTheLine_CharArray(out_string)

  write (out_string,"(A)") 'SCALARS trInd double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem = 1,nElems
      write (1) trInd(elem)
  end do
  write (1) FixTheLine_CharArray(' ')  
  
  write (out_string,"(A)") 'SCALARS KXRCF double'
  write (1) FixTheLine_CharArray(out_string)
  write (out_string,"(A)") 'LOOKUP_TABLE default'
  write (1) FixTheLine_CharArray(out_string)
  do elem = 1,nElems
    if (trouble(elem)) then
      write (1) 1D0
    else
      write (1) 0D0
    end if
  end do
  write (1) FixTheLine_CharArray(' ')    
  
  close (unit=1,status='Keep')
  write (*,*) 'File ',current_sol_name,' has been created'
  write (*,*) 't=',t_current
  write (*,*) '************************'
end subroutine save_vtk_solution_binary_big


end module solution_export
