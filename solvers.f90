! решатели

module solvers

  use names
  use servicefunc
  use boundc
  use flux
  use basisRKDG

  implicit none

contains

subroutine troubleIndicator(troubArr, solArr)
  real (kind=precis), Intent(In) :: solArr(:,:,:)
  logical, Intent(InOut) :: troubArr(:)
  
  integer :: nInd, indi, surf, elem, node, nNegat, nPosit, surfGlob, comp
  integer :: indikator(8), negat(3), posit(3)
  real (kind=precis), dimension(1:3) :: p,q,r,s,d,g,vel,globals,normal
  logical :: flag_trouble
  real (kind=precis) :: A, uSc, negArea, krivoInd, uu(5)
  real (kind=precis), dimension(1:2) :: intQ, normQ  
  real (kind=precis), dimension(1:3) :: b
  real (kind=precis), dimension(1:3):: cellCenterLoc=(/0.25D0,0.25D0,0.25D0/)
    
  nInd = 2
  indikator(1) = 1
  indikator(2) = 5
  ! krivodonova trouble cell indicator
  
  !$OMP PARALLEL DEFAULT(none) Shared(troubArr,trInd,nInd,indikator,solArr,elems_elems,nElems,surf_normal,surfaces,limitThreshold,elems_vol,surf_area,surf_center,nodes,elems_surfs) PRIVATE(flag_trouble,intQ,normQ,negArea,indi,elem,surf,surfGlob,normal,nNegat,nPosit,node,globals,vel,b,negat,posit,p,q,r,s,d,g,krivoInd,A)
  !$OMP DO SCHEDULE(auto)
  do elem = 1, nElems
    troubArr(elem) = .false.
    flag_trouble = .false.

    trInd(elem) = 0D0
    intQ = 0D0
    normQ = 0D0
    negArea = 0D0
    do indi = 1, nInd
      normQ(indi) = getMedianOnCellComp(elem,indikator(indi),solArr)
    end do

    do surf = 1,4
      if (elems_elems(surf, elem)<=nElems) then
        surfGlob = elems_surfs(surf,elem)
        normal = surf_normal(:,surfGlob)
        nNegat = 0
        nPosit = 0
        do node = 1,3
          globals = nodes(:, surfaces(node,surfGlob))
          vel = getVel(elem,globals,solArr)
          b(node) = dot_product(vel,normal)
          if (b(node)<0D0) then
            nNegat = nNegat + 1
            negat(nNegat) = node
          else 
            nPosit = nPosit + 1
            posit(nPosit) = node
          end if
        end do !node
            
        if (nNegat == 3) then
            globals = surf_center(:, surfGlob)
            negArea = surf_area(surfGlob)
        end if
        if (nNegat == 2) then
            r = nodes(:, surfaces(negat(1),surfGlob))
            q = nodes(:, surfaces(negat(2),surfGlob))
            d = nodes(:, surfaces(posit(1),surfGlob))
            s = d + (r-d)*b(posit(1))/(b(posit(1))-b(negat(1)))
            p = d + (q-d)*b(posit(1))/(b(posit(1))-b(negat(2)))
            globals = (p+q+r+s)*0.25D0
            negArea = 0.5D0*norm2_3(vect_mult(s-q,p-r))
        end if
        if (nNegat == 1) then
            p = nodes(:, surfaces(negat(1),surfGlob))
            d = nodes(:, surfaces(posit(1),surfGlob))
            g = nodes(:, surfaces(posit(2),surfGlob))
            q = d + (p-d)*b(posit(1))/(b(posit(1))-b(negat(1)))
            r = d + (g-d)*b(posit(2))/(b(posit(2))-b(negat(1)))
            globals = (p+q+r)/3.0D0
            negArea = 0.5D0*norm2_3(vect_mult(p-q,r-q))
        end if
            
        if (nNegat>0) then
          do indi = 1, nInd
            intQ(indi) = negArea * abs(getUComp(elem,indikator(indi),globals,solArr)-getUComp(elems_elems(surf,elem),indikator(indi),globals,solArr))
            krivoInd = 0D0
            if (negArea>1D-8) then
              A = (elems_vol(elem)**(1D0/3D0))*1.24902
              krivoInd = intQ(indi)/(normQ(indi)*negArea*A)
            end if
            if (krivoInd > limitThreshold) then
              flag_trouble = .true.
            end if
            trInd(elem) = max(krivoInd,trInd(elem))
          end do !indi
        end if
      end if      
    end do !surf
    troubArr(elem) = flag_trouble
    if (flag_trouble) then
      do surf=1,4
        if (elems_elems(surf,elem)<=nElems) then
          troubArr(elems_elems(surf,elem)) = flag_trouble
        end if
      end do
    end if
  end do !elem  

  !$OMP END DO nowait
  !$omp end parallel  
  
!   do elem = 1, nElems
!     if (troubArr(elem)) then
!       do surf=1,4
!         if (elems_elems(surf,elem)<=nElems) then
!           troubArr(elems_elems(surf,elem)) = troubArr(elem)
!         end if
!       end do
!     end if
!   end do
end subroutine

subroutine limiterHWENO(troubArr, solArr)
  real (kind=precis), Intent(InOut) :: solArr(:,:,:)
  logical, Intent(In) :: troubArr(:)
  real (kind=precis), dimension(9,3) :: cand_coefs
  real (kind=precis) :: medians(4), cooMatr(3,3), grads(3), RH(3)
  real (kind=precis) :: limGmm, limWei(9), limSumWei, limEps, limGrad(3)
  integer :: elem, comp, bn, surf, node
  logical :: flag_trouble  
  integer :: nCand, neib, sten,cand
  real (kind=precis), dimension(1:3) ::  randVec  
  integer :: ii
  logical :: stat
  logical :: presErr
  real (kind=precis), allocatable :: limSol(:,:,:)
  
  integer :: nWall, wall
  integer :: wallNormals(3,3)
  
  limGmm = 1.5
  limEps = 1D-3  
 

  allocate(limSol(nComponents,nElems,nBF))
  limSol = solArr
 
  !hweno limiter --- luo 2007
  
  !$OMP PARALLEL DEFAULT(none) Shared(limSol,troubArr,solArr,elems_elems,nElems,surf_normal,elems_surfs,limGmm,limEps,elems_center,hmin,locNeibSurf,nBF,nComponents,bc_type,baseT) PRIVATE(flag_trouble,elem,surf,node,nWall,wall,wallNormals,limWei,limSumWei,limGrad,stat,cooMatr,sten,randVec,RH,neib,medians,cand_coefs,nCand,presErr)
  !$OMP DO SCHEDULE(auto)
    
  do elem = 1, nElems 
    presErr = .false.
    do node = 1,4
      if(getUPresLoc(elem,baseT(:,node),solArr)<0d0 .or. getUComp(elem,1,baseT(:,node),solArr)<0d0) then
        presErr = .true.
        exit
      end if
    end do
    
    if (troubArr(elem) .or. presErr) then
      do node = 1,4 
        nWall = 0
        if (elems_elems(node,elem)>nElems) then
          if (bc_type(elems_elems(node,elem)-nElems) == 2) then
            nWall = nWall + 1
            wallNormals(:,nWall) = surf_normal(:,elems_surfs(node,elem))
          end if
        end if
      end do
      
      if (presErr) then
        do comp=1,nComponents
          limSol(comp,elem,2:nBF) = 0d0
        end do
      else
        do comp=1,nComponents
          nCand = 1
          ! наклон собственный
          cand_coefs(nCand,:) = getGradComp(elem,comp,solArr)
          
          medians(1) = getMedianOnCellComp(elem,comp,solArr)
          
          do surf=1,4
            ! наклоны соседей
            if (elems_elems(surf,elem)<=nElems) then
              nCand = nCand + 1
              cand_coefs(nCand,:) = getGradComp(elems_elems(surf,elem),comp,solArr)
            end if
            
            ! наклоны линейные по шаблонам
            flag_trouble = .false.
            do sten=1,3
              neib = elems_elems(locNeibSurf(sten,surf), elem)
              if (neib<=nElems) then
                medians(sten+1) = getMedianOnCellComp(neib,comp,solArr)
                RH(sten) = medians(sten+1) - medians(1)
                call RANDOM_NUMBER(randVec)
                randVec = 0d0 !randVec*hmin*0.05
                cooMatr(sten,:) = elems_center(:,neib) - (elems_center(:,elem)+randVec)
  !               print *, cooMatr(sten,:)
  !               print *, neib, elems_center(:,neib) 
  !               print *, elems_center(:,elem)
              else
                flag_trouble = .true.
                exit
              end if
            end do ! sten
            if (.not. flag_trouble) then
              nCand = nCand + 1
  !             do ii=1,3
  !               print *, cooMatr(ii,:)
  !             end do

              cand_coefs(nCand,:) = matmul(inverseM3(cooMatr, stat), RH)
              if (.not. stat) then
                nCand = nCand - 1
              end if
            end if
          end do ! surf
          
  !         print *, "00-----------------------------"
          limGrad = 0D0
          limSumWei = 0D0
          do cand=1,nCand
            limWei(cand) = norm2_3(cand_coefs(cand,:))
            limWei(cand) = (limEps+limWei(cand))**(-limGmm)
            limSumWei = limSumWei + limWei(cand)
            limGrad = cand_coefs(cand,:)*limWei(cand)
          end do ! cand
          limGrad = limGrad/limSumWei
          
          
          if (nWall>0 .and. (comp == 1 .or. comp == 5)) then
            do wall = 1, nWall
              limGrad = limGrad - wallNormals(:,wall)*dot_product(limGrad,wallNormals(:,wall))
            end do
          end if
          
  !         print *, limWei(1:nCand)/limSumWei
  !         
  !         
  !         print *, nCand
  !         do cand=1,nCand
  !           print *, cand_coefs(cand, :)
  !         end do
  !         print *, "----"
  !         print *, limGrad
          
          do bn=2,nBF
            cooMatr(:,bn-1) = gradBasis(elem,bn)
          end do ! bn
  !         U_hat(comp,elem,1) = solArr(comp,elem,1)
          limSol(comp,elem,2:nBF) = matmul(inverseM3(cooMatr, stat),limGrad)
          
  !         print *, norm2_3(U_hat(comp,elem,2:nBF)), norm2_3(solArr(comp,elem,2:nBF))
          
  !         if (norm2_3(U_hat(comp,elem,2:nBF))>norm2_3(solArr(comp,elem,2:nBF))) then
  !           print *, "aaaa"
  !         end if
          
  !         U_hat(comp,elem,2:nBF) = 0D0
  !         print *, solArr(comp,elem,2:nBF)
  !         print *, U_hat(comp,elem,2:nBF)
  !         
  !         print *, "pppppp"
  !         print *, getGradComp(elem,comp,U_hat), limGrad
  !         read *

          
  !         U_hat(comp,elem,2:nBF) = solArr(comp,elem,2:nBF)
        end do ! comp
      end if ! presErr
    end if ! troubleArr
  end do ! elem

  !$OMP END DO nowait
  !$omp end parallel    


  do elem = 1,nElems
    do comp=1,nComponents
      solArr(comp,elem,2:nBF)=limSol(comp,elem,2:nBF)	
    end do
  end do  
  deallocate (limSol)
!   do elem = 1, nElems
!     if (troubArr(elem)) then
!       do comp=1,nComponents
!         solArr(comp,elem,2:nBF) = U_hat(comp,elem,2:nBF)
!       end do ! comp
!     end if
!   end do ! elem    
end subroutine


subroutine RKStep(solSt, solRH, solEnd, solSurf, step)
  real (kind=precis), Intent(In) :: solSt(:,:,:), solRH(:,:,:), step
  real (kind=precis), Intent(InOut) :: solEnd(:,:,:), solSurf(:,:,:,:)
  
  integer :: elem, surf, point, bn
  real (kind=precis), dimension(1:5) :: Ul, Ur, gradInt
  real (kind=precis) :: cosz, cosy, sinz, siny, xyproek, A  
  real (kind=precis), dimension(1:5) :: flow_surf, flow_elem, uu
  real (kind=precis), dimension(1:3,1:3) :: T, Ti
  real (kind=precis), dimension(1:3) :: normal, globals, locals
  
  !$OMP PARALLEL DEFAULT(none) Shared(nElems,trianGaussN,elems_surfs,trianGaussP,solSurf,elems_surfs_num,solRH) PRIVATE(elem,surf,point,globals)
  !$OMP DO SCHEDULE(auto)  
  
  do elem = 1, nElems
    do surf = 1,4
      do point = 1,trianGaussN
        globals = locTri2glob(elems_surfs(surf,elem), trianGaussP(:,point))
        solSurf(:,point,elems_surfs_num(surf,elem),elems_surfs(surf,elem)) = getU(elem,globals,solRH)
      end do
    end do
  end do
  
  !$OMP END DO nowait
  !$omp end parallel      
  
  call update_bc(solSurf)

  
  !$OMP PARALLEL DEFAULT(none) Shared(nSurfaces,surf_normal,trianGaussN,solSurf,gas_surf_flux) PRIVATE(surf,normal,xyproek,cosz,sinz,cosy,siny,T,Ti,point,Ul,Ur,flow_surf)
  !$OMP DO SCHEDULE(auto)    
  
  do surf=1, nSurfaces
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
    
    do point=1,trianGaussN
      Ul(:) = solSurf(:,point,1,surf)
      Ur(:) = solSurf(:,point,2,surf)

      Ul(2:4) = matmul(Ti,Ul(2:4))
      Ur(2:4) = matmul(Ti,Ur(2:4))
      
      ! ==================================
      ! вычисляем поток из распадника
      flow_surf=flux_num_hllc(Ul(1:5),Ur(1:5))

      ! обратный поворот
      flow_surf(2:4) = matmul(T,flow_surf(2:4))

      gas_surf_flux(:,point,surf)= flow_surf
    end do
  end do
  !$OMP END DO nowait
  !$omp end parallel   
  
  !$OMP PARALLEL DEFAULT(none) Shared(nElems,nBF,elems_surfs_normSign,surf_area,elems_surfs,elems_vol,trianGaussN,gas_surf_flux,trianGaussP,trianGaussW,tetraGaussN,tetraGaussP,tetraGaussW,solRH,solEnd,solSt,step) PRIVATE(elem,bn,flow_elem,surf,A,point,globals,locals,gradInt,uu)
  !$OMP DO SCHEDULE(auto)      
    
  do elem = 1,nElems
    do bn=1,nBF
      flow_elem = 0D0
      do surf = 1,4
        A = elems_surfs_normSign(surf,elem)* 2D0*surf_area(elems_surfs(surf,elem)) /sqrt(elems_vol(elem))
        do point=1,trianGaussN
        ! ================================================================
        ! добавляем поток через данную грань в полный поток через границу ячейки
          globals = locTri2glob(elems_surfs(surf,elem), trianGaussP(:,point))
          locals = glob2loc(elem,globals)
          flow_elem = flow_elem + A*gas_surf_flux(:,point,elems_surfs(surf,elem))*basis(bn,locals)*trianGaussW(point)
        end do
      end do
      
      gradInt = 0D0
      do point = 1, tetraGaussN
        globals = loc2glob(elem,tetraGaussP(:,point))
        uu = getU(elem,globals,solRH)
        gradInt = gradInt + tetraGaussW(point)*matmul(FL_matrix(uu), gradBasis(elem,bn))
      end do
      gradInt = gradInt* 6D0*elems_vol(elem)
      ! ====================================================
      ! вычисляем новое значение в ячейке

      solEnd(:,elem,bn) = solSt(:,elem,bn) + step*(gradInt-flow_elem)
    end do
  end do  
  !$OMP END DO nowait
  !$omp end parallel        

end subroutine

subroutine RKStepBreak(solNew,solOld, SNC)
  real (kind=precis), Intent(InOut) :: solNew(:,:,:), solOld(:,:,:)
  integer :: elem, it, node
  logical, Intent(InOut) :: SNC
    
  SNC = .false.
  
  do elem=1,nElems
    do node=1,4
      if (getUPresLoc(elem,BaseT(:,node),solNew)<0d0) then
        SNC = .true.
        exit
      end if
    end do
  end do    
    
  if (SNC) then
    ntime = ntime - 1
    t_current = t_current - tau
    do elem=1,nElems
      do it = 1, nBF
        solNew(:,elem,it) = solOld(:,elem,it)
      end do
    end do
    CFL = CFLFactor*CFL
    CFLChangeStep = ntime
    CFLChanged = .true.    
  end if
end subroutine

subroutine gas_rkdg
  integer :: elem
  real (kind=precis), dimension(1:3):: cellCenterLoc=(/0.25D0,0.25D0,0.25D0/)
  logical :: StepNotCorrect
  
  call RKStep(U_til,U_til,U_hat,U_surf,0.5*tau)
  call troubleIndicator(trouble, U_hat)
  call limiterHWENO(trouble, U_hat)
  call RKStepBreak(U_hat, U, StepNotCorrect) 
  if (StepNotCorrect) RETURN
  
  call RKStep(U_til,U_hat,U_hat,U_surf,tau)
  call troubleIndicator(trouble, U_hat)
  call limiterHWENO(trouble, U_hat)
  call RKStepBreak(U_hat, U, StepNotCorrect) 
  if (StepNotCorrect) RETURN
  
  
  
!   !$OMP PARALLEL DEFAULT(none) Shared(U_til,nSurfaces,gas_surf_flux,surf_normal,surf_area,surf_elems) PRIVATE(surf,elem1,elem2,flow_surf,Ul,Ur,normal,cosz,cosy,sinz,siny,xyproek,T,Ti)
!   !$OMP DO SCHEDULE(auto)
!   do surf = 1, nSurfaces
!     ! вычисление потока через поверхность
!     ! цикл по поверхностям
! !     elem1=surf_elems(1,surf)
! !     elem2=surf_elems(2,surf)
! ! 
! !     normal(:) = surf_normal(:,surf)
! !     Ul(:) = U_til(:,elem1)
! !     Ur(:) = U_til(:,elem2)
! ! 
! !     xyproek = sqrt(Sum(normal(1:2)**2))
! !     if (xyproek > 1e-7) then
! !       cosz = normal(1)/xyproek
! !       sinz = normal(2)/xyproek
! !     else
! !       cosz = 1D0
! !       sinz = 0D0
! !     endif
! ! 
! !     cosy = xyproek
! !     siny = -normal(3)
! ! 
! !     !поворот
! ! 
! !     T(1,1:3) = (/ cosy*cosz,  -sinz, cosz*siny /)
! !     T(2,1:3) = (/ sinz*cosy,   cosz, sinz*siny /)
! !     T(3,1:3) = (/     -siny,    0D0,      cosy /)
! ! 
! !     Ti = Transpose(T)
! ! 
! !     Ul(2:4) = matmul(Ti,Ul(2:4))
! !     Ur(2:4) = matmul(Ti,Ur(2:4))
! !     ! ==================================
! !     ! вычисляем поток из распадника
! !     flow_surf=flux_num_hllc(Ul(1:5),Ur(1:5))
! ! 
! ! 
! !     ! обратный поворот
! ! 
! !     flow_surf(2:4) = matmul(T,flow_surf(2:4))
! ! 
! !     gas_surf_flux(surf,:)= flow_surf*surf_area(surf)
! 
!   end do
!   !$OMP END DO nowait
!   !$omp end parallel

!         print *, "<---->"

!   !$OMP PARALLEL DEFAULT(none) Shared(U_hat,elems_vol,U_til,elems_surfs_normSign,elems_surfs,gas_surf_flux,elem,tau,nElems) PRIVATE(surf,flow_elem,A)
!   !$OMP DO SCHEDULE(auto)
!   do elem = 1,nElems
!     flow_elem = 0D0
!     do surf = 1,4
!       ! ================================================================
!       ! добавляем поток через данную грань в полный поток через границу ячейки
!       flow_elem = flow_elem + elems_surfs_normSign(surf,elem)*gas_surf_flux(elems_surfs(surf,elem),:)
!     end do
!     ! ====================================================
!     ! вычисляем новое значение в ячейке
!     ! A -- коэфиициент, зависит от шага сетки и площади ячейки
!     A = -tau/elems_vol(elem)
! 
!     U_hat(1:5,elem) = U_til(1:5,elem) + flow_elem(:)*A
!   end do
!   !$OMP END DO nowait
!   !$omp end parallel

end subroutine gas_rkdg


subroutine allocate_solution_arrays

  nComponents = 5

  allocate(U(5,nElems+nElems_b,nBF), stat = err)
  print *, 'U allocate stat=', err
  allocate(U_hat(5,nElems+nElems_b,nBF), stat = err)
  print *, 'U_hat allocate stat=', err
  allocate(U_til(5,nElems+nElems_b,nBF), stat = err)
  print *, 'U_til allocate stat=', err
  
  allocate(U_surf(5,trianGaussN,2,nSurfaces), stat = err)
  print *, 'U_surf allocate stat=', err
  
  allocate(gas_surf_flux(5,trianGaussN,nSurfaces), stat = err)
  print *, 'gas_surf_flux allocate stat=', err
  
  allocate (trouble(nElems), stat = err)
  print *, 'trouble allocate stat=', err  
  allocate (trInd(nElems), stat = err)
  print *, 'trInd allocate stat=', err  
  
  
end subroutine allocate_solution_arrays

end module solvers
