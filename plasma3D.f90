! основная программа

program plasma
  use names
  use mesh3D
  use solvers
  use time
  use initials
  use solvers

  use basisRKDG
  
  
  integer :: it
  
  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)
  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour

  
  call init_basis
  call init_quadrature

  ! ===============================
  ! 1. Инициализация описателей сетки и массивов для
  ! граничных условий
  call init_mesh  ! описатели для треугольных элементов
! call init_aleph  
  call hmin_compute




  ! ======================================
  ! 3. Задание начальных и граничных условий
  !    Инициализация временных параметров. Задание
  !    - начального момента времени
  !    - конечного момента времени
  !    - начального временного шага
  !    - начальное значение номера временного слоя

  call init_ic
  CFLInit = CFL
  CFLChangeStep = ntime
  CFLChanged = .false.

  call save_vtk_solution_binary !_big

  timestepping: do

    if (mod(ntime,50)==0) then
      print *, '***************'
      print *, 'Information'
      print *, 'ntime =', ntime
      print *, 't_current =', t_current
      print *, 'tau =', tau
      print *, 'max SR element'
      print *, elems_center(:,nmsr)
      
      print *, 'CFL = ', CFL
!       print *, 'radDecFails = ', radDecFails
      print *, '***************'
      print *
    end if
! 
    ntime = ntime + 1 ! номер очередного временного шага
    t_current = t_current + tau ! текущий момент времени

    
    do elem=1,nElems
      do it = 1, nBF
        U_til(:,elem, it) = U(:,elem, it)
      end do
    end do


    call gas_rkdg


    do elem=1,nElems
      do it = 1, nBF
        U(:,elem,it) = U_hat(:,elem,it)
      end do
    end do

    ! ====================================================
    ! сохраняем решение на текущем временном слое
    if (t_current-tau < tau+tau_out*floor((t_current-tau)/tau_out) .and. ntime>1) then
      print *
      print *, '***************'
      print *, 'Solution print'
      call save_vtk_solution_binary !_big
      print *, 'max SR element'
      print *, elems_center(:,nmsr)
      print *, 'tau = ', tau
      print *, 't_current =', t_current
      print *, 'CFL = ', CFL      
      print *, '***************'
      print *
    end if

    if (CFLChanged) call CFLRestore
    call time_step_compute

    ! ==========================================
    ! Проверяем условия окончания расчета по времени
    if (is_finished_time()) then
      print *, "Calculation over"
      print *, "t = ", t_current
      print *, "flag = ", is_finished_time()
      call save_vtk_solution_binary
      exit timestepping
    end if

  end do timestepping

  call DATE_AND_TIME(ddate, dtime, dzone, dvalues)
  dhour = dvalues(5)
  dmin = dvalues(6)
  dsec = dvalues(7)
  dmils = dvalues(8)

  ddtime = dmils + 1000*dsec+60000*dmin+3600000*dhour - ddtime

  print *, 'Full computation time (msec): ', ddtime

end program plasma
