!==============================================================
! имена!
!==============================================================


module names

  implicit none
  
  ! code parameters
  integer, parameter :: precis = 8
  integer, parameter :: max_filename_length = 50

  
  ! method parameters
  
  integer :: order = 1
  integer :: nBF
  real (kind=precis) :: bfkoefM(4,4)
  real (kind=precis) :: baseT(3,4)
  real (kind=precis) :: gradBasisLoc(3,4)
  
  
  real (kind=precis), parameter :: sq2 = 1.4142135623730950488 !sqrt(2D0)
  real (kind=precis), parameter :: sq3 = 1.7320508075688772935 !sqrt(3D0)      
  real (kind=precis), parameter :: sq5 = 2.2360679774997896964 !sqrt(5D0)      
  real (kind=precis), parameter :: sq6 = 2.4494897427831780982 !sqrt(6D0)
  real (kind=precis), parameter :: sq10 = 3.1622776601683793320 !sqrt(10D0)
  
  integer :: tetraGaussN
  real(kind=precis), allocatable :: tetraGaussW(:)
  real(kind=precis), allocatable :: tetraGaussP(:,:)  

  integer :: trianGaussN
  real(kind=precis), allocatable :: trianGaussW(:)
  real(kind=precis), allocatable :: trianGaussP(:,:)  

  
  
  !!!!!basic variables================================================
  real (kind=precis), dimension(1:3,1:3) :: T = 0.0, Ti = 0.0
  real (kind=precis), dimension(1:8)     :: Ul = 0.0, Ur = 0.0 ! данные для распадника

  real (kind=precis), dimension(1:3) :: vert = (/ 0D0, 0D0, 1D0/), inp_center=(/0.43,0,0/), center1=(/1.0,0.0,0.0/), center2=(/1.0,0.0,0.0/)
  
  
  ! массив -- Q(i) -- локальный номер узла треугольной ячейки,
  ! стоящего после узла с локальным номером _i_
  integer, parameter, dimension(1:6):: Q = (/ 2, 3, 1, 2, 3, 1/)

  real (kind=precis) :: epsil = 10D-15
  real (kind=precis), parameter :: pi= 3.14159265358979323846
  real (kind=precis) :: gmm = 1.4 ! показатель адиабаты

  ! time and timestep parameters
  real (kind=precis) :: tau_out = 0.01   ! output timesteps
  real (kind=precis) :: t_start = 0.0    ! начальные момент времени
  real (kind=precis) :: t_final = 1.0    ! конечный момент времени
  real (kind=precis) :: t_current = 0.0  ! текущий момент времени
  integer            :: ntime  ! Номер текущего временного слоя

  ! auto timestep
  real (kind=precis) :: tau = 0.000001 ! Начальное значение временного шага
  real (kind=precis) :: CFL = 0.5
  integer            :: nmsr

  ! computation timer
  character(LEN=10) :: ddate, dtime, dzone
  integer, dimension(1:8) :: dvalues
  integer :: dhour, dmin, dsec, dmils, ddtime


  ! mesh descriptors
  real (kind=precis) :: hmin ! псевдошаг
  integer, dimension(1:3,1:4) :: surfInd = (/1,2,3, 1,3,4, 3,2,4, 2,1,4/) ! локальные номера
                                                                          ! вершин в поверхностях
  integer, dimension(1:3,1:4) :: locNeibSurf = (/2,3,4, 1,3,4, 1,2,4, 1,2,3/) ! 
  integer, allocatable :: aleph(:,:,:)
  integer, allocatable :: alephComb(:,:)
  
  integer:: nNodes                            ! количество узлов в сетке
  real (kind=precis),allocatable:: nodes(:,:) ! массив координат узлов
  integer :: nElems                         ! количество ячеек в сетке
  integer :: nSurfaces
  integer :: nElems_b                         ! количество ячеек в сетке

  integer,            allocatable :: elems(:,:) ! список КЭ (ячеек)
  real (kind=precis), allocatable :: elems_vol(:) ! elems_vol массив объемов ячеек
  real (kind=precis), allocatable :: elems_center(:,:) ! центр элемента
  integer,            allocatable :: elems_subdom(:) ! подобласть элемента
  integer,            allocatable :: elems_elems(:,:) ! соседи элемента
  integer,            allocatable :: elems_surfs(:,:) ! номера поверхностей элемента
  integer,            allocatable :: elems_surfs_normSign(:,:)
  integer,            allocatable :: elems_surfs_num(:,:) ! какой номер имеет элемент у поверхности с номером elems_surfs
  
  integer,            allocatable :: surfaces(:,:) ! список поверхностей ячеек
  real (kind=precis), allocatable :: surf_area(:) ! площади поверхностей
  real (kind=precis), allocatable :: surf_normal(:,:) ! внешние нормали поверхностей
  real (kind=precis), allocatable :: surf_center(:,:) ! центры поверхностей
  integer,            allocatable :: surf_elems(:,:) ! элементы, контактирующие через поверхность
  integer,            allocatable :: surf_elemsLocNum(:,:) !локальный номер поверхности в элементе


  integer,            allocatable :: bc_type(:) ! массив типов граничного условия
  integer,            allocatable :: bc_elems(:) ! номер элемента на границе
  integer,            allocatable :: bc_surf(:) ! номер поверхности на границе
  real (kind=precis), allocatable :: bc_surf_center(:,:)
  integer,            allocatable :: bc_elem_surf(:) ! локальный номер граничного ребра
  integer,            allocatable :: bc_nodes(:,:) ! номера точек на границе

  !!!!!solution=======================================================
  ! gas variables
  integer :: nComponents
  real (kind=precis), allocatable :: U(:,:,:)    ! массив для решения
  real (kind=precis), allocatable :: U_hat(:,:,:) ! массив для решения
  real (kind=precis), allocatable :: U_til(:,:,:) ! массив для решения
  logical, allocatable :: trouble(:)
  real (kind=precis), allocatable :: trInd(:)
  real (kind=precis) :: limitThreshold=1
  
  
  real (kind=precis), allocatable :: U_surf(:,:,:,:) ! массив для решения на поверхностях
  real (kind=precis), allocatable :: gas_surf_flux(:,:,:)
  
  ! service vars
  integer :: err, CFLChangeStep
  real (kind=precis) :: CFLInit
  logical :: CFLChanged
  real (kind=precis) :: CFLFactor = 0.75d0

  !!!!!file variables=================================================
  CHARACTER(LEN=*), PARAMETER:: root = './data' ! корневой каталог данных
  CHARACTER(LEN=*), PARAMETER:: fname_sol = 'sol' ! имя файла с суперэлементом
  CHARACTER(LEN=*), PARAMETER:: delim = '-'
  CHARACTER(LEN=*), PARAMETER:: extentionVTK = '.vtk'
  CHARACTER(LEN = max_filename_length):: current_sol_name = ''

  CHARACTER(LEN = max_filename_length):: &
      current_time_step_name = './data/last_time_step.dat'
  
end module names
