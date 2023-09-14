! tree.f90
! Tree implementation in Fortran.
! Miroslav Broz (miroslav.broz@email.cz), Sep 14th 2023

!               [tree%root]                                 
!             /              \                              
!       [left]                [right]         ---> x, sorted
!      /      \              /       \                      
! [left]       [right] [left]         [right] ---> y, sorted

module tree_module
  use node_module, node_init => init, node_prn => prn
  implicit none

  type treetype
    integer :: k
    type(nodetype), pointer :: root => null()
  contains
    procedure :: init
    procedure :: prn
    procedure :: nearest_neighbor
  end type treetype

  double precision :: best
  double precision, dimension(2) :: best_q

contains

  subroutine init(self, k, points)
    class(treetype), intent(out) :: self
    integer, intent(in) :: k
    double precision, dimension(:,:), intent(inout) :: points
    self%k = k
    allocate(self%root)
    call build_tree(self, self%root, points, 0)
  end subroutine init

  recursive subroutine build_tree(self, node, points, depth)
    use qsort_module
    class(treetype), intent(in) :: self
    type(nodetype), pointer, intent(out) :: node
    double precision, dimension(:,:), intent(inout) :: points
    integer, intent(in) :: depth
    integer :: i, n, axis

    if (size(points,1) == 0) return
    allocate(node)
    axis = modulo(depth, self%k)+1
    call qsort(points, axis)
    n = size(points,1)
    i = n/2+1
    call node%init(points(i,:))
    if (i>1) call build_tree(self, node%left, points(:i-1,:), depth+1)
    if (i<n) call build_tree(self, node%right, points(i+1:,:), depth+1)
  end subroutine build_tree

  subroutine prn(self)
    class(treetype), intent(in) :: self
    call node_prn(self%root)
  end subroutine prn

  function nearest_neighbor(self, p)
    class(treetype), intent(in) :: self
    double precision, dimension(:), intent(in) :: p
    double precision, dimension(2) :: nearest_neighbor
    best = 1.d38
    call recursive_search(self, self%root, 0, p)
    nearest_neighbor = best_q
  end function nearest_neighbor

!   p +  [node.p]        
!       /        \       
! [left]          [right]
!                        
!             <------>   
!               diff     
!                        
! ---------------------->
!           axis = x or y
!
! if p in on the left from [node.p], open [left]
! only if diff < best, open [right]!
! otherwise, it can hardly be closer
! because left-node-right are closest, sorted in x (or y)

  recursive subroutine recursive_search(self, node, depth, p)
    class(treetype), intent(in) :: self
    type(nodetype), pointer, intent(in) :: node
    integer, intent(in) :: depth
    double precision, dimension(:), intent(in) :: p
    integer :: axis
    double precision :: tmp

    if (.not.associated(node)) return
    tmp = square_distance(node%p, p)
    if (tmp < best) then
      best = tmp
      best_q = node%p
    endif
    axis = modulo(depth, self%k)+1
    tmp = p(axis) - node%p(axis)
    if (tmp <= 0.d0) then
      call recursive_search(self, node%left, depth+1, p)
      if (tmp**2 < best) call recursive_search(self, node%right, depth+1, p)
    else
      call recursive_search(self, node%right, depth+1, p)
      if (tmp**2 < best) call recursive_search(self, node%left, depth+1, p)
    endif
  end subroutine recursive_search

  double precision function square_distance(a, b)
    double precision, dimension(:), intent(in) :: a, b
    integer :: i
    double precision :: s
    s = 0.d0
    do i = 1, size(a)
      s = s + (a(i)-b(i))**2
    enddo
    square_distance = s
  end function square_distance

  function brute_force(points, p)
    implicit none
    double precision, dimension(:,:), intent(in) :: points
    double precision, dimension(:), intent(in) :: p
    double precision, dimension(2) :: brute_force
    integer :: i, j
    double precision :: tmp, best
    best = 1.d38
    do i = 1, size(points,1)
      tmp = square_distance(p, points(i,:))
      if (tmp < best) then
        best = tmp
        j = i
      endif
    enddo
    brute_force = points(j,:)
  end function brute_force

end module tree_module


program test
  use tree_module
  implicit none
  double precision, dimension(2) :: p, q
  double precision, dimension(:,:), pointer :: points
  type(treetype), pointer :: tree
  integer :: i, n
  double precision :: tmp
  double precision :: t1, t2

  n = 100000
  allocate(points(n, 2))
  tmp = rand(-1)
  do i = 1, n
    points(i,1) = rand()
    points(i,2) = rand()
  enddo

  allocate(tree)
  call cpu_time(t1)
  call tree%init(2, points)
  call cpu_time(t2)
  write(*,*) 'tree%init  cpu_time = ', t2-t1

!  call tree%prn()

  p = [0.5d0, 0.5d0]
  call cpu_time(t1)
  q = brute_force(points, p)
  call cpu_time(t2)
  write(*,*) 'brute_force = ', q, ' cpu_time = ', t2-t1

  call cpu_time(t1)
  q = tree%nearest_neighbor(p)
  call cpu_time(t2)
  write(*,*) 'nearest_neighbor = ', q, ' cpu_time = ', t2-t1

end program test


