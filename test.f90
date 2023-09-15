! test.f90
! Test tree.
! Miroslav Broz (miroslav.broz@email.cz), Sep 15th 2023

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


