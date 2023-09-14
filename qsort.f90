! qsort.f90
! Reference: https://github.com/dongli/fortran-kdtree/blob/master/src/qsort_mod.F90

module qsort_module

contains

  recursive subroutine qsort(x, axis, left_, right_)

    real(8), intent(inout) :: x(:,:)
    integer, intent(in) :: axis
    integer, intent(in), optional :: left_
    integer, intent(in), optional :: right_

    integer :: left, right, i, j
    real(8) :: x0
    real(8), dimension(:), allocatable :: tmp

    left  = merge(left_, 1, present(left_))
    right = merge(right_, size(x,1), present(right_))

    x0 = x((left+right)/2, axis)
    i = left; j = right
    do
      do while (x(i,axis) < x0)
        i = i+1
      end do
      do while (x(j,axis) > x0)
        j = j-1
      end do
      if (i >= j) exit
      tmp = x(i,:); x(i,:) = x(j,:); x(j,:) = tmp
      i = i+1; j = j-1
    end do

    if (left < i-1)  call qsort(x, axis, left, i-1)
    if (j+1 < right) call qsort(x, axis, j+1, right)

  end subroutine qsort

end module qsort_module


