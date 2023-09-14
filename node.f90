! node.f90
! Node-related stuff.
! Miroslav Broz (miroslav.broz@email.cz), Sep 14th 2023

module node_module
  implicit none

  type nodetype
    double precision, dimension(2) :: p
    type(nodetype), pointer :: left => null(), right => null()
  contains
    procedure :: init
  end type nodetype

contains

  subroutine init(self, p)
    class(nodetype), intent(out) :: self
    double precision, dimension(:), intent(in) :: p
    self%p = p
    ! Note: no allocation here, because of associated() test!
  end subroutine 

  recursive subroutine prn(node)
    type(nodetype), intent(in) :: node
    write(*,*) 'node%p = ', node%p
    if (associated(node%left)) call prn(node%left)
    if (associated(node%right)) call prn(node%right)
  end subroutine prn

end module node_module


