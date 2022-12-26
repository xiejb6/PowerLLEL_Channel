module mod_variables_postproc
    use mod_type,  only: fp
    use mod_utils, only: abort, setZero

    implicit none
    ! make everything public unless declared private
    public

    type grad_tensor_t
        real(fp), allocatable, dimension(:,:,:) :: ux, uy, uz
        real(fp), allocatable, dimension(:,:,:) :: vx, vy, vz
        real(fp), allocatable, dimension(:,:,:) :: wx, wy, wz
    end type

    ! main arrays
    real(fp), allocatable, dimension(:,:,:), save :: u, v, w
    real(fp), allocatable, dimension(:,:,:), save :: vor, q, lambda2
    type(grad_tensor_t), save :: vel_grad
    
    contains
    subroutine allocVariable(nhalo, sz, vartag, var)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        character(*), intent(in) :: vartag
        real(fp), allocatable, dimension(:,:,:), intent(out) :: var

        integer :: istatus

        allocate(var(1-nhalo(1):sz(1)+nhalo(2), 1-nhalo(3):sz(2)+nhalo(4), 1-nhalo(5):sz(3)+nhalo(6)), stat=istatus)
        if (istatus /= 0) call abort(101, "allocVariable: Out of memory when allocating variable <"//vartag//">!")
        call setZero(var)
    
        return
    end subroutine allocVariable

end module mod_variables_postproc