program PowerLLEL_Postprocess
    use mod_type, only: fp
    use mod_param_postproc
    use mod_variables_postproc
    use mod_mesh_postproc, only: initMesh, freeMesh
    use mod_dataio_postproc, only: inputVelField
    use mod_bc_postproc, only: updateBoundVel
    use mod_vortex_postproc, only: calcVelGradTensor, outVorticity, outQ, outLambda2
    use mod_hdf5, only: initIO, freeIO
    use mod_mpi, only: st, sz, myrank, ierr, initMPI, freeMPI
    
    implicit none
    
    include 'mpif.h'
    
    integer :: nt = -1
    character(80) :: prefix
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if (myrank == 0) then
        write(*,'(A)') "================================================================================"
        write(*,'(A)') "                               PowerLLEL_Postprocess                            "
        write(*,'(A)') "================================================================================"
        write(*,'(A)') "PowerLLEL_Postprocess.NOTE: Start postprocessing ..."
    endif

    CALL get_command_argument(1, prefix)
    if (len_trim(prefix) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Postprocess.NOTE: Please pass an inst. file prefix as an argument to PowerLLEL_Postprocess!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    if (myrank == 0) write(*,'(A)') "PowerLLEL_Postprocess.NOTE: Inst. file prefix is '"//trim(prefix)//"'"

    ! read input parameters from file
    call readInputParam('param_postproc.in')

    ! initialize MPI
    call initMPI(nx, ny, nz, (/'PP','PP','NN'/), p_row, p_col, nhalo)

    ! initialize parallel IO
    call initIO(MPI_COMM_WORLD)

    ! initialize mesh
    call initMesh()

    ! allocate variables
    call allocVariable(nhalo, sz, 'u', u)
    call allocVariable(nhalo, sz, 'v', v)
    call allocVariable(nhalo, sz, 'w', w)
    call allocVariable(nhalo_zero, sz, 'dudx', vel_grad%ux)
    call allocVariable(nhalo_zero, sz, 'dudy', vel_grad%uy)
    call allocVariable(nhalo_zero, sz, 'dudz', vel_grad%uz)
    call allocVariable(nhalo_zero, sz, 'dvdx', vel_grad%vx)
    call allocVariable(nhalo_zero, sz, 'dvdy', vel_grad%vy)
    call allocVariable(nhalo_zero, sz, 'dvdz', vel_grad%vz)
    call allocVariable(nhalo_zero, sz, 'dwdx', vel_grad%wx)
    call allocVariable(nhalo_zero, sz, 'dwdy', vel_grad%wy)
    call allocVariable(nhalo_zero, sz, 'dwdz', vel_grad%wz)

    ! read velocity fields
    call inputVelField(trim(prefix), nt, st, sz, nhalo, u, v, w)
    call updateBoundVel(nhalo, u, v, w)

    ! calculate velocity gradient tensor
    call calcVelGradTensor()

    ! write vortex fields
    if (out_vorticity) call outVorticity(trim(prefix), nt)
    if (out_q)         call outQ(trim(prefix), nt)
    if (out_lambda2)   call outLambda2(trim(prefix), nt)

    deallocate(u)
    deallocate(v)
    deallocate(w)
    deallocate(vel_grad%ux)
    deallocate(vel_grad%uy)
    deallocate(vel_grad%uz)
    deallocate(vel_grad%vx)
    deallocate(vel_grad%vy)
    deallocate(vel_grad%vz)
    deallocate(vel_grad%wx)
    deallocate(vel_grad%wy)
    deallocate(vel_grad%wz)

    call freeMesh()
    call freeIO()
    call freeMPI()

    if (myrank == 0) then
        write(*,'(A)') "********************************************************************************"
        write(*,'(A)') "PowerLLEL_Postprocess.NOTE: Postprocessing ends successfully!"
    endif

    call MPI_FINALIZE(ierr)
    
end program PowerLLEL_Postprocess
