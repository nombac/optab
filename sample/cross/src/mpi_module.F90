MODULE mpi_module

  USE ISO_FORTRAN_ENV
  USE MPI
  IMPLICIT NONE

  INTEGER(INT32), PROTECTED :: nprc, kprc, lprc, mprc, jprc
  INTEGER(INT32), PROTECTED :: myrk, myrk_k, myrk_l, myrk_m, myrk_j, unit
  INTEGER(INT32), PROTECTED :: mpi_grid_world, nprc_grid, myrk_grid
  INTEGER(INT32), PROTECTED :: mpi_line_world, nprc_line, myrk_line
  INTEGER(INT32), PROTECTED :: mpi_jconst_world, nprc_jconst, myrk_jconst
  INTEGER(INT32), PRIVATE :: mpi_tmp_world, nprc_tmp, myrk_tmp, error

CONTAINS
  

  SUBROUTINE mpi_first

    INTEGER(INT32) :: irnk
    INTEGER(INT32) :: k, l, m, j
    INTEGER(INT32), ALLOCATABLE :: itab(:,:,:,:)

    ! initialization
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprc, error)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrk, error)
    unit = 1000 + myrk
    
    ! process allocation
    CALL alloc_proc(nprc, kprc, lprc, mprc, jprc)
    
    ! communication table
    ALLOCATE(itab(-1:kprc, -1:lprc, -1:mprc, -1:jprc))
    itab = MPI_PROC_NULL

    irnk = 0
    DO j = 0, jprc-1
    DO m = 0, mprc-1
    DO l = 0, lprc-1
    DO k = 0, kprc-1
       itab(k,l,m,j) = irnk
       IF(myrk == irnk) THEN
          myrk_k = k
          myrk_l = l
          myrk_m = m
          myrk_j = j
       END IF
       irnk = irnk + 1
    END DO
    END DO
    END DO
    END DO

    DEALLOCATE(itab)

    ! j=const world
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, myrk_j, myrk, mpi_jconst_world, error)
    CALL MPI_COMM_SIZE(mpi_jconst_world, nprc_jconst, error)
    CALL MPI_COMM_RANK(mpi_jconst_world, myrk_jconst, error)

    ! (m,l)-world for collecting lines
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, myrk_k, myrk, mpi_line_world, error)
    CALL MPI_COMM_SIZE(mpi_line_world, nprc_line, error)
    CALL MPI_COMM_RANK(mpi_line_world, myrk_line, error)

    ! (m,k)-world
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, myrk_l, myrk, mpi_tmp_world, error)
    CALL MPI_COMM_SIZE(mpi_tmp_world, nprc_tmp, error)
    CALL MPI_COMM_RANK(mpi_tmp_world, myrk_tmp, error)

    ! k-world for calculating mean opacities in grid space
    CALL MPI_COMM_SPLIT(mpi_tmp_world, myrk_m, myrk, mpi_grid_world, error)
    CALL MPI_COMM_SIZE(mpi_grid_world, nprc_grid, error)
    CALL MPI_COMM_RANK(mpi_grid_world, myrk_grid, error)

    RETURN

  END SUBROUTINE mpi_first

  
  SUBROUTINE para_range(n1, n2, nprocs, irank, ista, iend)
    INTEGER(INT32), INTENT(IN) :: n1
    INTEGER(INT32), INTENT(IN) :: n2
    INTEGER(INT32), INTENT(IN) :: nprocs
    INTEGER(INT32), INTENT(IN) :: irank
    INTEGER(INT32), INTENT(OUT) :: ista
    INTEGER(INT32), INTENT(OUT) :: iend

    INTEGER(INT32) :: iwork1, iwork2

    iwork1 = (n2 - n1 + 1) / nprocs
    iwork2 = MOD(n2 - n1 + 1, nprocs)
    ista = irank * iwork1 + n1 + MIN(irank, iwork2)
    iend = ista + iwork1 - 1
    IF(iwork2 > irank) THEN
       iend = iend + 1
    END IF

    RETURN

  END SUBROUTINE para_range
  

  SUBROUTINE alloc_proc(nprc, kprc, lprc, mprc, jprc)
    INTEGER(INT32), INTENT(IN) :: nprc
    INTEGER(INT32), INTENT(OUT) :: kprc
    INTEGER(INT32), INTENT(OUT) :: lprc
    INTEGER(INT32), INTENT(OUT) :: mprc
    INTEGER(INT32), INTENT(OUT) :: jprc

    NAMELIST /mpi_decomp/ kprc, lprc, mprc, jprc
    
    OPEN(1,FILE='input/fort.5',STATUS='old')
    READ(1,NML=mpi_decomp)
    CLOSE(1)

!!$    WRITE(6,'(5(A,I4))') &
!!$         '<mpi information> nprc:', nprc, ' kprc:', kprc, ' lprc:', lprc, ' mprc:', mprc, ' jprc', jprc
!!$
    IF(nprc /= kprc * lprc * mprc * jprc) THEN
       WRITE(6,*) '(ERROR) SUBROUTINE alloc_proc: nprc /= kprc * lprc * mprc * jprc'
       CALL MPI_FINALIZE(error) ; STOP
    END IF
   
    RETURN
  END SUBROUTINE alloc_proc
  
  
END MODULE mpi_module
