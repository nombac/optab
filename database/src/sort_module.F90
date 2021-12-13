MODULE sort_module
  USE ISO_FORTRAN_ENV

  PRIVATE 
  PUBLIC :: sorti

CONTAINS

  SUBROUTINE sorti(n,arr,indx)
    
    IMPLICIT NONE

    INTEGER(INT64), INTENT(IN) :: n
    REAL(REAL64), INTENT(IN) :: arr(:)
    INTEGER(INT64), INTENT(INOUT) :: indx(:)
    !************************************************************************
    !*  index-based sorting routine for s3r2t mu integration
    !*  adapted from the routine indexx by press etal (1986)
    !*           version 1.1 of 22/mar/91 changed by phh
    !*-- notes:
    !* original code taken from the routine 'indexx' by
    !* press etal: 1986, 'numerical recipes',
    !* cambridge university press
    !************************************************************************
    INTEGER(INT64) :: l, ir, indxt, i, j
    REAL(REAL64) :: q

    l = n / 2 + 1
    ir = n
    DO 
       IF(l > 1) THEN
          l = l - 1
          indxt = indx(l)
          q = arr(indxt)
       ELSE
          indxt = indx(ir)
          q = arr(indxt)
          indx(ir) = indx(1)
          ir = ir - 1
          IF(ir == 1)THEN
             indx(1) = indxt
             RETURN
          ENDIF
       ENDIF
       i = l
       j = l + l
       DO 
          IF(j <= ir) THEN
             IF(j  < ir) THEN
                IF(arr(indx(j)) < arr(indx(j+1))) j = j + 1
             ENDIF
             IF(q < arr(indx(j))) THEN
                indx(i) = indx(j)
                i = j
                j = j + j
             ELSE
                j = ir + 1
             ENDIF
          ELSE
             EXIT
          ENDIF
       END DO
       indx(i) = indxt

    END DO
    
  END SUBROUTINE sorti
    
END MODULE sort_module
