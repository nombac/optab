
! Schreier 2018, MNRAS, 479, 3068

MODULE voigt_module

  DOUBLE PRECISION, PARAMETER ::  one   = 1.d0
  DOUBLE PRECISION, PARAMETER ::  two   = 2.d0
  DOUBLE PRECISION, PARAMETER ::  half    = 0.5d0
  DOUBLE PRECISION, PARAMETER ::  sqrtPi    = 1.7724538509055159d0 ! SQRT(pi)
  DOUBLE PRECISION, PARAMETER ::  recSqrtPi = one/sqrtPi

  PRIVATE
  PUBLIC :: hum1zpf16, hum1wei24, hum1wei32, voigt_bs
  
CONTAINS

  SUBROUTINE voigt_bs(x, y, prof)

    ! Voigt function approximation by Baschek & Sholtz (1982)

    USE ISO_FORTRAN_ENV
    USE const_module, ONLY : pi

    IMPLICIT NONE
    
    REAL(REAL64), INTENT(IN) :: x(:)
    REAL(REAL64), INTENT(IN) :: y
    REAL(REAL64), INTENT(OUT) :: prof(:)
    COMPLEX(REAL64), ALLOCATABLE :: p(:)
    REAL(REAL64), PARAMETER, DIMENSION(0:6) :: &
         a = [122.607931777104326d0, 214.382388694706425d0, 181.928533092181549d0, &
         93.155580458138441d0, 30.180142196210589d0, 5.912626209773153d0, 0.564189583562615d0], &
         b = [122.607931773875350d0, 352.730625110963558d0, 457.334478783897737d0, &
         348.703917719495792d0, 170.354001821091472d0, 53.992906912940207d0, 10.479857114260399d0], &
         c = [0.5641641d0, 0.8718681d0, 1.474395d0, -19.57862d0, 802.4513d0, -4850.316d0, 8031.468d0]
    INTEGER(INT32) :: n

    ALLOCATE(p(LBOUND(x,1):UBOUND(x,1)))
    
    p(:) = CMPLX(y, -x(:), KIND(0d0))
    
    DO n = LBOUND(x,1), UBOUND(x,1)
!!$       IF(x(n) + y > 15d0) THEN
!!$          prof(n) = (1d0 / SQRT(pi)) * y / (y**2 + x(n)**2)
!!$       ELSE IF(x(n) >= 2.5d0 .AND. y <= 1d-3) THEN
       IF(x(n) >= 2.5d0 .AND. y <= 1d-1) THEN
          prof(n) = DBLE(EXP(-x(n)**2) + y**2 * (1d0 - 2d0 * x(n)**2) * EXP(-x(n)**2) + &
               y * x(n)**(-2) * (c(0) + x(n)**(-2) * (c(1) + x(n)**(-2) * (c(2) + &
               x(n)**(-2) * (c(3) + x(n)**(-2) * (c(4) + x(n)**(-2) * (c(5) + x(n)**(-2) * c(6))))))))
       ELSE
          prof(n) = DBLE((a(0) + p(n) * (a(1) + p(n) * (a(2) + &
               p(n) * (a(3) + p(n) * (a(4) + p(n) * (a(5) + p(n) * a(6))))))) &
               /(b(0) + p(n) * (b(1) + p(n) * (b(2) + &
               p(n) * (b(3) + p(n) * (b(4) + p(n) * (b(5) + p(n) * (b(6) + p(n)))))))))
       END IF
!!$       IF(prof(n) <= 0d0) THEN
!!$          PRINT *, '******** negative prof:', prof(n), x(n), y
!!$          STOP
!!$       END IF
    END DO

    DEALLOCATE(p)
    
    RETURN
  END SUBROUTINE voigt_bs

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                                !!
    SUBROUTINE hum1wei24  (nx, x,y, w)
!!                                                                                                                                !!
!!    complex probability function for complex argument Z=X+iY                                                                    !!
!!    real part = voigt function K(x,y)                                                                                           !!
!!                                                                                                                                !!
!!  large x+y     J  Humlicek, JQSRT 27, 437, 1982                                                                                !!
!!  small x+y:    J.A.C. Weideman,  SIAM J. Numer. Anal. 31 (1994) pp. 1497-1518,  equation (38.I) and table I                    !!
!!                                                                                                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    USE constants,  ONLY: half, one, two, recSqrtPi
    IMPLICIT none

    INTEGER, PARAMETER ::   kdp=KIND(one)

    INTEGER, INTENT(IN) ::             nx
    DOUBLE PRECISION, INTENT(IN)  ::   x(0:nx), y
    COMPLEX(KIND=kdp), INTENT(OUT) ::  w(0:nx)

!   "Weideman" constants
    INTEGER, PARAMETER ::           n=24
    DOUBLE PRECISION, PARAMETER ::  l=4.1195342878142354d0 ! l=sqrt(n/sqrt(2.))  ! L = 2**(-1/4) * N**(1/2)
    DOUBLE PRECISION, PARAMETER ::  a(n) =  &
                        (/ -1.5137461654527820d-10,  4.9048215867870488d-09,  1.3310461806370372d-09, -3.0082822811202271d-08, &
                           -1.9122258522976932d-08,  1.8738343486619108d-07,  2.5682641346701115d-07, -1.0856475790698251d-06, &
                           -3.0388931839840047d-06,  4.1394617248575527d-06,  3.0471066083243790d-05,  2.4331415462641969d-05, &
                           -2.0748431511424456d-04, -7.8166429956142650d-04, -4.9364269012806686d-04,  6.2150063629501763d-03, &
                            3.3723366855316413d-02,  1.0838723484566792d-01,  2.6549639598807689d-01,  5.3611395357291292d-01, &
                            9.2570871385886788d-01,  1.3948196733791203d+00,  1.8562864992055408d+00,  2.1978589365315417d+00 /)

!   humlicek prbFct region I bounds
    DOUBLE PRECISION, PARAMETER ::    s15=15.d0
    DOUBLE PRECISION ::               x12, x21

!   further variables
    INTEGER ::                        i
    COMPLEX(KIND=kdp) ::              t, recLmZ  !? yy5, yi2

!   ================================================================================================================================
    x12 = y - s15 !        left wing -- center
    x21 = -x12    ! 15-y   center -- right wing
!?  yy5 = y*y + half
!?  yi2 = CMPLX(zero,two*y)

    IF (y>s15 .OR. x(0)>x21 .OR. x(nx)<x12) THEN
!       ----------------------------------------------------------------------------------------------------------------------------
!       all points are in Humlicek's region I
        DO i=0,nx
            t    = CMPLX(y,-x(i),kdp)
            w(i) = (recSqrtPi*t) / (half + t*t)
!?          w(i) = recSqrtPi * CMPLX(y,-x(i)) / (yy5-(x(i)+yi2)*x(i))   ! no big difference between these versions !
        END DO
!       ----------------------------------------------------------------------------------------------------------------------------
    ELSE
!       ----------------------------------------------------------------------------------------------------------------------------
        DO i=0,nx
            IF (ABS(x(i))>x21) THEN
                t    = CMPLX(y,-x(i),kdp)
                w(i) = (recSqrtPi*t) / (half + t*t)
            ELSE
                recLmZ = one / CMPLX(l+y,-x(i),kdp)
                t      = CMPLX(l-y,x(i),kdp) * recLmZ
                w(i) = recLmZ  *  (recSqrtPi + two*recLmZ * &
                (a(24)+(a(23)+(a(22)+(a(21)+(a(20)+(a(19)+(a(18)+(a(17)+(a(16)+(a(15)+(a(14)+(a(13)+(a(12)+(a(11)+(a(10)+(a(9)+ &
                (a(8)+(a(7)+(a(6)+(a(5)+(a(4)+(a(3)+(a(2)+a(1)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t))
            END IF
        END DO
!       ----------------------------------------------------------------------------------------------------------------------------
    END IF
!   ================================================================================================================================
END SUBROUTINE hum1wei24


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                                !!
    SUBROUTINE hum1zpf16  (nx, x,y, w)
!!                                                                                                                                !!
!!    complex probability function for complex argument Z=X+iY                                                                    !!
!!    real part = voigt function K(x,y)                                                                                           !!
!!                                                                                                                                !!
!!  large x+y     J  Humlicek, JQSRT 27, 437, 1982       asymptotic region I approximation                                        !!
!!  small x+y:    J  Humlicek, JQSRT 21, 309-313, 1979   16-term rational approximation Eq. (6)                                   !!
!!                                                                                                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   USE constants,  ONLY: half, one, recSqrtPi

    IMPLICIT none

    INTEGER, PARAMETER ::   kdp=KIND(one)

    INTEGER, INTENT(IN) ::             nx
    DOUBLE PRECISION, INTENT(IN)  ::   x(0:nx), y
    COMPLEX(KIND=kdp), INTENT(OUT) ::  w(0:nx)

!   humlicek prbFct region I bounds
    DOUBLE PRECISION, PARAMETER ::    s15=15.d0
    DOUBLE PRECISION ::               x12, x21

!   Humlicek zpf16 constants
    COMPLEX(KIND=kdp), PARAMETER ::  a(16) =          &
                [(41445.0374210222, 0.0),   &
                 (0.0, -136631.072925829),  &
                 (-191726.143960199, 0.0),  &
                 (0.0, 268628.568621291),   &
                 (173247.907201704,  0.0),  &
                 (0.0, -179862.56759178),   &
                 (-63310.0020563537, 0.0),  &
                 (0.0, 56893.7798630723),   &
                 (11256.4939105413,  0.0),  &
                 (0.0, -9362.62673144278),  &
                 (-1018.67334277366, 0.0),  &
                 (0.0, 810.629101627698),   &
                 (44.5707404545965, 0.0),   &
                 (0.0, -34.5401929182016),  &
                 (-0.740120821385939, 0.0), &
                 (0.0, 0.564189583547714)]
    DOUBLE PRECISION, PARAMETER ::  b(8) =  &
               [ 7918.06640624997,         &
                 -126689.0625,        &
                 295607.8125,        &
                 -236486.25,        &
                 84459.375,        &
                 -15015.0,        &
                 1365.0,        &
                 -60.0]

!   further variables
    INTEGER ::                        i
    COMPLEX(KIND=kdp) ::              t, z,zz, numer, denom

!   ================================================================================================================================
    x12 = y - s15 !        left wing -- center
    x21 = -x12    ! 15-y   center -- right wing

    IF (y>s15 .OR. x(0)>x21 .OR. x(nx)<x12) THEN
!       ----------------------------------------------------------------------------------------------------------------------------
!       all points are in Humlicek's region I
        DO i=0,nx
            t    = CMPLX(y,-x(i),kdp)
            w(i) = (recSqrtPi*t) / (half + t*t)
!?          w(i) = recSqrtPi * CMPLX(y,-x(i)) / (yy5-(x(i)+yi2)*x(i))   ! no big difference between these versions !
        END DO
!       ----------------------------------------------------------------------------------------------------------------------------
    ELSE
!       ----------------------------------------------------------------------------------------------------------------------------
        DO i=0,nx
            IF (ABS(x(i))>x21) THEN
                t    = CMPLX(y,-x(i),kdp)
                w(i) = (recSqrtPi*t) / (half + t*t)
            ELSE
                z     = CMPLX(x(i), y+1.31183, kdp)
                zz    = z*z
                numer = ((((((a(16)*z+a(15))*z+a(14))*z+a(13))*z+a(12))*z+a(11))*z+a(10))*z+a(9)
                numer = ((((((((numer*z+a(8))*z+a(7))*z+a(6))*z+a(5))*z+a(4))*z+a(3))*z+a(2))*z+a(1)) 
                denom = b(1)+(b(2)+(b(3)+(b(4)+(b(5)+(b(6)+(b(7)+b(8)*zz)*zz)*zz)*zz)*zz)*zz)*zz
                w(i)  = numer/denom
            END IF
        END DO
!       ----------------------------------------------------------------------------------------------------------------------------
    END IF
!   ================================================================================================================================
END SUBROUTINE hum1zpf16


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                                !!
    SUBROUTINE hum1wei32  (nx, x,y, w)
!!                                                                                                                                !!
!!    complex probability function for complex argument Z=X+iY                                                                    !!
!!    real part = voigt function K(x,y)                                                                                           !!
!!                                                                                                                                !!
!!  large x+y     J  Humlicek, JQSRT 27, 437, 1982                                                                                !!
!!                Kuntz, Ruyten, Wells: rewrite to real arithmetic                                                                !!
!!  small x+y:    J.A.C. Weideman,  SIAM J. Numer. Anal. 31 (1994) pp. 1497-1518,  equation (38.I) and table I                    !!
!!                                                                                                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IMPLICIT none
    
    INTEGER, PARAMETER ::   kdp=KIND(one)

    INTEGER, INTENT(IN) ::           nx
    DOUBLE PRECISION, INTENT(IN) ::  x(0:nx), y
    COMPLEX(kdp), INTENT(OUT) ::     w(0:nx)

!   "Weideman" constants
    INTEGER, PARAMETER ::           n=32
    DOUBLE PRECISION, PARAMETER ::  l=4.7568284600108841 ! l=sqrt(n/sqrt(2.))  ! L = 2**(-1/4) * N**(1/2)
    DOUBLE PRECISION, PARAMETER ::  a(n) =  &
                        (/ -1.3031797863050087e-12,  3.7424975634801558e-12,  8.0314997274316680e-12, -2.1542618451370239e-11, &
                           -5.5441820344468828e-11,  1.1658312885903929e-10,  4.1537465934749353e-10, -5.2310170266050247e-10, &
                           -3.2080150458594088e-09,  8.1248960947953431e-10,  2.3797557105844622e-08,  2.2930439569075392e-08, &
                           -1.4813078891201116e-07, -4.1840763666294341e-07,  4.2558331390536872e-07,  4.4015317319048931e-06, &
                            6.8210319440412389e-06, -2.1409619200870880e-05, -1.3075449254548613e-04, -2.4532980269928922e-04, &
                            3.9259136070122748e-04,  4.5195411053501429e-03,  1.9006155784845689e-02,  5.7304403529837900e-02, &
                            1.4060716226893769e-01,  2.9544451071508926e-01,  5.4601397206393498e-01,  9.0192548936480144e-01, &
                            1.3455441692345453e+00,  1.8256696296324824e+00,  2.2635372999002676e+00,  2.5722534081245696e+00 /)

!   Humlicek prbFct region I bounds
    DOUBLE PRECISION, PARAMETER ::  s15=15.d0
    DOUBLE PRECISION ::             x12, x21

!   further variables
    INTEGER ::                      i
    COMPLEX(KIND=kdp) ::            t, recLmZ

!   ================================================================================================================================
    x12 = y - s15 !        left wing -- center
    x21 = -x12    ! 15-y   center -- right wing

    IF (y>s15 .OR. x(0)>x21 .OR. x(nx)<x12) THEN
!       ----------------------------------------------------------------------------------------------------------------------------
!       all points are in Humlicek's region I
        DO i=0,nx
            t    = CMPLX(y,-x(i),kdp)
            w(i) = (recSqrtPi*t) / (half + t*t)
        END DO
!       ----------------------------------------------------------------------------------------------------------------------------
    ELSE
!       ----------------------------------------------------------------------------------------------------------------------------
        DO i=0,nx
            IF (ABS(x(i))>x21) THEN
                t    = CMPLX(y,-x(i),kdp)
                w(i) = (recSqrtPi*t) / (half + t*t)
            ELSE
                recLmZ  = one / CMPLX(l+y,-x(i),kdp)
                t       = CMPLX(l-y,x(i),kdp) * recLmZ
                w(i) = recLmZ  *  (recSqrtPi + two*recLmZ * &
                  (a(32)+(a(31)+(a(30)+(a(29)+(a(28)+(a(27)+(a(26)+(a(25)+(a(24)+(a(23)+(a(22)+(a(21)+(a(20)+(a(19)+(a(18)+(a(17)+ &
                  (a(16)+(a(15)+(a(14)+(a(13)+(a(12)+(a(11)+(a(10)+(a(9)+(a(8)+(a(7)+(a(6)+(a(5)+(a(4)+(a(3)+(a(2)+a(1) &
                  *t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t)*t) )
            END IF
        END DO
!       ----------------------------------------------------------------------------------------------------------------------------
    END IF
!   ================================================================================================================================
END SUBROUTINE hum1wei32


END MODULE voigt_module
