module gffree_module

  private
  public :: gffree

contains
!c
!c
!c
!c--- free free gaunt factors

!c       double precision function gffree(gam22,t,u)
       double precision function gffree(gam22,u)
!c      -------------------------------------------
!c
!c   ferland's fabulous functional fits
!c-- to karzas & latter ap. j. 1966
!c
!*     gam22 1.58e5*z*z/t
!*     t     temperature ( k)
!*     u     h*f/(k*t)
       implicit none
!c
      real*8 :: A
      real*8 :: B
      real*8 :: C
      real*8 :: COEFF
      real*8 :: D
      real*8 :: FRAC
      real*8 :: GAM22
      integer :: I
      integer :: J
      integer :: K
      integer :: M
      real*8 :: SUM1
      real*8 :: SUM2
      real*8 :: ULOG
!c
       double precision gam2
!c       double precision t
       double precision u
       dimension coeff(28), a(7)
!c
       data coeff/ &
            1.102d0       ,-0.1085d0     ,0.09775d0     ,-0.01125d0,     &
            1.2d0         ,-0.24016667d0 ,0.07675d0     ,-0.01658333d0  ,&
            1.26d0        ,-0.313166667d0,0.15075d0     ,0.00241667d0   ,&
            1.29d0        ,-0.4518333d0  ,0.12925d0     ,0.00258333d0   ,&
            1.27d0        ,-0.579d0      ,0.092d0       ,-0.003d0       ,&
            1.16d0        ,-0.707333d0   ,0.112d0       ,0.0053333333d0 ,&
            0.883d0       ,-0.76885d0    ,0.190175d0    ,0.022675d0     /
       data a/100.d0, 10.d0, 3.d0, 1.d0, 0.3d0, 0.1d0, 0.001d0/
!c
       save a,coeff
!c
!c      u = 1.44e+8 / (wl*t)
       ulog = log10(u)
!c      gam2 = 1.58e+5 * z*z/t
       gam2 = gam22 ! gam2 later changed!
       if (gam2.gt.a(7)) go to 10
         i = 7
         j = 7
         k = 7
         frac = 0.5d0
         go to 60
   10  continue
       if (gam2.lt.a(1)) go to 20
         i = 1
         j = 1
         k = 1
         frac = 0.5d0
         go to 60
   20  continue
       do 30 i = 2, 7
           if (gam2.gt.a(i)) go to 40
   30  continue
   40  continue
       k = i - 1
!   50  continue
       b = log10(a(k))
       c = log10(a(i))
       gam2 = log10(gam2)
       frac = abs ((gam2-b) / (b-c))
   60  continue
       k = (k-1)*4
       sum1 = coeff(k+1)
       d = 1.0d0
       do 70 m = 2, 4
           d = d*ulog
           sum1 = sum1 + coeff(k+m)*d
   70  continue
       sum1 = sum1 * (1.0d0 - frac)
       i = (i-1)*4
       sum2 = coeff(i+1)
       d = 1.0d0
       do 80 m = 2, 4
           d = d*ulog
           sum2 = sum2 + coeff(i+m)*d
   80  continue
       sum2 = sum2 * frac
!c
       gffree = sum1 + sum2

     end function gffree

     end module gffree_module
