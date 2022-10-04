MODULE const_module

  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

  REAL(REAL64), PARAMETER :: m_ele             = 9.10938d-28
  REAL(REAL64), PARAMETER :: clight            = 2.99792d10
  REAL(REAL64), PARAMETER :: c                 = 2.99792d10
  REAL(REAL64), PARAMETER :: amu               = 1.660538921d-24
  REAL(REAL64), PARAMETER :: pi                = 3.1415926535898d0
  REAL(REAL64), PARAMETER :: hbar              = 1.0545726397789448d-27
  REAL(REAL64), PARAMETER :: h                 = 6.6260755d-27
  REAL(REAL64), PARAMETER :: k_bol             = 1.3806503d-16
  REAL(REAL64), PARAMETER :: alpha             = 1d0 / 137.035999084d0
  REAL(REAL64), PARAMETER :: e2                = alpha * hbar * clight
  REAL(REAL64), PARAMETER :: c2                = clight * h / k_bol
  REAL(REAL64), PARAMETER :: a_0               = hbar**2 / (m_ele * e2)
  REAL(REAL64), PARAMETER :: m_hyd             = 1.008d0 * amu
  REAL(REAL64), PARAMETER :: ene_hyd           = 0.5d0 * alpha**2 * m_ele * clight**2 ! [erg]
  REAL(REAL64), PARAMETER :: cross_Thomson     = 6.6524587321d-25 ! [cm^2]
  REAL(REAL64), PARAMETER :: megabarn_to_cm2   = 1d-18 ! [cm^2]
  REAL(REAL64), PARAMETER :: barn_to_cm2       = 1d-24 ! [cm^2]
  REAL(REAL64), PARAMETER :: ev_to_erg         = 1.602176634d-12
  REAL(REAL64), PARAMETER :: erg_to_ev         = 1d0 / ev_to_erg
  REAL(REAL64), PARAMETER :: cm_to_ang         = 1d8
  REAL(REAL64), PARAMETER :: wn_to_ev          = (h * clight) * erg_to_ev
  REAL(REAL64), PARAMETER :: crs_thomson       = 8d0 * pi / 3d0 * (alpha * hbar / (m_ele * clight))**2
  REAL(REAL64), PARAMETER :: h_electron_aff_ev = 0.754195d0
  REAL(REAL64), PARAMETER :: rydberg           = 0.5d0 * m_ele * (alpha * clight)**2
  REAL(REAL64), PARAMETER :: helium_ion_wavenum= 198310.66637d0 ! [cm^-1]
  REAL(REAL64), PARAMETER :: h2_ion_ev         = 15.425927d0
  REAL(REAL64), PARAMETER :: atm_to_bar        = 1.01325d0
  REAL(REAL64), PARAMETER :: hartree           = m_ele * clight**2 * alpha**2
  REAL(REAL64), PARAMETER :: ryd_to_wnum       = rydberg / (h * clight)
END MODULE const_module
