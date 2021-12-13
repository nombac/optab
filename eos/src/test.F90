PROGRAM test
  IMPLICIT NONE
  
  TYPE species_id
     INTEGER(4) :: id
     CHARACTER(LEN=32) :: species
  END type species_id
  
  TYPE(species_id) :: fastchem(2) = [species_id(1,'a'),&
                                     species_id(2,'b')]
  PRINT *, fastchem(1)%id, TRIM(fastchem(2)%species)
  
END PROGRAM test
