#!/bin/bash
#
# Retrieve energy levels and ionization energy from NIST Atomic Spectra Database
#

dir=./
echo ${dir}

lynx -dump -nolist "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all" \
    | sed -e '1,4d' \
    | head -n -3 \
    | sed -e 's/Atomic Symbol = D/Atomic Symbol = H/g' \
    | sed -e 's/Atomic Symbol = T/Atomic Symbol = H/g' \
	  > $dir/nist_isotope.txt
