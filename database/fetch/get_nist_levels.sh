#!/bin/bash
#
# Retrieve energy levels and ionization energy from NIST Atomic Spectra Database
#

dir=levels
echo ${dir}

atom_min=1
atom_max=92
#atom_min=69
#atom_max=69
for atom in `seq ${atom_min} ${atom_max}` ; do
    charge_min=0
    charge_max=`expr ${atom} - 1`
#    charge_min=14
#    charge_max=14
    for charge in `seq ${charge_min} ${charge_max}` ; do
	name=nist_`printf %03d ${atom}`.`printf %02d ${charge}`
	file_tmp=${dir}/${name}.tmp
	file_states=${dir}/${name}.states
	file_ionize=${dir}/${name}.ionize

	URL="https://physics.nist.gov/cgi-bin/ASD/energy1.pl?encodedlist=XXT2&de=0&spectrum=Z%3D${atom}+${charge}&submit=Retrieve+Data&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&conf_out=on&term_out=on&level_out=on&g_out=on&temp="

	wget -q -O ${file_tmp} ${URL}
	echo ${name}

	if [ ${atom} -eq 1 ]; then
	    cat ${file_tmp} | grep -e \"1s\" -e \"\"\"\"\" | sed 's/,,/,0,/g' | sed '/Limit/d' | sed 's/\?//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/(//g' | sed 's/)//g' | sed 's/\&dagger\;//g'  | awk -F "," '{print $(NF-2),$(NF-1),$(NF-3),$(NF-4)}' | sed 's/\"=\"\"/ /g'| sed 's/\"\"\"//g' > ${file_states} #| egrep '[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?' > ${file_states}
	else
	    cat ${file_tmp} | sed -e '1d' | sed 's/,,/,0,/g' | sed '/Limit/d' | sed 's/\?//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/(//g' | sed 's/)//g' | sed 's/\&dagger\;//g'  | awk -F "," '{print $(NF-2),$(NF-1),$(NF-3),$(NF-4)}' | sed 's/\"=\"\"/ /g' | sed 's/\"\"\"//g' | awk 'gsub("[^0-9.]*","",$2)' | sed 's/\//\\/g' > ${file_states} #> ${file_states} #| egrep '[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?' > ${file_states}
	fi

	cat ${file_tmp} | grep -e Limit -e limit | sed -n 1P | sed 's/\?//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/(//g' | sed 's/)//g' | sed 's/\&dagger\;//g' | awk -F "," '{print $(NF-1)}'  | sed 's/\"=\"\"/ /g'| sed 's/\"\"\"//g' > ${file_ionize}

	rm ${file_tmp}
    done
done

touch ${dir}/levels_status.dat
