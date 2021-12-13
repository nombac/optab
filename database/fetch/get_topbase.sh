#!/bin/bash
#
# Retrieve energy levels and photoionization cross sections from TOPbase
#

dir=./
#echo ${dir}

#for nz in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 16 18 20 26 ; do
for nz in 11 ; do
#    for ne in `seq 1 ${nz}` ; do
    for ne in 9 ; do

	ename=${dir}/elevel.`printf %02d ${nz}`.`printf %02d ${ne}`.dat
	URL="http://cdsweb.u-strasbg.fr/cgi-bin/topbase/topbase.sh?com=dt&ent=e&nz1=${nz}&nz2=${nz}&ne1=${ne}&ne2=${ne}&is1=1&is2=9&il1=0&il2=9&ip1=0&ip2=1&lv1=0&lv2=0&en1=0.0&en2=0.0&so=s2&te=on&gi=on"
	wget -O ${ename} ${URL}
	grep 'ERROR' ${ename}
	if [ "$?" -eq 0 ]
	then
	    echo "*** ERROR ***"
	    head ${ename}
	    rm ${ename}
	else
	    sed -i '1,8d' ${ename}
	fi

	xname=${dir}/xsectn.`printf %02d ${nz}`.`printf %02d ${ne}`.dat
	URL="http://cdsweb.u-strasbg.fr/cgi-bin/topbase/topbase.sh?com=dt&ent=p&nz1=${nz}&nz2=${nz}&ne1=${ne}&ne2=${ne}&is1=1&is2=9&il1=0&il2=9&ip1=0&ip2=1&lv1=0&lv2=0&en1=0.0&en2=0.0&so=s2"
	wget -O ${xname} ${URL}
	grep 'ERROR' ${xname}
	if [ "$?" -eq 0 ]
	then
	    echo "*** ERROR ***"
	    head ${xname}
	    rm ${xname}
	else
	    sed -i '1,8d' ${xname}
	fi
	
    done
done
