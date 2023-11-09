#!/bin/bash
#
# Retrieve energy levels and photoionization cross sections from TOPbase
#

dir=./atoms/
#echo ${dir}
mkdir -p ${dir}

nz_min=1
nz_max=92
for nz in `seq ${nz_min} ${nz_max}` ; do
#for nz in 11 ; do
    ne_min=0
    ne_max=`expr ${nz} - 1`
    for ne in `seq ${ne_min} ${ne_max}` ; do

	number=`printf %02d ${nz}``printf %02d ${ne}`
	ename=gf${number}.gam
	URL="http://kurucz.harvard.edu/atoms/${number}/${ename}"
	wget --spider ${URL} 2>&1 | grep "Not Found"
	if [ "$?" -eq 0 ]
	then
	    echo ${ename} does not exist...
	else
	    echo ${ename} is fetched...
	    wget -q -O ${dir}${ename} ${URL}
	fi
	
    done
done
