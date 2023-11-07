#!/bin/sh

dir=./

# meta data
echo 'copy meta data...'
w3m -dump -cols 120 "https://hitran.org/docs/iso-meta/" \
    | sed -e '1,41d' \
    | sed -e 's/ × 10^/E/' \
    | sed -e 's/\[2\]/2/g' \
    | sed -e 's/\[3\]/3/g' \
    | sed -e 's/\[4\]/4/g' \
    | sed -e 's/\[6\]/6/g' \
    | sed -e 's/ ^/ /g' \
    | sed -e 's/\^/-/g' \
    | sed -e '/txt/s/H/-1H/g' \
    | sed -e '/txt/s/D/-2H/g' \
    | sed -e 's/ -/ /g' \
    | sed -e 's/I-2H/ID /g' \
    | sed -e 's/-+/_p/g' \
    | sed -e 's/mol--1/mol^-1/g' \
    | sed -e 's/:/ :/g' \
	  > $dir/hitran_meta.txt

sed -i -e 's/ × 10/E/' $dir/hitran_meta.txt

#partition functions
echo 'copy partition functions...'
mkdir -p Q/
for n in `seq 1 148` ; do
   fname="q${n}.txt"
   m=`echo ${n} | awk '{printf("%03d\n"),$1}'`
   oname="q${m}.txt"
   URL=https://hitran.org/data/Q/${fname}
   wget --spider ${URL} 2>&1 | grep "Not Found"
   if [ "$?" -ne 0 ]
   then
	echo ${URL} is fetched...
	wget -q -O ${dir}/Q/${fname} ${URL}
   fi
done
wget --quiet --no-directories -r --no-parent -A 'q*.txt' https://hitran.org/data/Q/ -P $dir/Q
