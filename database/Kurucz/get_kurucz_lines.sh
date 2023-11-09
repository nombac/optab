mkdir -p linelists

echo "Retrieving gfall08oct17.dat..."
wget http://kurucz.harvard.edu/linelists/gfnew/gfall08oct17.dat -O linelists/gfall08oct17.dat
ls linelists/gfall08oct17.dat > list_convert.txt
echo "Converting gfall08oct17.dat..."
../src/convert_lines_h5

echo "Retrieving gfpred26apr18.dat..."
wget http://kurucz.harvard.edu/linelists/gfpred/gfpred26apr18.dat.gz -O linelists/gfpred26apr18.dat.gz
gunzip linelists/gfpred26apr18.dat.gz
ls linelists/gfpred26apr18.dat > list_convert.txt
echo "Converting gfpred26apr18.dat..."
../src/convert_lines_h5
