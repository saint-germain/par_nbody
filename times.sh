# run inside sim folder
ls | awk '{print "tail -n 2 "$1"\/Nmass.txt"}'  | /bin/sh | awk '{print $1}' > ../endtime.txt
ls > ../flist.txt
paste ../flist.txt ../endtime.txt > ../maxtimes.txt
rm ../endtime.txt
rm ../flist.txt
