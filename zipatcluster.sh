ls | awk '{print "ls "$1"\/" " | tail -n 1 | cut -c 12-16"}' | /bin/sh > ../latest.txt
ls | xargs -n 1 > ../flist.txt
paste ../flist.txt ../latest.txt > ../zipix.txt
cat ../zipix.txt |awk '{print ""$1"\/Nmass.txt "$1"\/OrbitalElements"$2".txt " $1"\/PltAllGraph"$2".txt " $1"\/OrbitalElements00001.txt " $1"\/PltAllGraph00001.txt"}' | xargs -n 1 > ../ziplist.txt
zip ../archive.zip -@ < ../ziplist.txt
