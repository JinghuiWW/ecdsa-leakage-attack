n=256
s=3
m1=90
m2=0
m3=0
m4=512
x=0
e=0
t=24
g=0

folder_name="minerva-data/athena/${n}_${s}"
f1="${folder_name}/lines.txt"
f2="${folder_name}/msb.txt"
f3="${folder_name}/sk.txt"

python solveECDSAfromMSB.py -n $n -s $s -m1 $m1 -m2 $m2 -m3 $m3 -m4 $m4 -x $x -e $e -t $t -g $g -f1 $f1 -f2 $f2 -f3 $f3 
# mprof run solveECDSAfromMSB.py -n $n -s $s -m1 $m1 -m2 $m2 -m3 $m3 -m4 $m4 -x $x -e $e -t $t -g $g -f1 $f1 -f2 $f2 -f3 $f3 

