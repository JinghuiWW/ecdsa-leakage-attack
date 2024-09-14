n=256
s=3
m=86
x=0
e=0
f=1
t=24
g=0

folder_name="Instances/${n}_${s}"
mkdir -p "${folder_name}"

f1="${folder_name}/lines.txt"
f2="${folder_name}/lsb.txt"
f3="${folder_name}/sk.txt"

python solveECDSAfromLSB.py -n $n -s $s -m $m -x $x -e $e -t $t -g $g -f $f -f1 $f1 -f2 $f2 -f3 $f3
# mprof run solveECDSAfromLSB.py -n $n -s $s -m $m -x $x -e $e -t $t -g $g -f $f -f1 $f1 -f2 $f2 -f3 $f3
