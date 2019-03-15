# ./bin/rcflp -i data/orLibrary/cap101 -v 1 -t 1
# cat solution.txt >> summary.txt
#
./bin/rcflp -i data/orLibrary/cap101 -v 2 -t 1
cat solution.txt >> summary.txt
#
./bin/rcflp -i data/orLibrary/cap101 -v 3 -t 1
cat solution.txt >> summary.txt

for epsi in `seq 0.0 0.1 1`;
do
    echo $epsi > parameters/paramsBox.txt
    ./bin/rcflp -i data/orLibrary/cap101 -v 4 -t 1 -u 1
    cat solution.txt >> summary.txt
done
for epsi in `seq 0.0 0.1 1`;
do
    for delta in `seq 0.0 0.1 1.0`;
    do
        for gamma in `seq 0.0 0.05 0.3`;
        do
            echo $epsi > parameters/paramsBudget.txt
            echo $delta >> parameters/paramsBudget.txt
            echo $gamma >> parameters/paramsBudget.txt
            echo 1 >> parameters/paramsBudget.txt
            ./bin/rcflp -i data/orLibrary/cap101 -v 4 -t 1 -u 2
            cat solution.txt >> summary.txt
        done
    done
done
