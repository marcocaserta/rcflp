# ./bin/rcflp -i data/orLibrary/cap101 -v 1 -t 1
# cat solution.txt >> summary.txt

# Nominal Version
./bin/rcflp -i data/orLibrary/cap101 -v 2 -t 1
cat solution.txt >> summary.txt

# Ellipsoidal
for omega in `seq 1 1 10`;
do
    echo 0.9 > parameters/paramsEllipsoidal.txt
    echo $omega >> parameters/paramsEllipsoidal.txt
    ./bin/rcflp -i data/orLibrary/cap101 -v 3 -t 1
    cat solution.txt >> summary.txt
done

# Box
for epsi in `seq 0.0 0.1 1`;
do
    echo $epsi > parameters/paramsBox.txt
    ./bin/rcflp -i data/orLibrary/cap101 -v 4 -t 1 -u 1
    cat solution.txt >> summary.txt
done

# Budget
for gamma in `seq 0.0 0.05 0.3`;
do
    echo 0.0 > parameters/paramsBudget.txt
    echo 1.0 >> parameters/paramsBudget.txt
    echo $gamma >> parameters/paramsBudget.txt
    echo 1 >> parameters/paramsBudget.txt
    ./bin/rcflp -i data/orLibrary/cap101 -v 4 -t 1 -u 2 -r 0
    cat solution.txt >> summary.txt
    # store B_l set for replicability
    cp support/cap101.budget support/cap101.budget.$gamma

    for epsi in `seq 0.1 0.1 1`;
    do
        for delta in `seq 0.8 0.05 0.95`;
        do
                echo $epsi > parameters/paramsBudget.txt
                echo $delta >> parameters/paramsBudget.txt
                echo $gamma >> parameters/paramsBudget.txt
                echo 1 >> parameters/paramsBudget.txt
                # NOTE: B_l is read from disk
                ./bin/rcflp -i data/orLibrary/cap101 -v 4 -t 1 -u 2 -r 1
                cat solution.txt >> summary.txt
            done
        done
done
