g++ ./main.cpp -O3 -fopenmp -o main.out
sh mk-results.sh

COUNT=0
while true
do
    if [ $COUNT -lt 1000 ];
    then
        echo COUNT
        ./main.out
        COUNT=$(( $COUNT+1 ))
    else
        break
    fi
done