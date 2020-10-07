identify: main.o mcc.o
        g++ -fopenmp -g -o identify main.o mcc.o
main.o: main.cpp mcc.h
        g++ -fopenmp -g -c main.cpp
mcc.o: mcc.cpp mcc.h
        g++ -fopenmp -g -c mcc.cpp
clean:
        rm *.o identify
