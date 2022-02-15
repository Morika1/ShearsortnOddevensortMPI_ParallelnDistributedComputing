Running Instructions:

0. Program name: testMpi (c file), input file name: Rectangles:dat.txt, output file name: Results.dat.txt
1. Open Ubuntu command line
2. Compile the files using: mpicc <fileName.c> -lm -o <fileName>
3. Execute with:mpiexec -n <numOfProcesses> ./<fileName>
4. the program will provide a list of id's after sorting them and writing them snakewise.The sorted list of id's will be printed both in command line and file results.dat
