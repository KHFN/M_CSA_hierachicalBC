## M_CSA_hierachicalBC
Monitoring the Supply Chain of Controlled Substances with Privacy-Preserving Hierarchical Blockchain

# instruction

this code uses MITSUNARI Shigeo's MCL

install [herumi/mcl](https://github.com/herumi/mcl)

    git clone https://github.com/herumi/mcl;
    cd mcl;
    make -j4;
    
    mkdir build
    cd build
    cmake ..
    make

After clone this repository to the same directory containing mcl and
how to implement how to mainalgoritms: 

    g++ -c -o filename.o test.cpp -I../mcl/include -lmcl -L../mcl/lib;
    g++ -o a.exe filename.o -Iinclue -I../mcl/include -lmcl -L../mcl/lib
    ./a.exe

how to implement TxGen: 

    g++ -c -o filename.o test.cpp -I../mcl/include -lmcl -L../mcl/lib;
    g++ -o a.exe filename.o -Iinclue -I../mcl/include -lmcl -L../mcl/lib;
    ./a.exe
