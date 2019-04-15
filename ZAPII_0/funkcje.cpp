#include "funkcje.h"

void run()
{
    Matrix matA;//,matB,matX,l,u,c;

    matA.open("A.txt");
    matA.fill("A.txt");
    matA.dop();
    matA.save("SaveA.txt");
    Matrix A(0.5,10,0.0000785,250,100);
    A.save("TestA.txt");
    Matrix B(0.5,10,0.0000785,100,500,500,300);
    B.save("TestB.txt");

//    matB.open("B.txt");
//    matB.fill("B.txt");
//    matA.lu(l,u);
//    l.invert();
//    l.save("l.txt");
//    u.invert();
//    u.save("u.txt");
//    c.multiple(l,matB);
//    c.save("c.txt");
//    matX.multiple(u,c);
//    matX.save("SaveX.txt");


//    mat1.open("test.txt");
//    mat1.fill("test.txt");
//    mat2.open("savetest2.txt");
//    mat2.fill("savetest2.txt");
//    mat.multiple(mat1,mat2);
//    mat.save("saveTest3.txt");
//    mat1.transp();
//    mat1.save("saveTransp.txt");
//    mat1.lu(l,u);
//    l.save("saveL.txt");
//    u.save("saveU.txt");

}

