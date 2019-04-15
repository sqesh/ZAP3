#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <stdio.h>
using namespace std;

class Matrix
{
private:
    double **p_Matrix_;
    double det_;
    int h_;
    int w_;


public:
    Matrix();
    Matrix(int h,int w);
    Matrix(double l, int n, double p, double k, double a);
    Matrix(double l, int n, double p, double a, double ta, double tb, double to);
    void create(int h, int w);
    void open(string file);
    void fill(string file);
    void save(string file);
    void multiple(Matrix const &mat1, Matrix const &mat2);
    void lu(Matrix &l, Matrix &u);
    void transp();
    void dop();
    void invert();
    ~Matrix();
};

#endif // MATRIX_H

//Matrix mat = new Matrix(3,4);
//mat.Oblicz
