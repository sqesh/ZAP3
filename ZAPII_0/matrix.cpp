#include "matrix.h"

Matrix::Matrix(){}
Matrix::Matrix(int h, int w):h_(h),w_(w)
{
    if(h>1 && w>1)
    {
        p_Matrix_ = new double *[h];
        for(int i=0; i<h ; i++)
        {
            p_Matrix_[i] = new double [w];
        }
    }
}
Matrix::Matrix(double l, int n, double p, double k, double a)
{
    this->create(n,n);
    double dx = l/n;
    double pom1 = (k*p)/dx;
    double pom2 = 2*a*sqrt(3.14*p)*dx;


    for(int h=0; h<h_; h++)
    {
        for(int w=0; w<w_ ; w++)
        {
            if((h==0 && w==0) || (h==h_-1 && w== w_-1)) this->p_Matrix_[h][w] = 1;
            else if((w == h-1 || w == h+1)&&(h!=0&&h!=h_-1)) this->p_Matrix_[h][w] = -pom1;
            else if(w==h) this->p_Matrix_[h][w] = pom1 + pom2;
            else  this->p_Matrix_[h][w] = 0;
        }
    }
}
Matrix::Matrix(double l, int n, double p, double a, double ta, double tb, double to)
{
    this->create(n,1);
    double dx = l/n;
    double pom2 = 2*a*sqrt(3.14*p)*dx;
    for (int h=0 ; h<h_ ; h++)
    {
        if(h==0) p_Matrix_[h][0] = ta;
        else if(h==h_-1) p_Matrix_[h][0] = tb;
        else p_Matrix_[h][0] = pom2*to;
    }

}

void Matrix::create(int h, int w)
{
    h_ = h;
    w_ = w;
//    if(h>1 && w>1)
//    {
        p_Matrix_ = new double *[h];
        for(int i=0; i<h ; i++)
        {
            p_Matrix_[i] = new double [w];
        }
//    }
//    else if(h>1 && w==1)
//        p_Matrix_ = new double [h];
}
void Matrix::open(string file)
{
    ifstream inputFile;
    inputFile.open(&file[0]);
    if(inputFile.is_open())
    {
    /*   while(!inputFile.eof())
        {
            char znak(0);
            int h=0;
            int w=0;
            int tmp=0;
            for(int h=1; !inputFile.eof(); h++)
            {
                for(int w=1; znak != '\n'; w++)
                {
                    znak(0);
                    while(znak !='\t')inputFile.get(znak);
                    tmp = w;
                }
            }
        }*/
        int h=0, w=1;
        for(h=0; !inputFile.eof(); h++)
        {
            char znak(0);
            while(znak != '\n' && !inputFile.eof())
            {
                inputFile.get(znak);
                if(znak == '\t'&&h==1)w++;
            }
        }
        inputFile.close();
        create(h,w);
    }
}
void Matrix::fill(string file)
{
    ifstream inputFile;
    inputFile.open(&file[0]);
    if(inputFile.is_open())
    {
        for(int h=0; h<h_; h++)
        {
            for(int w=0; w<w_; w++)
            {
                inputFile >> p_Matrix_[h][w];
            }
        }
        inputFile.close();
    }
}
void Matrix::save(string file)
{
    ofstream outputFile;
    outputFile.open(&file[0]);
    if(outputFile.is_open())
    {
        for(int h=0; h<h_; h++)
        {
            for(int w=0; w<w_; w++)
            {
                outputFile << p_Matrix_[h][w]; // <<'\t';
                if(w<(w_-1))outputFile << '\t';
            }
            if(h<(h_-1))
            {
                outputFile << endl;
            }
        }
        outputFile.close();
    }
}
void Matrix::multiple(Matrix const &mat1, Matrix const &mat2)
{
    if(mat1.w_ == mat2.h_)
    {
        this->create(mat1.h_, mat2.w_);
        for(int h=0; h< this->h_ ; h++)
        {
            for(int w=0; w< this->w_;w++)
            {
                double sum=0;
                for(int i=0; i< mat1.w_;i++)
                {
                    sum += mat1.p_Matrix_[h][i] * mat2.p_Matrix_[i][w];
                }
                this->p_Matrix_[h][w] = sum;
            }
        }
    }
}

void Matrix::lu(Matrix &l, Matrix &u)
{
    if(this->h_ == this->w_)
    {
        l.create(this->h_, this->w_);
        u.create(this->h_, this->w_);
        l.det_ = 1;
        u.det_ = 1;
        for (int i = 0; i < this->h_; i++)
        {
            for (int j = 0; j < this->h_; j++)
            {
                if (j < i) l.p_Matrix_[j][i] = 0;
                else
                {
                    l.p_Matrix_[j][i] = this->p_Matrix_[j][i];
                    for (int k = 0; k < i; k++)
                    {
                        l.p_Matrix_[j][i] = l.p_Matrix_[j][i] - l.p_Matrix_[j][k] * u.p_Matrix_[k][i];
                    }
                }
            }
            for (int j = 0; j < this->h_; j++)
            {
                if (j < i)
                    u.p_Matrix_[i][j] = 0;
                else if (j == i)
                    u.p_Matrix_[i][j] = 1;
                else
                {
                    u.p_Matrix_[i][j] = this->p_Matrix_[i][j] / l.p_Matrix_[i][i];
                    for (int k = 0; k < i; k++)
                    {
                        u.p_Matrix_[i][j] = u.p_Matrix_[i][j] - ((l.p_Matrix_[i][k] * u.p_Matrix_[k][j]) / l.p_Matrix_[i][i]);
                    }
                }
            }
        }
        for(int i=0; i<l.h_; i++)
        {
            l.det_ *= l.p_Matrix_[i][i];
        }
        for(int i=0; i<u.h_; i++)
        {
            u.det_ *= u.p_Matrix_[i][i];
        }
        this->det_ = l.det_ * u.det_;
    }
}
void Matrix::transp()
{
    double pom;
    for(int h=0; h<this->h_; h++)
    {
        for(int w=0;w<h;w++)
        {
            pom = this->p_Matrix_[h][w];
            this->p_Matrix_[h][w] = this->p_Matrix_[w][h];
            this->p_Matrix_[w][h] = pom;
        }
    }
}
void Matrix::dop()
{
    Matrix pom2(this->h_,this->w_);
    for(int h=0; h<h_; h++)
    {
        for(int w=0; w<w_; w++)
        {
            Matrix pom,l,u;
            pom.create(h_-1,w_-1);
            for(int h1=0; h1<(h_-1); h1++)
            {
                for(int w1=0; w1<(w_-1); w1++)
                {
                    if(h1<h)
                    {
                        if(w1<w)
                        {
                            pom.p_Matrix_[h1][w1]=this->p_Matrix_[h1][w1];
                        }
                        else
                        {
                            pom.p_Matrix_[h1][w1]=this->p_Matrix_[h1][w1+1];
                        }
                    }
                    else
                    {
                        if(w1<w)
                        {
                            pom.p_Matrix_[h1][w1]=this->p_Matrix_[h1+1][w1];
                        }
                        else
                        {
                            pom.p_Matrix_[h1][w1]=this->p_Matrix_[h1+1][w1+1];
                        }
                    }
                }
            }
            string pomt = "pom[ ][ ].txt";  // debug
            pomt[4]=h+0x30;                 // debug
            pomt[7]=w+0x30;                 // debug
            pom.save(pomt);                 // debug

            pom.lu(l,u);

            string lt = "l[ ][ ].txt";      // debug
            lt[2]=h+0x30;                   // debug
            lt[5]=w+0x30;                   // debug
            string ut = "u[ ][ ].txt";      // debug
            ut[2]=h+0x30;                   // debug
            ut[5]=w+0x30;                   // debug
            l.save(lt);                     // debug
            u.save(ut);                     // debug

            pom2.p_Matrix_[h][w]=pom.det_;
        }
    }
    pom2.save("pom2.txt");                  //debug
}

void Matrix::invert()
{
    this->det_ = 1;
    for(int i=0; i<this->h_; i++)
    {
        this->det_ *= this->p_Matrix_[i][i];
    }
    if(this->det_ !=0)
    {
        this->dop();
        this->transp();
        double inv = 1 / this->det_;
        for(int h=0; h< this->h_ ; h++)
        {
            for(int w=0; w< this->w_;w++)
            {
                this->p_Matrix_[h][w] *= inv;
            }
        }
    }
}

Matrix::~Matrix()
{
//    if(h_>1 && w_>1)
//    {
        for(int i=0; i<h_ ; i++)
        {
            delete[]p_Matrix_[i];
        }
        delete[]p_Matrix_;
//    }
//    else if(h_>1 && w_==1) delete[]p_Matrix_;
}
