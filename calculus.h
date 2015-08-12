#ifndef CALCULUS_H_INCLUDED
#define CALCULUS_H_INCLUDED

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include "consts.h"

//Romberg
template<typename F = complex<double> (*)(double),int INITIAL_NUM=20,int MAX_STEP=20>
complex<double> Intergate_Romberg(F const &func,double low,double up,double epsilon)
{
    complex<double> **R=new complex<double>* [INITIAL_NUM+1];
    for(int i=0;i<INITIAL_NUM+1;i++)
        R[i]=new complex<double> [INITIAL_NUM+1];

    //To be developed,should be triple storage
    /*
    R[0][0]
    R[0][1] R[1][1]
    R[0][2] R[1][2] R[2][2]
    R[0][3] R[1][3] R[2][3] R[3][3]
    ...                             ...

    */
    double h = up - low;
    //const complex<double> h(up-low,0);
    const complex<double> temp_low(low,0);
    const complex<double> temp_up(up,0);

    R[1][1]=(0.5*h*(func(up)+func(low)));

    for(int k=2;k<MAX_STEP+1;k++)
    {
        complex<double> sum = 0;
        for(int i=1;i<(1<<(k-2))+1;i++)
        {
            sum+=func(real(temp_low+(2*i-1)*h/(1<<(k-1))));

        }
        R[1][k]=0.5*(R[1][k-1]+(h/(1<<(k-2)))*sum);

        for(int j=2;j<k+1;j++)
        {
            R[j][k] = R[j-1][k]+(double)(1/((1<<(2*j-2))-1))*(R[j-1][k]-R[j-1][k-1]);
        }
        //errors
        if((abs(R[k][k])-abs(R[k-1][k-1])<epsilon)&&(abs(R[k][k])-abs(R[k-1][k-1])>-epsilon))
            return R[k][k];


    }
    return -10000;
}

template<typename F = complex<double> (*)(double,double),int INITIAL_NUM=256,int MAX_STEP = 20>
complex<double> Intergate_Romberg_first(F const &func,double y2,double low,double up,double epsilon)
{
    complex<double> **R=new complex<double>* [INITIAL_NUM+1];
    for(int i=0;i<INITIAL_NUM+1;i++)
        R[i]=new complex<double> [INITIAL_NUM+1];

    //To be developed,should be triple storage
    /*
    R[0][0]
    R[0][1] R[1][1]
    R[0][2] R[1][2] R[2][2]
    R[0][3] R[1][3] R[2][3] R[3][3]
    ...                             ...

    */
    double h = up - low;
    //const complex<double> h(up-low,0);
    const complex<double> temp_low(low,0);
    const complex<double> temp_up(up,0);

    R[1][1]=(0.5*h*(func(up,y2)+func(low,y2)));

    for(int k=2;k<MAX_STEP+1;k++)
    {
        complex<double> sum = 0;
        for(int i=1;i<(1<<(k-2))+1;i++)
        {
            sum+=func(real(temp_low+(2*i-1)*h/(1<<(k-1))),y2);

        }
        R[1][k]=0.5*(R[1][k-1]+(h/(1<<(k-2)))*sum);

        for(int j=2;j<k+1;j++)
        {
            R[j][k] = R[j-1][k]+((double)(1/((1<<(2*j-2))-1)))*(R[j-1][k]-R[j-1][k-1]);
        }
        //errors
        if((abs(R[k][k])-abs(R[k-1][k-1])<epsilon)&&(abs(R[k][k])-abs(R[k-1][k-1])>-epsilon))
            return R[k][k];


    }
    return -10000;
}


template<typename F = complex<double> (*)(double,double),int INITIAL_NUM=20,int MAX_STEP = 20>
complex<double> Intergate_Romberg_double(F const &f,double low1,double up1,double low2,double up2,double epsilon)
{
    complex<double> res;

    auto First_Kernel=[low1,up1,low2,up2,epsilon,f](double y2)
    {
        complex<double> res;
        res = Intergate_Romberg_first(f,y2,low1,up1,epsilon);
        return res;
    };

    res = Intergate_Romberg(First_Kernel,low2,up2,epsilon);
    return res;
}

#endif // CALCULUS_H_INCLUDED
