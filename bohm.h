#ifndef BOHM_H_INCLUDED
#define BOHM_H_INCLUDED

#include "calculus.h"
#include "basic_vec.h"
#include <fstream>

class Split
{
private:
    double L;//source position
    double h;//split height
    double lambda;//wave length
    double d;//splits spacing
    double omega;//split width
    double working_precision;
public:
    Split(double _d,double _omega,double _h,double _lambda,double _L,double _working_precision)
    {
        d      = _d;
        omega  = _omega;
        h      = _h;
        lambda = _lambda;
        L      = _L;
        working_precision = _working_precision;
    }
//-------------------------------------------------------------------------------------------------------------
//Fresnel-Huygens Model
//-------------------------------------------------------------------------------------------------------------
    //gamma
    complex<double> Gamma(double D)
    {
        complex<double> res(cos(2*PI/lambda*(L+D)/(L*D)),sin(2*PI/lambda*(L+D)/(L*D)));
        return res;
    }
    //C_Z
    complex<double> C_Z(double D)
    {
        complex<double> res;

        auto f = [this,D](double z)
        {
            complex<double> res(cos(2*PI/lambda)*(z*z*0.5/L+z*z*0.5/D),sin(2*PI/lambda)*(z*z*0.5/L+z*z*0.5/D));
            return res;
        };

        res = Intergate_Romberg(f,-h,h,working_precision);
        return res;
    }

    //K1
    //Kernel with one split open
    complex<double> K_Single(double x,double low,double up,double D)
    {
        complex<double> res;

        auto f = [x,low,up,this,D](double y)
        {
            complex<double> res(cos(2*PI/lambda*(y*y*0.5/L+(y-x)*(y-x)*0.5/D)),sin(2*PI/lambda*(y*y*0.5/L+(y-x)*(y-x)*0.5/D)));
            return res;
        };

        res = (-1/(lambda*lambda))*Intergate_Romberg(f,low,up,working_precision);
        return res;
    }

    //K2
    //Kernel with two split open
    complex<double> K_double(double x,double low1,double up1,double low2,double up2,double D)
    {
        complex<double> res;

        auto f=[x,low1,up1,low2,up2,this,D](double y1,double y2)
        {
            complex<double> res(sqrt(abs(y2-y1))*cos(2*PI/lambda*(y1*y1*0.5/L+abs(y2-y1)+(x-y2)*(x-y2)*0.5/D)),sqrt(abs(y2-y1))*sin(2*PI/lambda*(y1*y1*0.5/L+abs(y2-y1)+(x-y2)*(x-y2)*0.5/D)));
            return res;
        };

        res = Intergate_Romberg_double(f,low1,up1,low2,up2,working_precision);
        complex<double> t(cos(3*PI/4),sin(3*PI/4));
        res = 1/(sqrt(pow(lambda,2.5)))*t*res;
        return res;
    }

    //K
    //total Kernel
    complex<double> K(double x,double D)
    {
        complex<double> res;
            res = (K_Single(x,-d-0.5*omega,-d+0.5*omega,D)+K_Single(x,-0.5*omega,+0.5*omega,D)+K_Single(x,d-0.5*omega,d+0.5*omega,D)+K_double(x,-d-omega*0.5,-d+omega*0.5,-omega*0.5,omega*0.5,D)+K_double(x,-omega*0.5,omega*0.5,-d-omega*0.5,-d+omega*0.5,D)+K_double(x,-omega*0.5,omega*0.5,d-omega*0.5,d+omega*0.5,D)+K_double(x,d-omega*0.5,d+omega*0.5,-omega*0.5,omega*0.5,D)+K_double(x,-d-omega*0.5,-d+omega*0.5,d-omega*0.5,d+omega*0.5,D)+K_double(x,d-omega*0.5,d+omega*0.5,-d-omega*0.5,-d+omega*0.5,D));
            return res;
    }

    //abs(K)
    double K_abs(double x,double D)
    {
        double res;
        res = abs(K(x,D))/abs(K(0,D));
        return res;
    }


    vec<double> speed(double D,double time_step,double space_step)
    {
        vec<complex<double>>  nabla_K(1/space_step*0.5*(K(space_step,D)-K(-space_step,D)),1/space_step*(K(0,D+space_step)-K(0,D)),0);
        vec<double> cur_speed  = Const_h_ba/Const_e_m * imag((1/(abs(K(0,D))*abs(K(0,D))))*conj(K(0,D))*nabla_K);
        return cur_speed;
    }

    vec<double> next_timestep(vec<double> &cur_p,double time_step,double space_step)
    {
        //cout<<K(cur_p.x[0]+space_step,cur_p.x[1])-K(cur_p.x[0]-space_step,cur_p.x[1])<<endl;
        //cout<<K(cur_p.x[0]+space_step,cur_p.x[1]);
        cout<<K(cur_p.x[0]-space_step,cur_p.x[1]);
        vec<complex<double>>  nabla_K(1/space_step*0.5*(K(cur_p.x[0]+space_step,cur_p.x[1])-K(cur_p.x[0]-space_step,cur_p.x[1])),1/space_step*(K(cur_p.x[0],cur_p.x[1]+space_step)-K(cur_p.x[0],cur_p.x[1])),0);
        vec<double> cur_speed  = -Const_h_ba/Const_e_m * imag((1/(abs(K(cur_p.x[0],cur_p.x[1]))*abs(K(cur_p.x[0],cur_p.x[1]))))*conj(K(cur_p.x[0],cur_p.x[1]))*nabla_K);

        //cout<<cur_speed.x[0]<<endl;
        //ofstream file("test.txt");
        //if(cur_speed.x[1]>0) cout<<cur_speed<<endl;
        vec<double> next_pos = cur_p+time_step*cur_speed;
        return next_pos;
    }

    //Phi_Field calculate a Phi Field
    //in the range of
    //D_low to D_up on z axis
    //-X_RANGE to X_RANGE

    complex<double> **Kernel_Field(double D_low,double D_up,double X_RANGE,double X_STEP,double D_precision)
    {
        int X_NUM = (int)(floor((2*X_RANGE)/X_STEP)+1);
        int D_RANGE = (int)floor(((D_up-D_low)/D_precision))+1;
        complex<double> **res = new complex<double>* [X_NUM];
        for(int i=0;i<X_NUM;i++)
        {
            res[i] = new complex<double> [D_RANGE];
        }

        int counter_i,counter_j;
        counter_i=0;counter_j=0;

        ofstream file("field.txt");
        if(!file.is_open())
        {
            cout<<"Error"<<endl;
            exit(1);
        }

        for(double i=-X_RANGE;i<=X_RANGE;i = i+X_STEP)
        {//to calculate Phi Field on range -X_RANGE to X_RANGE with step length X_STEP
            counter_j = 0;
            for(double j = D_low;j<=D_up;j = j+D_precision)
            {//calculate the Phi field on range D_low to D_up
                res[counter_i][counter_j] = K(i,j);
                counter_j++;
            }
            counter_i++;
        }
        return res;
    }
};


#endif // BOHM_H_INCLUDED
