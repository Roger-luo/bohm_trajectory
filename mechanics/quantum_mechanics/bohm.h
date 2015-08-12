#ifndef BOHM_H_INCLUDED
#define BOHM_H_INCLUDED
#include"../../algebra.h"
#include"../../const.h"
#include"../../Integration.h"
#include<fstream>

namespace bohm
{
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
    Split(double _d,double _omega,double _h,double _lambda,double _L,double _working_precision);
//-------------------------------------------------------------------------------------------------------------
//Fresnel-Huygens Model
//-------------------------------------------------------------------------------------------------------------
        //Gamma
        std::complex<double> Gamma(double D) const;
        //C_Z
        std::complex<double> C_Z(double D) const;
        //K
        //Kernel with one split open
        std::complex<double> K_Single(double x,double low,double up,double D) const;

        //K2
        //Kernel with two split open
        std::complex<double> K_double(double x,double low1,double up1,double low2,double up2,double D) const;
        //K
        //total Kernel
        std::complex<double> K(double x,double D) const;
        //K for double Split
        std::complex<double> K_double_split(double x,double D) const;
        //S for double Split
        double S_double(double x,double D) const;
        //S
        double S(double x,double D) const;
        linear_algebra::vector<double> nabla_S(double x,double D)const;
};//Split

//print S Field in range:
//x from range[0] to range[1] by step range[2]
//D from range[3] to range[4] by step range [5]
void S_Field_Print(const Split &K,const std::vector<double> &range);

class photon
{
private:
        //slit parameter
        double d;//width
        double D;//Distance of Screen
        double eta;//half-spacing of two splits
        //wave parameter
        double lambda;//wave length
        double k;//wave vector
        double E;
        double B;

    linear_algebra::vector<double> q;//position of a photon
public:
        photon();
        photon(linear_algebra::vector<double> initial_pos,double lambda,double _d,double _E,double _B,double _eta);
        linear_algebra::vector<double> next_timestep(double dt);
        linear_algebra::vector<double> cur_pos();
        linear_algebra::vector<double> wave_par();
        linear_algebra::vector<double> slit_par();
};//photon

void
photon_trajectory(double lambda,double d,double eta,double dt,double STOP_DISTANCE);

double
photon_tracer(photon &t_photon,double dt,const double MAX_DISTANCE);

double
bohm_trac_classic_trac
(
linear_algebra::vector<double> initial_pos,
double lambda,double _d,double _eta,
double MAX_DISTANCE,double error
);


#include"bohm_impl.h"
}
#endif // BOHM_H_INCLUDED
