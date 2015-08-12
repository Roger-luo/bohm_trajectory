#ifndef CLASSIC_DYNAMICS_H_INCLUDED
#define CLASSIC_DYNAMICS_H_INCLUDED

#include"../algebra.h"
#include<functional>

class particle
{
private:
linear_algebra::vector<double> velocity;
linear_algebra::vector<double> position;
double                                              static_mass;
double                                              cur_mass;
double                                              charge;
double                                              spin;

public:
//constructor & destructor
        particle();
        particle(const linear_algebra::vector<double> &initial_pos,const double &static_m);
        particle(const linear_algebra::vector<double> &initial_pos,const linear_algebra::vector<double> &initial_vel,const double &static_m);
        particle(const linear_algebra::vector<double> &initial_pos,const linear_algebra::vector<double> &initial_vel,const double &static_m,const double &_charge);
        particle(const linear_algebra::vector<double> &initial_pos,const linear_algebra::vector<double> &initial_vel,const double &static_m,const double &_charge,const double &_spin);
        ~particle();

//other
        linear_algebra::vector<double> cur_pos();
        linear_algebra::vector<double> cur_vel();
        double cur_m() const;
        double static_m() const;
        void reset_pos(const linear_algebra::vector<double> &_cur);
        void reset_vel(const linear_algebra::vector<double> &_cur);
};//particle

class classic_dynamic:public particle
{
public:
        classic_dynamic(const linear_algebra::vector<double> &initial_pos,const double &static_m):particle(initial_pos,static_m){};
        classic_dynamic(const linear_algebra::vector<double> &initial_pos,const linear_algebra::vector<double> &initial_vel,const double &_static_m):
        particle(initial_pos,initial_vel,_static_m){}
        //template<typename Force_Type>
        //linear_algebra::vector<double> next_timestep(const std::function<Force_Type(linear_algebra::vector<double> &)> &f,const double &dt);
        void next_timestep(const linear_algebra::vector<double> vel,const double dt);
};

#include"classic_dynamics_impl.h"
#endif // CLASSIC_DYNAMICS_H_INCLUDED
