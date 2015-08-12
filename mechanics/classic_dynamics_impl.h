#ifndef CLASSIC_DYNAMICS_IMPL_H_INCLUDED
#define CLASSIC_DYNAMICS_IMPL_H_INCLUDED
particle::particle():
static_mass(0),cur_mass(0),charge(0),spin(0){}

particle::particle(const linear_algebra::vector<double> &initial_pos,const double &_static_m)
{
        position = initial_pos;
        velocity  = 0;
        static_mass = _static_m;
        cur_mass     = _static_m;
        charge          = 0;
        spin               = 0;
}

particle::particle(const linear_algebra::vector<double> &initial_pos,const linear_algebra::vector<double> &initial_vel,const double &_static_m)
{
        position = initial_pos;
        velocity  = initial_vel;
        static_mass = _static_m;
        cur_mass     = _static_m/std::sqrt(1-linear_algebra::dot(initial_vel,initial_vel));
        charge          = 0;
        spin               = 0;
}

particle::particle(const linear_algebra::vector<double> &initial_pos,const linear_algebra::vector<double> &initial_vel,const double &_static_m,const double &_charge)
{
        position = initial_pos;
        velocity  = initial_vel;
        static_mass = _static_m;
        cur_mass     = _static_m/std::sqrt(1-linear_algebra::dot(initial_vel,initial_vel));
        charge          = _charge;
        spin               = 0;
}

particle::particle(const linear_algebra::vector<double> &initial_pos,const linear_algebra::vector<double> &initial_vel,const double &_static_m,const double &_charge,const double &_spin)
{
        position = initial_pos;
        velocity  = initial_vel;
        static_mass = _static_m;
        cur_mass     = _static_m/std::sqrt(1-linear_algebra::dot(initial_vel,initial_vel));
        charge          = _charge;
        spin               = _spin;
}

particle::~particle()
{}

inline linear_algebra::vector<double> particle::cur_pos()
{
        return position;
}

inline linear_algebra::vector<double> particle::cur_vel()
{
        return velocity;
}

double particle::cur_m() const
{
        return cur_mass;
}

double particle::static_m() const
{
        return static_mass;
}
void particle::reset_pos(const linear_algebra::vector<double> &_cur)
{
        position = _cur;
}

void particle::reset_vel(const linear_algebra::vector<double> &_cur)
{
        velocity = _cur;
}
/*
template<typename Force_Type>
linear_algebra::vector<double> classic_dynamic::next_timestep(const std::function<Force_Type(linear_algebra::vector<double> &)> &f,const double &dt)
{
        Force_Type F = f(position);

}
*/

void classic_dynamic::next_timestep(const linear_algebra::vector<double> vel,const double dt)
{
        //linear_algebra::vector<double> cur_position(cur_pos()+(dt*vel));
        reset_pos(cur_pos()+(dt*vel));
        reset_vel(vel);
        return ;
}
#endif // CLASSIC_DYNAMICS_IMPL_H_INCLUDED
