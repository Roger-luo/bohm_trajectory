#ifndef BOHM_IMPL_H_INCLUDED
#define BOHM_IMPL_H_INCLUDED
Split::Split(double _d,double _omega,double _h,double _lambda,double _L,double _working_precision)
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
//Gamma
std::complex<double>
Split::Gamma(double D) const
{
        std::complex<double> res(cos(2*PI/lambda*(L+D)/(L*D)),sin(2*PI/lambda*(L+D)/(L*D)));
        return res;
}
//C_Z
std::complex<double>
Split::C_Z(double D) const
{
        auto f = [this,D](double z)
        {
                std::complex<double> res(cos(2*PI/lambda)*(z*z*0.5/L+z*z*0.5/D),sin(2*PI/lambda)*(z*z*0.5/L+z*z*0.5/D));
                return res;
        };

        std::complex<double> res;
        res = Intergate<std::complex<double>>(f,-h,h,working_precision);
        return res;
}

//K
//Kernel with one split open
std::complex<double>
Split::K_Single(double x,double low,double up,double D) const
{
        auto f =[x,low,this,D](double y,double z)
        {
                std::complex<double> res(std::cos(2*PI/lambda*(std::sqrt(y*y+L*L+z*z)+std::sqrt((x-y)*(x-y)+D*D+z*z))),std::sin(2*PI/lambda*(std::sqrt(y*y+L*L+z*z)+std::sqrt((x-y)*(x-y)+D*D+z*z))));
                return res;
        };

        return (-1/(lambda*lambda))*Intergate<std::complex<double>>(f,low,up,-h,h,working_precision);
}
/*

//K
//Kernel with one split open
std::complex<double> Split::K_Single(double x,double low,double up,double D) const
{
        auto f =[x,low,this,D](double y)
        {
                std::complex<double> res(std::cos(2*PI/lambda*(y*y*0.5/L+(y-x)*(y-x)*0.5/D)),std::sin(2*PI/lambda*(y*y*0.5/L+(y-x)*(y-x)*0.5/D)));
                return res;
        };

        return (-1/(lambda*lambda))*Intergate<std::complex<double>>(f,low,up,working_precision);
}
*/

//K2
    //Kernel with two split open
std::complex<double>
Split::K_double(double x,double low1,double up1,double low2,double up2,double D) const
{
        auto f=[x,low1,up1,low2,up2,this,D](double y1,double y2)
        {
            std::complex<double> res(std::sqrt(std::abs(y2-y1))*cos(2*PI/lambda*(y1*y1*0.5/L+std::abs(y2-y1)+(x-y2)*(x-y2)*0.5/D)),std::sqrt(std::abs(y2-y1))*std::sin(2*PI/lambda*(y1*y1*0.5/L+std::abs(y2-y1)+(x-y2)*(x-y2)*0.5/D)));
            return res;
        };

        std::complex<double> res;
        res = Intergate<std::complex<double>>(f,low1,up1,low2,up2,working_precision);
        std::complex<double> t(std::cos(3*PI/4),std::sin(3*PI/4));
        res = 1/(std::sqrt(std::pow(lambda,2.5)))*t*res;
        return res;
}
/*

//K2
    //Kernel with two split open
std::complex<double> Split::K_double(double x,double low1,double up1,double low2,double up2,double D) const
{
        auto f=[x,low1,up1,low2,up2,this,D](double y1,double y2)
        {
            std::complex<double> res(std::sqrt(std::abs(y2-y1))*cos(2*PI/lambda*(y1*y1*0.5/L+std::abs(y2-y1)+(x-y2)*(x-y2)*0.5/D)),std::sqrt(std::abs(y2-y1))*std::sin(2*PI/lambda*(y1*y1*0.5/L+std::abs(y2-y1)+(x-y2)*(x-y2)*0.5/D)));
            return res;
        };

        std::complex<double> res;
        res = Intergate<std::complex<double>>(f,low1,up1,low2,up2,working_precision);
        std::complex<double> t(std::cos(3*PI/4),std::sin(3*PI/4));
        res = 1/(std::sqrt(std::pow(lambda,2.5)))*t*res;
        return res;
}
*/
//K
//total Kernel
std::complex<double>
Split::K(double x,double D) const
{
        std::complex<double> res;
            res = Gamma(D)*C_Z(D)*(K_Single(x,-d-0.5*omega,-d+0.5*omega,D)+
            K_Single(x,-0.5*omega,+0.5*omega,D)+
            K_Single(x,d-0.5*omega,d+0.5*omega,D)+
            K_double(x,-d-omega*0.5,-d+omega*0.5,-omega*0.5,omega*0.5,D)+
            K_double(x,-omega*0.5,omega*0.5,-d-omega*0.5,-d+omega*0.5,D)+
            K_double(x,-omega*0.5,omega*0.5,d-omega*0.5,d+omega*0.5,D)+
            K_double(x,d-omega*0.5,d+omega*0.5,-omega*0.5,omega*0.5,D)+
            K_double(x,-d-omega*0.5,-d+omega*0.5,d-omega*0.5,d+omega*0.5,D)+
            K_double(x,d-omega*0.5,d+omega*0.5,-d-omega*0.5,-d+omega*0.5,D));
            return res;
}

std::complex<double>
Split::K_double_split(double x,double D) const
{
        std::complex<double> res;
        res = Gamma(D)*C_Z(D)*(K_Single(x,0.5*(-d-omega),0.5*(-d+omega),D)+K_Single(x,0.5*(d-omega),0.5*(d+omega),D));//+K_double(x,0.5*(-d-omega),0.5*(-d+omega),0.5*(d-omega),0.5*(d+omega),D)+K_double(x,0.5*(d-omega),0.5*(d+omega),0.5*(-d-omega),0.5*(-d+omega),D);
        return res;
}

double
Split::S_double(double x,double D) const
{
        return std::arg(K_double_split(x,D))*Dirac_const;
}
double Split::S(double x,double D) const
{
        return std::arg(K(x,D))*Dirac_const;
}

linear_algebra::vector<double>
Split::nabla_S(double x,double D)const
{
        double h = working_precision*10.0;
        linear_algebra::vector<double> res(2);
        res[0] = (S(x+h,D)-S(x-h,D))/(h*2.0);
        res[1] = (S(x,D+h)-S(x,D))/h;
        return res;
}

void
S_Field_Print(const Split &K,const std::vector<double> &range)
{
        std::ofstream S_Field("S_Field.dat");
        if(!S_Field.is_open()) throw "Bohm_E0001";//can not open data file of S Field
        for(double x=range[0];x<range[1];x = x+range[2])
        {
                for(double D=range[3];D<range[4];D = D+range[5])
                {
                        S_Field<<x<<"\t"<<D<<"\t"<<K.S(x,D)<<std::endl;
                }
        }
        S_Field.close();
}

//-----------------------------------------------------------
//Class Photon
//-----------------------------------------------------------
photon::photon()
{
        k=lambda=d=D=0;
        q = linear_algebra::vector<double> (2,0);
}

photon::photon(linear_algebra::vector<double> initial_pos,double _lambda,double _d,double _E,double _B,double _eta)
{
        lambda = _lambda;
        d            = _d;
        E            = _E;
        B            = _B;
        eta        = _eta;
        k            = 2*PI/lambda;
        q            = initial_pos;
        D            = 0;
}

linear_algebra::vector<double>
photon::next_timestep(double dt)
{
        double g_A,g_B,theta_A,theta_B,r_A,r_B;
        const double light_speed = vacuum_light_speed*1e6;//miu m/s

        D = q[0];

        theta_A = std::abs(std::atan((q[1]-eta)/q[0]));
        theta_B = std::abs(std::atan((q[1]+eta)/q[0]));

        r_A         = std::sqrt(q[0]*q[0]+(q[1]-eta)*(q[1]-eta));
        r_B         = std::sqrt(q[0]*q[0]+(q[1]+eta)*(q[1]+eta));

        g_A         = std::sin(k*(q[1]-eta)*d*0.5/D)/(k*(q[1]-eta)*d*0.5/D);
        g_B         = std::sin(k*(q[1]+eta)*d*0.5/D)/(k*(q[1]+eta)*d*0.5/D);

        double phi_phi = (E*E+B*B)*(g_A*g_A+g_B*g_B)+2*g_A*g_B*(E*E*std::cos(theta_A+theta_B)+B*B)*std::cos(k*(r_A-r_B));
        double v_x        = 2*E*B/phi_phi*(g_A*g_A*std::sin(theta_A)+g_B*g_B*std::sin(theta_B)+g_A*g_B*std::cos(k*(r_A-r_B))*(std::cos(theta_A)+std::cos(theta_B)));
        double v_y        = 2*E*B/phi_phi*(-g_A*g_A*std::sin(theta_A)+g_B*g_B*std::sin(theta_B)+g_A*g_B*std::cos(k*(r_A-r_B))*(std::sin(theta_A)-std::sin(theta_B)));

        linear_algebra::vector<double> v = {v_x,v_y};

        linear_algebra::vector<double> res =  (v*dt)+q;
        q = res;
        return res;
}//next_timestep

linear_algebra::vector<double>
photon::cur_pos()
{
        return q;
}

linear_algebra::vector<double>
photon::wave_par()
{
        linear_algebra::vector<double> res  = {lambda,k,E,B};
        return res;
}

linear_algebra::vector<double>
photon::slit_par()
{
        linear_algebra::vector<double> res = {d,D,eta};
        return res;
}

void
photon_trajectory(double lambda,double d,double eta,double dt,double STOP_DISTANCE)
{
        const int MAX_HALF_PHOTON_NUM = 8;
        double y_step = 1;

        std::ofstream bohm_photon_t("bohm_photon_t.dat");

        for(int i=1;i<MAX_HALF_PHOTON_NUM;i++)
        {
                std::cout<<i<<std::endl;
                photon up_photon({10,eta-d/2+i*y_step},lambda,d,1,1,eta);
                linear_algebra::vector<double> temp = {0.0,0.0};
                while(temp[0]<STOP_DISTANCE)
                {
                        temp = up_photon.cur_pos();
                        bohm_photon_t<<temp<<std::endl;
                        up_photon.next_timestep(dt);
                }

                photon down_photon({10,-eta+d/2-i*y_step},lambda,d,1,1,eta);
                temp = {0.0,0.0};
                while(temp[0]<STOP_DISTANCE)
                {
                        temp = down_photon.cur_pos();
                        bohm_photon_t<<temp<<std::endl;
                        down_photon.next_timestep(dt);
                }
        }
}

double
photon_tracer(photon &t_photon,const double dt,const double STOP_DISTANCE)
{
        double trace = 0;
        linear_algebra::vector<double> last_pos = {0.0,0.0};
        while(last_pos[0]<STOP_DISTANCE)
        {
                last_pos = t_photon.cur_pos();
                t_photon.next_timestep(dt);
                trace += std::sqrt(linear_algebra::dot(t_photon.cur_pos()-last_pos,t_photon.cur_pos()-last_pos));
        }
        return trace;
}

double
bohm_trac_classic_trac
(linear_algebra::vector<double> initial_pos,
double lambda,double _d,double _eta,
double MAX_DISTANCE,double error)
{
        bohm::photon p(initial_pos,lambda,_d,1,1,_eta);
        double dt = 1e-11;
        double temp1 = 1000;
        double temp2 = 0;
        while(abs(temp1-temp2)>error*lambda)
        {
                temp1 = bohm::photon_tracer(p,dt,MAX_DISTANCE)-sqrt(linear_algebra::dot(p.cur_pos()-initial_pos,p.cur_pos()-initial_pos));
                temp2 = bohm::photon_tracer(p,dt*1e-1,MAX_DISTANCE)-sqrt(linear_algebra::dot(p.cur_pos()-initial_pos,p.cur_pos()-initial_pos));
                dt*=1e-1;
        }
        std::cout<<p.cur_pos()<<std::endl;
        std::cout<<temp1<<std::endl;
        return temp1-std::sqrt(linear_algebra::dot(p.cur_pos()-initial_pos,p.cur_pos()-initial_pos));
}
#endif // BOHM_IMPL_H_INCLUDED
