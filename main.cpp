#include<iostream>
#include<functional>
#include<cstdlib>
#include"mechanics/mechanics.h"
using namespace std;

int main()
{
/*
        linear_algebra::vector<double> initial_pos = {1e4,2300};
        //cout<<bohm::bohm_trac_classic_trac(initial_pos,0.943,80,2345,1e7,0.1);

        bohm::photon_trajectory(0.5,20,200,10,1e3);
        system("gnuplot plot.gnu");
        system("eog bohm_photo.jpeg");

       // bohm::photon test(initial_pos,0.943,80,1,1,2345);
        //cout<<bohm::photon_tracer(test,1e-11,1e7);
*/
        bohm::Split A(0.272,0.062,4,5e-5,30.5e4,1e-3);
        const vector<double> range = {-100,100,1,12e4,24e4,1e3};

        ofstream file("S_Field.dat");
        for(double x=1e5;x<24e5;x=x+1e4)
        {
                file<<x<<"\t"<<A.S_double(0,x)<<endl;
        }

        system("gnuplot plot.gnu");
        system("eog S_Field.jpeg");

       // bohm::S_Field_Print(A,range);
/*
        linear_algebra::vector<double> initial_pos = {0,1e4,0.0};

        linear_algebra::vector<double> initial_vel = {0,0,0};
        cout<<initial_pos<<endl;
        classic_dynamic p(initial_pos,initial_vel,Electron_Mass);

        ofstream file("S_D.dat");
        for(double D=12e4;D<24e4;D = D+1e3)
        {
                file<<D<<"\t"<<A.S(0,D)<<endl;
        }

        function<linear_algebra::vector<double>(linear_algebra::vector<double>)> speed = [A,p](linear_algebra::vector<double> pos)
        {d
                return A.nabla_S(pos[0],pos[1])/p.static_m();
        };

        for(int i=0;i<1000;i++)
        {
                cout<<"Position:"<<p.cur_pos()<<endl;
                cout<<"Speed:"<<p.cur_vel()<<endl;
                linear_algebra::vector<double> cur_speed = speed(p.cur_pos());
                p.next_timestep(cur_speed,1e-2);
        }
        */
        return 0;
}
