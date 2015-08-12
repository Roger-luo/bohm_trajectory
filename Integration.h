#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED
#include<iostream>
#include<cmath>
#include<random>
#include<functional>
#include<vector>
//------------------------------------------------------------------------------------------------------------------------------------------
//Romberg Intergration
//------------------------------------------------------------------------------------------------------------------------------------------
template<typename Data ,int INITIAL_INDEX = 64>
class Triangle
{
public:
        Data **dp;
        int data_size;

        Triangle(Data INITIAL_DATA)
        {
                dp = new Data *[INITIAL_INDEX];
                for(int i=0;i<INITIAL_INDEX;i++)
                {
                        dp[i] = new Data [INITIAL_INDEX-i];
                        for(int j=0;j<INITIAL_INDEX-i;j++)
                        {
                                dp[i][j]  = INITIAL_DATA;
                        }
                }
                data_size = INITIAL_INDEX;
        }

        Triangle(Triangle &_T)
        {
                dp = new Data *[_T.data_size];
                for(int i=0;i<_T.data_size;i++)
                {
                        dp[i] = new Data[_T.data_size-i];
                        for(int j=0;j<_T.data_size-i;j++)
                        {
                                dp[i][j] = _T.dp[i][j];
                        }
                }
        }

        void Enlarge_Size(int new_size)
        {
                if(new_size<data_size) return;
                Data **temp_dp;
                temp_dp = new Data* [new_size];
                for(int i=0;i<new_size;i++)
                {
                        temp_dp[i] = new Data [new_size-i];
                }

                for(int i=0;i<data_size;i++)
                        for(int j=0;j<data_size-i;j++)
                                temp_dp[i][j] = dp[i][j];

                for(int i=0;i<data_size;i++)
                {
                        delete []dp[i];
                }
                delete []dp;

                data_size = new_size;
                dp             = temp_dp;
        }

        ~Triangle()
        {
                for(int i=0;i<data_size;i++)
                {
                        delete []dp[i];
                }
                delete []dp;
        }
};


template<typename TYPE ,int MAX_STEP = 20, typename F = TYPE (*)(double) >
TYPE Intergate(F const &func,double low,double up,double working_precision)
{
        double h = up - low;
        Triangle<TYPE,32> R(0);

        R.dp[0][0] = (0.5*h*(func(up)+func(low)));

        bool flag = 1;
        int cur = 1;
        while(flag)
        {
                if(cur>MAX_STEP) throw 10000;//"10000" is out of range
                else if(cur>R.data_size) R.Enlarge_Size(R.data_size<<1);
                //sum
                TYPE sum = 0;
                for(int i=1;i<(1<<(cur-1))+1;i++)
                {
                        sum += func(low+(2*i-1)/(double)(1<<cur)*h);
                }

                R.dp[0][cur] = 0.5*(R.dp[0][cur-1] + h/(double)((1<<cur)-1)*sum);

                for(int j=1;j<cur+1;j++)
                {
                        R.dp[j][cur-j] = R.dp[j-1][cur-j+1] + 1/(double)((1<<(j<<1))-1)*
                                                   (R.dp[j-1][cur-j+1]-R.dp[j-1][cur-j]);
                }

                if(abs(R.dp[cur][0]-R.dp[cur-1][0]) < working_precision) flag = 0;
                else cur++;
        }
        return R.dp[cur][0];
}

template<typename TYPE ,int MAX_STEP = 100, typename F >
TYPE Intergate(const F &func,double y,double low,double up,double working_precision)
{
        double h = up - low;
        Triangle<TYPE,32> R(0);

        R.dp[0][0] = (0.5*h*(func(up,y)+func(low,y)));

        bool flag = 1;
        int cur = 1;
        while(flag)
        {
                if(cur>MAX_STEP) throw 10000;//"10000" is out of range
                else if(cur>R.data_size) R.Enlarge_Size(R.data_size<<1);
                //sum
                TYPE sum = 0;
                for(int i=1;i<(1<<(cur-1))+1;i++)
                {
                        sum += func(low+(2*i-1)/(double)(1<<cur)*h,y);
                }

                R.dp[0][cur] = 0.5*(R.dp[0][cur-1] + h/(double)((1<<cur)-1)*sum);

                for(int j=1;j<cur+1;j++)
                {
                        R.dp[j][cur-j] = R.dp[j-1][cur-j+1] + 1/(double)((1<<(j<<1))-1)*
                                                   (R.dp[j-1][cur-j+1]-R.dp[j-1][cur-j]);
                }

                if(abs(R.dp[cur][0]-R.dp[cur-1][0]) < working_precision) flag = 0;
                else cur++;
        }
        return R.dp[cur][0];
}

template<typename TYPE,int MAX_STEP = 20,typename F>
TYPE Intergate(const F &f,double low1,double up1,double low2,double up2,double working_precision)
{
       std:: function<TYPE(double)> First_Kernel=[low1,up1,low2,up2,working_precision,f](double y)
        {
                TYPE res = Intergate<TYPE>(f,y,low1,up1,working_precision);
                return res;
        };

        TYPE res = Intergate<TYPE>(First_Kernel,low2,up2,working_precision);
        return res;
}

//Monte Carlo
//f has to be TYPE (vector )
template<typename TYPE,typename F,int MAX_STEP = 100>
TYPE Intergate(const F &f,std::vector<double> low,std::vector<double>up,double working_precision)
{
        double Volume = 1;
        if(low.size()!=up.size()) throw 10001;//Error
        for(unsigned int i=0;i<low.size();i++)
        {
                Volume *= up[i]-low[i];
        }

        //random vector
        std::vector<double> sample;
        std::random_device rd;
        std::mt19937 gen(rd());


        bool flag = 1; int N = 1;
        TYPE sum = 0;TYPE mean = 0; double Var = 0;
        std::vector<TYPE> sample_value;
        while(flag)
        {
                sample.clear();
                for(unsigned int i=0;i<low.size();i++)
                {
                        std::uniform_real_distribution<double> dis(low[i],up[i]);
                        sample.push_back(dis(gen));
                }

                sample_value.push_back(f(sample));
                sum += f(sample);
                mean = sum/N;

                double temp_sum = 0;
                for(unsigned int i=0;i<sample_value.size();i++)
                {
                        temp_sum +=std::abs((sample_value[i]-mean))*std::abs((sample_value[i]-mean));
                }
                Var = temp_sum/(sample_value.size()-1);


                double Error_Bar = Volume*sqrt(Var/N);


                if(Error_Bar<working_precision) flag = 0;
                else N++;
        }

        TYPE res = Volume*mean;
        return res;
}
#endif // INTEGRATION_H_INCLUDED
