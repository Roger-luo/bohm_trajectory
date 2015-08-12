#ifndef BASIC_VEC_H_INCLUDED
#define BASIC_VEC_H_INCLUDED

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <complex>

using namespace std;

/*---------------------------------------------------------------------------------------------------------------------*/
//This is a basic vector operation file.
//Vector Operations:
//addition: vec + vec = vec
//dot product:vec & vec = double
//cross product:vec * vec = vec
//scalar multiplication:double * vec = vec
/*---------------------------------------------------------------------------------------------------------------------*/

template<typename element_type>
class vec
{
public:
    element_type *x;
    vec(element_type _a,element_type _b,element_type _c)
    {
        x = new element_type[3];
        x[0] = _a;x[1] = _b;x[2] = _c;
    }
    vec(element_type _a,element_type _b)
    {
        x = new element_type[3];
        x[0] = _a;x[1] = _b;x[2] = 0;
    }
    vec(element_type _a)
    {
        x = new element_type[3];
        x[0] = _a;x[1] = 0;x[2] = 0;
    }

    vec()
    {
        x = new element_type[3];
        x[0]=x[1]=x[2]=0;
    }

    vec(const vec &_vec)
    {
        x = new element_type[3];
        x[0] = _vec.x[0];x[1] = _vec.x[1];x[2] = _vec.x[2];
    }

    friend vec operator + (vec a,vec b)
    {
        vec res(a.x[0]+b.x[0],a.x[1]+b.x[1],a.x[2]+b.x[2]);
        return res;
    }

    friend vec operator - (vec a,vec b)
    {
        vec res(a.x[0]-b.x[0],a.x[1]-b.x[1],a.x[2]-b.x[2]);
        return res;
    }


    friend vec operator * (vec a,element_type _m)
    {
        vec res(_m*a.x[0],_m*a.x[1],_m*a.x[2]);
        return res;
    }

    friend vec operator * (element_type _m,vec a)
    {
        vec res(_m*a.x[0],_m*a.x[1],_m*a.x[2]);
        return res;
    }

    friend vec operator * (vec a,vec b)
    {
        vec res(a.x[1]*b.x[2]-a.x[2]*b.x[1],a.x[2]*b.x[0]-a.x[0]*b.x[2],a.x[0]*b.x[1]-a.x[1]*b.x[0]);
        return res;
    }

    friend double operator & (vec a,vec b)
    {
        return a.x[0]*b.x[0]+a.x[1]*b.x[1]+a.x[2]*b.x[2];
    }

    friend bool operator == (vec a,vec b)
    {
        if((a.x[0]==b.x[0])&&(a.x[1]==b.x[1])&&(a.x[2]==b.x[2]))
            return true;
        else return false ;
    }

    double abs()
    {
        double res;
        res = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
        return res;
    }


    vec &operator = (const vec &_vec)
    {
        if(this==&_vec) return *this;
        else
        {
            delete []x;
            x = new element_type[3];
            x[0] = _vec.x[0];
            x[1] = _vec.x[1];
            x[2] = _vec.x[2];
        }
        return *this;
    }

    friend ostream &operator << (ostream &out_stream,const vec &_vec)
    {
        out_stream<<_vec.x[0]<<"\t"<<_vec.x[1];//<<"\t"<<_vec.x[2];//the 3rd parameter is no use in this program,shall not be output
        return out_stream;
    }

    void display()
    {
        cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
    }



    ~vec()
    {
        delete []x;
    }
};

template<typename T>
double abs(vec<T> _vec)
{
    double res;
    res = abs(_vec.x[0])*abs(_vec.x[0])+abs(_vec.x[1])*abs(_vec.x[1])+abs(_vec.x[2])*abs(_vec.x[2]);
    return res;
}

template<typename element_type>
vec<int> floor(vec<element_type> _vec)
{
    vec<int> res(floor(_vec.x[0]),floor(_vec.x[1]),floor(_vec.x[2]));
    return res;
}

template<typename element_type>
vec<int> ceil(vec<element_type> _vec)
{
    vec<int> res(ceil(_vec.x[0]),ceil(_vec.x[1]),ceil(_vec.x[2]));
    return res;
}

template<typename element_type>
vec<element_type> imag(vec<complex<element_type> > complex_vec)
{
    vec<element_type> res(imag(complex_vec.x[0]),imag(complex_vec.x[1]),imag(complex_vec.x[2]));
    return res;
}
#endif // BASIC_VEC_H_INCLUDED
