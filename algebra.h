#ifndef ALGEBRA_H_INCLUDED
#define ALGEBRA_H_INCLUDED
#include<vector>
#include<iostream>
#include<cmath>
#include<complex>
#include <cassert>

namespace linear_algebra
{
template<typename T>
class vector
{
public:
        T *element;
        int size;
        typedef                 T* iterator;
        typedef const     T* const_iterator;
        //constructor & destrucstor
        vector();
        vector(T a,T b, T c);
        vector(int len,const T &x = T(0));
        vector(T *_element,int length);
        vector(const vector<T> &V);
        vector(std::initializer_list<T> _IL);
        ~vector();

        //assignments
        vector<T> &operator = (const vector<T> &V);
        vector<T> &operator = (const T &x);

        //accessors
        T& operator[](int i);
        const T& operator[](int i) const;

        //iterators
        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;

        //type conversion
        operator T*();
        operator const T*() const;

        //others
        vector<T>& resize( int length );
        //compute assignment
        vector<T>& operator+=(const T&);
        vector<T>& operator-=(const T&);
        vector<T>& operator*=(const T&);
        vector<T>& operator/=(const T&);
        vector<T>& operator+=(const vector<T>&);
        vector<T>& operator-=(const vector<T>&);
        vector<T>& operator*=(const vector<T>&);
        vector<T>& operator/=(const vector<T>&);
private:
        void CopyFromArray(const T *_element,const int len);
        void destroy()
        {
                if(element == NULL) return;
                delete []element;
        }
};

        // input and output
        template<typename T>
        std::ostream& operator<<( std::ostream&, const vector<T>& );
        template<typename T>
        std::istream& operator>>( std::istream&, vector<T>& );

        // arithmetic operators
        template<typename T>
        vector<T> operator-( const vector<T>& );
        template<typename T>
        vector<T> operator+( const vector<T>&, const T& );
        template<typename T>
        vector<T> operator+( const T&, const vector<T>& );
        template<typename T>
        vector<T> operator+( const vector<T>&, const vector<T>& );
        template<typename T>
        vector<T> operator-( const vector<T>&, const T& );
        template<typename T>
        vector<T> operator-( const T&, const vector<T>& );
        template<typename T>
        vector<T> operator-( const vector<T>&, const vector<T>& );
        template<typename T>
        vector<T> operator*( const vector<T>&, const T& );
        template<typename T>
        vector<T> operator*( const T&, const vector<T>& );
        template<typename T>
        vector<T> operator*( const vector<T>&, const vector<T>& );
        template<typename T>
        vector<T> operator/( const vector<T>&, const T& );
        template<typename T>
        vector<T> operator/( const T&, const vector<T>& );
        template<typename T>
        vector<T> operator/( const vector<T>&, const vector<T>& );

        //dot & cross
        template<typename T>
        T dot(const vector<T>&,const vector<T>&);
        template<typename T>
        vector<T> cross(const vector<T>&,const vector<T>&);

        //utilities
        template<typename T>
        T sum(const vector<T>&);
        template<typename T>
        T mean(const vector<T> &);
        template<typename T>
        vector<T> abs(const vector<T> &);
        template<typename Type>
        vector<Type> abs(const vector<std::complex<Type>> &v);
        template<typename T>
        vector<T> real(const vector<std::complex<T>> &);
        template<typename T>
        vector<T> imag(const vector<std::complex<T>> &);
        #include"algebra/algebra_vec_impl.h"
}


#endif // ALGEBRA_H_INCLUDED
