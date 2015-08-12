#ifndef ALGEBRA_VEC_IMPL_H_INCLUDED
#define ALGEBRA_VEC_IMPL_H_INCLUDED



template <typename T>
inline void vector<T>::CopyFromArray(const T *_element,const int len)
{
        if(_element == NULL) throw "E00001";//no value to copy
        else
        {
                size = len;
                element = new T[size];
                for(int i=0;i<size;i++)
                {
                        element[i] = _element[i];
                }
        }
}

template<typename T>
vector<T>::vector()
{
        element = nullptr;
        size = 0;
}

template<typename T>
vector<T>::vector(const T a,const T b,const T c)
{
        element = new T[3];
        element[0] = a;
        element[1] = b;
        element[2] = c;
        size = 3;
}

template<typename T>
vector<T>::vector(int len,const T&x)
{
        size = len;
        element = new T[len];
        for(int i=0;i<len;i++)
                element[i] = x;
}

template<typename T>
vector<T>::vector(T *_element,int length)
{
        element = new T[length];
        CopyFromArray(_element,length);
}

template<typename T>
vector<T>::vector(const vector<T> &V)
{
        element = new T[V.size];
        size        = V.size;
        for(int i=0;i<size;i++)
        {
                element[i] = V.element[i];
        }
}

template<typename T>
vector<T>::vector(std::initializer_list<T> _IL)
{
        element = new T[_IL.size()];
        double *temp = element;
        size        = _IL.size();
        typename std::initializer_list<T>::iterator it;
        for(it = _IL.begin();it!=_IL.end();++it)
        {
                *(temp++) = *it;
        }
}

template<typename T>
vector<T>::~vector()
{
        destroy();
}

template<typename T>
vector<T>& vector<T>::operator = (const vector<T> &V)
{
        if(element == V.element) return *this;
        if(size == V.size)
        {
                for(int i=0;i<size;i++)
                        element[i] = V.element[i];
        }
        else
        {
                destroy();
                element = new T[V.size];
                size = V.size;
                CopyFromArray(V.element,V.size);
        }
        return *this;
}

template<typename T>
vector<T>& vector<T>::operator = (const T &x)
{
        for(int i=0;i<size;i++)
                element[i] = x;
        return *this;
}

template<typename T>
inline T& vector<T>::operator[](int i)
{
#ifdef BOUNDS_CHECK
        assert(0<=i);
        assert(i<size);
#endif // BOUNDS_CHECK
        return element[i];
}

template<typename T>
inline const T& vector<T>::operator[](int i) const
{
#ifdef BOUNDS_CHECK
        assert(0<=i);
        assert(i<size);
#endif // BOUNDS_CHECK
        return element[i];
}

template<typename T>
inline typename vector<T>::iterator vector<T>::begin()
{
        return element;
}

template<typename T>
inline typename vector<T>::const_iterator vector<T>::begin() const
{
        return element;
}

template<typename T>
inline typename vector<T>::iterator vector<T>::end()
{
        return element + size;
}

template<typename T>
inline typename vector<T>::const_iterator vector<T>::end() const
{
        return element + size;
}

template <typename Type>
vector<Type>& vector<Type>::resize( int length )
{
	if( size == length )
		return *this;

        Type *temp = new Type[length];
        if(size > length)
        {
                for(int i=0;i<length;i++)
                        temp[i] = element[i];
        }
        else if(size < length)
        {
                for(int i=0;i<size;i++)
                        temp[i] = element[i];
                for(int i=size;i<length;i++)
                        temp[i] = 0;
        }
	destroy();
	element = temp;
	size = length;

	return *this;
}

/**
 * compound assignment operators +=
 */
template <typename T>
vector<T>& vector<T>::operator+=( const T &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ += x;

	return *this;
}

template <typename T>
vector<T>& vector<T>::operator+=( const vector<T> &rhs )
{
    assert( size == rhs.size );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ += *itrR++;

	return *this;
}


/**
 * compound assignment operators -=
 */
template <typename T>
vector<T>& vector<T>::operator-=( const T &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ -= x;

	return *this;
}

template <typename T>
vector<T>& vector<T>::operator-=( const vector<T> &rhs )
{
    assert( size == rhs.size );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ -= *itrR++;

	return *this;
}


/**
 * compound assignment operators *=
 */
template <typename T>
vector<T>& vector<T>::operator*=( const T &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ *= x;

	return *this;
}

template <typename T>
vector<T>& vector<T>::operator*=( const vector<T> &rhs )
{
    assert( size == rhs.size );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ *= *itrR++;

	return *this;
}


/**
 * compound assignment operators /=
 */
template <typename T>
vector<T>& vector<T>::operator/=( const T &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ /= x;

	return *this;
}

template <typename Type>
vector<Type>& vector<Type>::operator/=( const vector<Type> &rhs )
{
    assert( size == rhs.size );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ /= *itrR++;

	return *this;
}


/**
 * Overload the output stream function.
 */
template <typename Type>
std::ostream& operator<<( std::ostream &out, const vector<Type> &v )
{
	int N = v.size;
	for( int i=0; i<N; ++i )
		out << v[i] << "\t";

        out<<std::endl;
	return out;
}


/**
 * Overload the input stream function.
 */

/**
 * get negative vector
 */
template <typename Type>
vector<Type> operator-( const vector<Type> &v )
{
        vector<Type> tmp(v.size);
        typename vector<Type>::iterator itrL = tmp.begin();
        typename vector<Type>::const_iterator itrR = v.begin();

        while( itrL != tmp.end() )
                *itrL++ = -(*itrR++);

        return tmp;
}


/**
 * vector-scalar addition.
 */
template <typename Type>
inline vector<Type> operator+( const vector<Type> &v, const Type &x )
{
	vector<Type> tmp( v );
	return tmp += x;
}

template <typename Type>
inline vector<Type> operator+( const Type &x, const vector<Type> &v )
{
	return v+x;
}


/**
 * vector-scalar substraction.
 */
template <typename Type>
inline vector<Type> operator-( const vector<Type> &v, const Type &x )
{
	vector<Type> tmp( v );
	return tmp -= x;
}

template <typename Type>
inline vector<Type> operator-( const Type &x, const vector<Type> &v )
{
	vector<Type> tmp( v );
	return -tmp += x;
}


/**
 * vector-scalar multiplication.
 */
 template <typename Type>
inline vector<Type> operator*( const vector<Type> &v, const Type &x )
{
	vector<Type> tmp = v;
	for(int i=0;i<tmp.size;i++)
	{
                tmp[i] = tmp[i]*x;
        }
	return tmp ;
}

template <typename Type,typename T>
inline vector<Type> operator*( const vector<Type> &v, const T &x )
{
	vector<Type> tmp = v;
	for(int i=0;i<tmp.size;i++)
	{
                tmp[i] = tmp[i]*x;
        }
	return tmp ;
}

template <typename Type>
inline vector<Type> operator*( const Type &x, const vector<Type> &v )
{
	vector<Type> tmp( v );
	return tmp *= x;
}

template <typename T,typename Type>
inline vector<Type> operator*( const T &x, const vector<Type> &v )
{
	vector<Type> tmp( v );
	return tmp *= x;
}


/**
 * vector-scalar division.
 */
 template <typename Type>
inline vector<Type> operator/(const vector<Type> &v, const Type &x )
{
	vector<Type> tmp( v );
	return tmp /= x;
}

template <typename Type,typename T>
inline vector<Type> operator/(const vector<Type> &v, const T &x )
{
	vector<Type> tmp( v );
	return tmp /= x;
}

template <typename Type>
inline vector<Type> operator/( const Type &x, const vector<Type> &v )
{
	int N = v.size;
	vector<Type> tmp( N );

	for( int i=0; i<N; ++i )
		tmp[i] = x / v[i];

	return tmp;
}

template <typename Type,typename T>
inline vector<Type> operator/( const T &x, const vector<Type> &v )
{
	int N = v.size;
	vector<Type> tmp( N );

	for( int i=0; i<N; ++i )
		tmp[i] = x / v[i];

	return tmp;
}


/**
 * vector-vector addition.
 */
template <typename Type>
inline vector<Type> operator+( const vector<Type> &v1, const vector<Type> &v2 )
{
    vector<Type> tmp( v1 );
	return tmp += v2;
}


/**
 * vector-vector substraction.
 */
template <typename Type>
inline vector<Type> operator-( const vector<Type> &v1, const vector<Type> &v2 )
{
   vector<Type> tmp( v1 );
	return tmp -= v2;
}


/**
 * vector-vector multiplication.
 */
template <typename Type>
inline vector<Type> operator*( const vector<Type> &v1, const vector<Type> &v2 )
{
    vector<Type> tmp( v1 );
	return tmp *= v2;
}

template<typename T>
T dot(const vector<T>& v1,const vector<T>& v2)
{
        T res = 0;
        if(v1.size!=v2.size) throw "V00002";//Vector does not match
        for(int i=0;i<v1.size;i++)
        {
                res += v1[i]*v2[i];
        }
        return res;
}

template<typename T>
vector<T> cross(const vector<T>& v1,const vector<T>& v2)
{
        if(v1.size>3||v2.size>3) throw "V00003";//do not have cross product

        vector<T> res(3);
        vector<T> tmp_v1(v1);
        vector<T> tmp_v2(v2);

        if(tmp_v1.size<3)
        {
               tmp_v1.resize(3);
        }
        else if(tmp_v2.size<3)
        {
                tmp_v2.resize(3);
        }
        res[0] = tmp_v1[1]*tmp_v2[2]-tmp_v1[2]*tmp_v2[1];
        res[1] = tmp_v1[2]*tmp_v2[0]-tmp_v1[0]*tmp_v2[2];
        res[2] = tmp_v1[0]*tmp_v2[1]-tmp_v1[1]*tmp_v2[0];
        return res;
}

template<typename Type>
Type sum(const vector<Type>& v)
{
        Type res = 0;
        for(int i=0;i<v.size;i++)
        {
                res+= v[i];
        }
        return res;
}

template<typename Type>
Type mean(const vector<Type> &v)
{
        Type res = sum(v)/v.size;
        return res;
}

template<typename Type>
Type abs(const vector<Type> &v)
{
        return dot(v,v);
}

template<typename Type>
vector<Type> abs(const vector<std::complex<Type>> &v)
{
        vector<Type> tmp(v.size);
        typename vector<Type>::iterator itrL = tmp.begin();
        typename vector< std::complex<Type> >::const_iterator itrR = v.begin();

        while( itrL != tmp.end() )
                *itrL++ = std::abs(*itrR++);

        return tmp;
}

template<typename Type>
vector<Type> real(const vector<std::complex<Type>> &v)
{
        vector<Type> tmp( v.size );
        typename vector<Type>::iterator itrL = tmp.begin();
        typename vector<std::complex<Type> >::const_iterator itrR = v.begin();

        while( itrL != tmp.end() )
                *itrL++ = (*itrR++).real();

        return tmp;
}

template<typename Type>
vector<Type> imag(const vector<std::complex<Type>> &v)
{
        vector<Type> tmp( v.dim() );
        typename vector<Type>::iterator itrL = tmp.begin();
        typename vector<std::complex<Type> >::const_iterator itrR = v.begin();

        while( itrL != tmp.end() )
                *itrL++ = (*itrR++).imag();

        return tmp;
}
#endif // ALGEBRA_VEC_IMPL_H_INCLUDED
