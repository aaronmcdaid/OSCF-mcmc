#ifndef MACROS_HPP__
#define MACROS_HPP__

#include <iostream>
#define PP1(x)                 std :: cout << #x << ":" << (x) << std :: endl
#define PP2(x,y)               std :: cout << #x << ',' << #y                                           << ":\t" << (x) << " , " << (y) << std :: endl
#define PP3(x,y,z)             std :: cout << #x << ',' << #y << ',' << #z                              << ":\t" << (x) << " , " << (y) << " , " << (z) << std :: endl
#define PP4(x,y,z,w)           std :: cout << #x << ',' << #y << ',' << #z << ',' << #w                 << ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << std :: endl
#define PP5(x,y,z,w,v)         std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v    << ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << std :: endl
#define PP6(x,y,z,w,v,u)       std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v << ',' << #u      \
	<< ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << " , " << (u) << std :: endl
#define PP7(x,y,z,w,v,u,t)     std :: cout << #x << ',' << #y << ',' << #z << ',' << #w << ',' << #v << ',' << #u << ',' << #t      \
	<< ":\t" << (x) << " , " << (y) << " , " << (z) << " , " << (w) << " , " << (v) << " , " << (u) << " , " << (t) << std :: endl
#define unless(x) if(!(x))
#define DYINGWORDS(x) for (int klsdjfslkfj = (x) ? 0 : 1; klsdjfslkfj!=0; klsdjfslkfj--, ({ assert (x); }) )
#define VERYCLOSE(a,b) (1e-07 > fabs((a)-(b)))
#define VERYCLOSE2(a,b) (1e-02 > fabs((a)-(b)))
#define VERYCLOSE3(a,b) (1e-03 > fabs((a)-(b)))
#define VERYCLOSE4(a,b) (1e-04 > fabs((a)-(b)))
#define VERYCLOSE5(a,b) (1e-05 > fabs((a)-(b)))
#define For(it, container) for( typeof((container).begin()) it = (container).begin(); it != (container).end(); ++it)
#define ELAPSED() (double(clock()) / CLOCKS_PER_SEC)
#define assertVERYCLOSE(a,b) assert(VERYCLOSE(a,b))
//#define assertVERYCLOSE2(a,b) assert(VERYCLOSE2(a,b))
#define assertVERYCLOSE2(a_,b_) do{ auto a=a_; auto b=b_; unless(VERYCLOSE2(a,b)) { PP(a,b,a-b); } ;assert(VERYCLOSE2(a,b)); } while(0)
#define assertVERYCLOSE3(a_,b_) do{ auto a=a_; auto b=b_; unless(VERYCLOSE3(a,b)) { PP(a,b,a-b); } ;assert(VERYCLOSE3(a,b)); } while(0)
#define assertVERYCLOSE4(a,b) assert(VERYCLOSE4(a,b))
#define assertVERYCLOSE5(a,b) assert(VERYCLOSE5(a,b))
#define assertEQ(a,b)        assert((a)==(b))

template<typename T>
static inline void Ignore(const T&) { }


#define GET_ARG_11(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,N,...) N
#define COUNT_MACRO_ARGS(...) GET_ARG_11( __VA_ARGS__, 10,9,8,7,6,5,4,3,2,1,I_CANNOT_SEE_ZERO_ARGS)

#define SELECT_PP_IMPL( n ) PP ## n
#define SELECT_PP( n ) SELECT_PP_IMPL(n)
#define PP(...) SELECT_PP( COUNT_MACRO_ARGS(__VA_ARGS__) )(__VA_ARGS__)

#define PPL(...) do{ cout << __FILE__ ":" << __LINE__ << "\t   "; PP(__VA_ARGS__); } while(0)

#define staticcast(T,v) static_cast<T>(v)

#endif

//#include "dummy_assert.hpp"
