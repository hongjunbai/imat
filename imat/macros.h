#ifndef IMAT_MACROS_H
#define IMAT_MACROS_H
#include <iostream>

namespace imat {

#define MFOR(i, n) for (long M_TMP = (n), (i) = 0; (i) < M_TMP; ++(i))

#define MRANGE(i, a, b) for ( long M_TMP = (b), (i) = (a); (i) < M_TMP; ++(i))
#define MRANGE_I(i, a, b) for (int M_TMP = (b), (i) = (a); (i) < M_TMP; ++(i))

#define MLOOP_CONST(t, i, c) for(t::const_iterator (i) = (c).begin(); (i) != (c).end(); ++(i))
#define MLOOP(t, i, c)       for(t::iterator       (i) = (c).begin(); (i) != (c).end(); ++(i))

// Used in template functions, since typename is required for the compiler to parse the
// template functions
#define MLOOP_CONST_T(t, i, c) for(typename t::const_iterator (i) = (c).begin(); (i) != (c).end(); ++(i))
#define MLOOP_T(t, i, c) for(typename t::const_iterator (i) = (c).begin(); (i) != (c).end(); ++(i))

#define SMALLER(a, b) (a) < (b) ? (a) : (b)
#define LARGER(a, b) (a) > (b) ? (a) : (b)

#define MSHOW(x)       do { std::cout << #x << " = " << (x) << std::endl; } while(false)
#define MESHOW(x)      do { std::cerr << (x) << std::endl; } while(false)
#define MSSHOW(x)      do { std::cout << (x) << std::endl; } while(false)

}

#endif // IMAT_MACROS_H
