#pragma once
#ifndef PRECISION_H_
#define PRECISION_H_

#include <math.h>
#include <complex.h>

/*
  Change this to modify the precision of the real numbers used
    1 -> float
    2 -> double
    3 -> long double
*/
# define PRECISION 2
# define DECIMAL_PLACES 5 // max: 17 (MinGWx64-gcc)
# define DECIMAL_PLACES_S "5" // same as before, but as a string
# define TOLERANCE 0.000000000001L // Any number n between -TOLERANCE and TOLERANCE will be considered n = 0
# define NOTATION "g" // f for normal behaviour, g for scientific notation

# if PRECISION==1
  # define REAL float
  # define REAL_STRING_FORMAT "%." DECIMAL_PLACES_S NOTATION
  # define SIN sinf
  # define COS cosf
  # define TAN tanf
  # define ASIN asinf
  # define ACOS acosf
  # define ATAN atanf
  # define ATAN2 atan2f
  # define SQRT sqrtf
  # define LOG2 log2f
  # define MOD cabsf
  # define ARG cargf
  # define RE crealf
  # define IM cimagf
  # define CONJ conjf
  # define POW powf
# elif PRECISION==2
  # define REAL double
  # define REAL_STRING_FORMAT "%." DECIMAL_PLACES_S "l" NOTATION
  # define SIN sin
  # define COS cos
  # define TAN tan
  # define ASIN asin
  # define ACOS acos
  # define ATAN atan
  # define ATAN2 atan2
  # define SQRT sqrt
  # define LOG2 log2
  # define MOD cabs
  # define ARG carg
  # define RE creal
  # define IM cimag
  # define CONJ conj
  # define POW pow
# elif PRECISION==3
  # define REAL long double
  # define REAL_STRING_FORMAT "%." DECIMAL_PLACES_S "L" NOTATION
  # define SIN sinl
  # define COS cosl
  # define TAN tanl
  # define ASIN asinl
  # define ACOS acosl
  # define ATAN atanl
  # define ATAN2 atan2l
  # define SQRT sqrtl
  # define LOG2 log2l
  # define MOD cabsl
  # define ARG cargl
  # define RE creall
  # define IM cimagl
  # define CONJ conjl
  # define POW powl
# else
  #error Impossible precision! Please choose 1, 2 or 3!
# endif
/*
  Change this to modify the type of the index
    1 -> unsigned short
    2 -> unsigned int
    3 -> unsigned long
    4 -> unsigned long long
*/
# define MAXSIZE 2

# if MAXSIZE==1
  # define INTEGER short
  # define INTEGER_STRING_FORMAT "%hd"
  # define UINTEGER unsigned short
  # define UINTEGER_STRING_FORMAT "%hu"
# elif MAXSIZE==2
  # define INTEGER int
  # define INTEGER_STRING_FORMAT "%d"
  # define UINTEGER unsigned int
  # define UINTEGER_STRING_FORMAT "%u"
# elif MAXSIZE==3
  # define INTEGER long
  # define INTEGER_STRING_FORMAT "%ld"
  # define UINTEGER unsigned long
  # define UINTEGER_STRING_FORMAT "%lu"
# elif MAXSIZE==4
  # define INTEGER long long
  # define INTEGER_STRING_FORMAT "%lld"
  # define UINTEGER unsigned long long
  # define UINTEGER_STRING_FORMAT "%llu"
# else
  #error Impossible max size! Please choose 1, 2, 3 or 4!
# endif

# define Q_PI 3.14159265358979323846264338327950288419716939937510582097494459230781641L
# define Q_TAU 6.28318530717958647692528676655900576839433879875021164194988918461563281L // Let Q_TAU = 2 * Pi (3.14159265...)
# define Q_SQRT2 1.41421356237309504880168872420969807856967187537694807317667973799073248L
# define Q_SQRT1_2 0.70710678118654752440084436210484903928483593768847403658833986899536624L

# endif
