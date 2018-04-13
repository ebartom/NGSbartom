#ifndef _PEAR_
#define _PEAR_
#include "config.h"

#define PROGRAM_NAME            "PEAR"
#define PROGRAM_VERSION         "0.9.10"
#define VERSION_DATE            "May 30, 2016"
#define LICENCE                 "Creative Commons Licence"
#define CONTACT                 "Tomas.Flouri@h-its.org and Jiajie.Zhang@h-its.org"

#if (defined(HAVE_BZLIB_H) && defined(HAVE_ZLIB_H))
#define COMPILE_INFO            " - [+bzlib +zlib]"
#elif (defined(HAVE_BZLIB_H) && !defined(HAVE_ZLIB_H))
#define COMPILE_INFO            " - [+bzlib]"
#elif (!defined(HAVE_BZLIB_H) && defined(HAVE_ZLIB_H))
#define COMPILE_INFO            " - [+zlib]"
#else
#define COMPILE_INFO           ""
#endif

#endif
