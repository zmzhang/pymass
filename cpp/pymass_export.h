
#ifndef PYMASS_EXPORT_H
#define PYMASS_EXPORT_H

#ifdef PYMASS_STATIC_DEFINE
#  define PYMASS_EXPORT
#  define PYMASS_NO_EXPORT
#else
#  ifndef PYMASS_EXPORT
#    ifdef pymass_EXPORTS
        /* We are building this library */
#      ifdef _MSC_VER
#      		define PYMASS_EXPORT __declspec(dllexport)
#      else
#			define PYMASS_EXPORT
#      endif
#    else
        /* We are using this library */
#      ifdef _MSC_VER
#      define PYMASS_EXPORT __declspec(dllimport)
#      else
#			define PYMASS_EXPORT
#      endif
#    endif
#  endif

#  ifndef PYMASS_NO_EXPORT
#    define PYMASS_NO_EXPORT 
#  endif
#endif

#ifndef PYMASS_DEPRECATED
#  define PYMASS_DEPRECATED __declspec(deprecated)
#endif

#ifndef PYMASS_DEPRECATED_EXPORT
#  define PYMASS_DEPRECATED_EXPORT PYMASS_EXPORT PYMASS_DEPRECATED
#endif

#ifndef PYMASS_DEPRECATED_NO_EXPORT
#  define PYMASS_DEPRECATED_NO_EXPORT PYMASS_NO_EXPORT PYMASS_DEPRECATED
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define PYMASS_NO_DEPRECATED
#endif

#endif
