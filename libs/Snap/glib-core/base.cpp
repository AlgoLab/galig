#include "stdafx.h"

#include "base.h"

void BaseTralala(){
  printf("Active defines:\n");
  #ifdef GLib_WIN
  printf("  GLib_WIN\n");
  #endif
  #ifdef GLib_WIN32
  printf("  GLib_WIN32\n");
  #endif
  #ifdef GLib_WIN64
  printf("  GLib_WIN64\n");
  #endif
  #ifdef GLib_UNIX
  printf("  GLib_UNIX\n");
  #endif
  #ifdef GLib_LINUX
  printf("  GLib_LINUX\n");
  #endif
  #ifdef GLib_SOLARIS
  printf("  GLib_SOLARIS\n");
  #endif
  #ifdef GLib_MSC
  printf("  GLib_MSC\n");
  #endif
  #ifdef GLib_CYGWIN
  printf("  GLib_CYGWIN\n");
  #endif
  #ifdef GLib_BCB
  printf("  GLib_BCB\n");
  #endif
  #ifdef GLib_GCC
  printf("  GLib_GCC\n");
  #endif
  #ifdef GLib_MACOSX
  printf("  GLib_MACOSX\n");
  #endif
  #ifdef GLib_64Bit
  printf("  GLib_64Bit\n");
  #endif
  #ifdef GLib_32Bit
  printf("  GLib_32Bit\n");
  #endif
  #ifdef GLib_GLIBC
  printf("  GLib_GLIBC\n");
  #endif
  #ifdef GLib_POSIX_1j
  printf("  GLib_POSIX_1j\n");
  #endif
}

#if defined(GLib_UNIX) && ! defined(GLib_CYGWIN)
int _daylight = 0;
#endif

#if defined(WIN32_LEAN_AND_MEAN)
int gettimeofday(struct timeval * tp, struct timezone * tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime( &system_time );
    SystemTimeToFileTime( &system_time, &file_time );
    time =  ((uint64_t)file_time.dwLowDateTime )      ;
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec  = (long) ((time - EPOCH) / 10000000L);
    tp->tv_usec = (long) (system_time.wMilliseconds * 1000);
    return 0;
}
#endif

#include "bd.cpp"
#include "fl.cpp"
#include "dt.cpp"
#include "ut.cpp"
#include "hash.cpp"

#include "unicode.cpp"
#include "unicodestring.cpp"
#include "tm.cpp"
#include "os.cpp"

#include "bits.cpp"
#include "env.cpp"
#include "wch.cpp"
#include "xfl.cpp"
#include "xmath.cpp"

#include "blobbs.cpp"
#include "lx.cpp"
#include "url.cpp"

#include "http.cpp"
#include "html.cpp"
#include "md5.cpp"
#include "ss.cpp"
#include "ssmp.cpp"
#include "xml.cpp"
#include "json.cpp"
//#include "prolog.cpp"

#include "zipfl.cpp"

