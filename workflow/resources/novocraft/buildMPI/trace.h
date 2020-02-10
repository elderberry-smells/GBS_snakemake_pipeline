#undef tracel
#undef tracef
#include <time.h>
#if defined(TRACE) // && !defined(tracel)
    #include <stdio.h>
    #define tracef(format, ...) {\
        fprintf(stderr,"%s:%s:%d %ld ", __FILE__, __func__, __LINE__, time(NULL));\
        fprintf(stderr, format, ## __VA_ARGS__);\
        fputc('\n', stderr);\
        }
    #define tracel {\
        fprintf(stderr,"%s:%s:%d %ld\n", __FILE__, __func__, __LINE__, time(NULL));\
        }
    #define trace tracel
#else
//#if !defined(tracel)
    #define tracef(...)
    #define tracel
    #define trace tracel
//#endif
#endif
