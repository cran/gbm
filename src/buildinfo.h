//  GBM by Greg Ridgeway  Copyright (C) 2003
//  License:    GNU GPL (version 2 or later)

#ifndef BUILDINFO_H
#define BUILDINFO_H

    #include <R.h>
    #define RSWERRORPRINTF Rprintf(
    #define ErrorTrace(p)
    #define delete_item(p)
    #define new_item(p)
    #define mem_shut_down()
    #define FAILED(hr) ((unsigned long)hr != 0)

    typedef unsigned long HRESULT;
    #define S_OK 0
    #define E_FAIL 1
    #define E_INVALIDARG 2
    #define E_OUTOFMEMORY 3
    #define E_INVALID_DATA 4
    #define E_NOTIMPL 5

    #define LEVELS_PER_CHUNK ((unsigned long) 1)

    typedef unsigned long ULONG;
    typedef char *PCHAR;

    // #define NOISY_DEBUG

#endif // BUILDINFO_H
