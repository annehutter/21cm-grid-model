#ifndef SETTINGS_H
#define SETTINGS_H
#endif 

#define debug_printf(DEBUG, fmt, ...) \
            do { if (DEBUG) fprintf(stderr, fmt, ##__VA_ARGS__); } while (0)
                
#ifdef DEBUGFLAG_REDSHIFTLIST
#define DEBUG_REDSHIFTLIST 1
#else
#define DEBUG_REDSHIFTLIST 0
#endif
                
#ifdef DEBUGFLAG_21CM_TB_CALCULATION
#define DEBUG_21CM_TB_CALCULATION 1
#else
#define DEBUG_21CM_TB_CALCULATION 0
#endif
                
#ifdef DEBUGFLAG_3CM_TB_CALCULATION
#define DEBUG_3CM_TB_CALCULATION 1
#else
#define DEBUG_3CM_TB_CALCULATION 0
#endif

#ifdef DEBUGFLAG_TS_CALCULATION
#define DEBUG_TS_CALCULATION 1
#else
#define DEBUG_TS_CALCULATION 0
#endif

#ifdef DEBUGFLAG_TS_HE_CALCULATION
#define DEBUG_TS_HE_CALCULATION 1
#else
#define DEBUG_TS_HE_CALCULATION 0
#endif
                
#ifdef DEBUGFLAG_XRAY_FILTER
#define DEBUG_XRAY_FILTER 1
#else
#define DEBUG_XRAY_FILTER 0
#endif

#ifdef DEBUGFLAG_LYA_FILTER
#define DEBUG_LYA_FILTER 1
#else
#define DEBUG_LYA_FILTER 0
#endif
                
#ifdef DEBUGFLAG_LYA_SPECTRUM
#define DEBUG_LYA_SPECTRUM 1
#else
#define DEBUG_LYA_SPECTRUM 0
#endif
