#ifndef INPUTS_H
#define INPUTS_H
#endif

typedef struct
{
    int num;
    
    double *redshift;
    
    int *dens_snap;
    int *xray_snap;
    int *lya_snap;
    
    double *lumX;
    double *alphaX;
    double *nuX_min;
    
    double *lumLya;
    double *alphaLya;
    double *nuLya_min;
} inputlist_t;

inputlist_t *read_inputlist(char *filename);
inputlist_t *initInputlist(int num);
void deallocateInputlist(inputlist_t *thisInputlist);
