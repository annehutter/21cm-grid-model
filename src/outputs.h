#ifndef OUTPUTS_H
#define OUTPUTS_H
#endif

typedef struct
{
    int num;
    
    double *redshift;
}outputlist_t;

outputlist_t *read_outputlist(char *filename);
outputlist_t *initOutputlist(int num);
void deallocateOutputlist(outputlist_t *thisOutputlist);
