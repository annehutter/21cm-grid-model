/*
 *  confObj.h
 *  uvff
 *
 *  Created by Adrian Partl on 4/15/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONFOBJ_H
#define CONFOBJ_H

/*--- Includes ----------------------------------------------------------*/
#include "parse_ini.h"
#include <stdint.h>
#include <stdbool.h>


/*--- ADT handle --------------------------------------------------------*/
typedef struct confObj_struct *confObj_t;


/*--- Implemention of main structure ------------------------------------*/
struct confObj_struct {
	//General
	int             input_doubleprecision;
    char            *input_list;
	char 			*igm_density_file;
    char            *xray_file;
    char            *lya_file;
	
	int             grid_size;
	double			box_size;
	
    double          zstart;
    double          zend;
    double          dz_max;
	
	double			h;
	double			omega_b;
	double 			omega_m;
	double 			omega_l;
	double			sigma8;
	double			Y;
    
    int             solve_3He;
    
    double          igm_temperature_start;
};


/*--- Prototypes of exported functions ----------------------------------*/
extern confObj_t
readConfObj(char *fileName);

extern confObj_t
confObj_new(parse_ini_t ini);

extern void
confObj_del(confObj_t *config);


#endif
