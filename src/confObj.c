/*
 *  confObj.c
 *  uvff
 *
 *  Created by Adrian Partl on 4/15/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

/*--- Includes ----------------------------------------------------------*/
#include "confObj.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include "xmem.h"

/*--- Defines for the Ini structure -------------------------------------*/


/*--- Prototypes of local functions -------------------------------------*/


/*--- Implementations of exported functios ------------------------------*/
extern confObj_t
confObj_new(parse_ini_t ini)
{
	confObj_t config;
	assert(ini != NULL);
	
	config = xmalloc(sizeof(struct confObj_struct));
	
	//reading mandatory stuff
	getFromIni(&(config->input_doubleprecision), parse_ini_get_int32,
	           ini, "inputFilesAreInDoublePrecision", "General");
	getFromIni(&(config->input_list), parse_ini_get_string,
	           ini, "inputList", "General");
	getFromIni(&(config->igm_density_file), parse_ini_get_string,
	           ini, "inputIgmDensityFile", "General");
	getFromIni(&(config->xray_file), parse_ini_get_string,
	           ini, "inputXraySourcesFile", "General");
	getFromIni(&(config->lya_file), parse_ini_get_string,
	           ini, "inputLyaSourcesFile", "General");
    
	getFromIni(&(config->grid_size), parse_ini_get_int32,
            ini, "gridsize", "General");
	getFromIni(&(config->box_size), parse_ini_get_double,
            ini, "boxsize", "General");
    
	getFromIni(&(config->zstart), parse_ini_get_double,
		   ini, "initialRedshift", "General");
	getFromIni(&(config->zend), parse_ini_get_double,
		   ini, "finalRedshift", "General");
    getFromIni(&(config->dz_max), parse_ini_get_double,
		   ini, "dz_max", "General");
	
	getFromIni(&(config->h), parse_ini_get_double,
	           ini, "h", "General");
	getFromIni(&(config->omega_b), parse_ini_get_double,
	           ini, "omega_b", "General");
	getFromIni(&(config->omega_m), parse_ini_get_double,
	           ini, "omega_m", "General");
	getFromIni(&(config->omega_l), parse_ini_get_double,
	           ini, "omega_l", "General");
	getFromIni(&(config->sigma8), parse_ini_get_double,
	           ini, "sigma8", "General");
	getFromIni(&(config->Y), parse_ini_get_double,
		   ini, "Y", "General");

	return config;
}

extern void
confObj_del(confObj_t *config)
{
	assert(config != NULL);
	assert(*config != NULL);
	
    xfree((*config)->input_list);
	xfree((*config)->igm_density_file);
    xfree((*config)->xray_file);
	xfree((*config)->lya_file);
	xfree(*config);
	*config = NULL;
}

extern confObj_t
readConfObj(char *fileName) {
	confObj_t newConfObj;
	parse_ini_t ini;
	
	ini = parse_ini_open(fileName);
	if (ini == NULL) {
		fprintf(stderr, "FATAL:  Could not open %s for reading.\n",
				fileName);
		exit(EXIT_FAILURE);
	}
	
	newConfObj = confObj_new(ini);
	
	parse_ini_close(&ini);
	
	return newConfObj;
}

/*--- Implementations of local functions --------------------------------*/
