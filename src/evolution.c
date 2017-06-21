#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "xray_heating_lya_coupling/phys_const.h"
#include "xray_heating_lya_coupling/chem_const.h"
#include "xray_heating_lya_coupling/cosmology.h"
#include "xray_heating_lya_coupling/collisional_coupling.h"
#include "xray_heating_lya_coupling/xray_heating.h"
#include "xray_heating_lya_coupling/wouthuysen_effect.h"
#include "xray_heating_lya_coupling/21cm_grid.h"
#include "xray_heating_lya_coupling/spin_temperature.h"
#include "xray_heating_lya_coupling/3cm_grid.h"
#include "xray_heating_lya_coupling/spin_temperature_He.h"
#include "xray_heating_lya_coupling/compute_grid.h"

#include "ion_and_temp_evolution/cell.h"
#include "ion_and_temp_evolution/recomb_rates.h"
#include "ion_and_temp_evolution/evolution_loop.h"
#include "ion_and_temp_evolution/solve_ionfraction.h"
#include "ion_and_temp_evolution/solve_temperature.h"
#include "ion_and_temp_evolution/evolution_loop.h"
#include "ion_and_temp_evolution/solve_ion_temp.h"

#include "confObj.h"
#include "inputs.h"
#include "evolution.h"
#include "utils.h"
#include "settings.h"

void read_update_density_grid(int snap, int double_precision, grid_21cm_t *this21cmGrid, confObj_t simParam)
{
    char igm_density_file[256];
    char snap_string[8];
    
    for(int i=0; i<256; i++)
    {
        igm_density_file[i] = '\0';
    }
    
    sprintf(snap_string,"%03d",snap); 
    strcat(igm_density_file, simParam->igm_density_file);
    strcat(igm_density_file, "_");
    strcat(igm_density_file, snap_string);
    
    if(file_exist(igm_density_file) == 1)
    {
        read_density_21cmgrid(this21cmGrid, igm_density_file, double_precision);
    }else{
        fprintf(stderr, "No density file available, or names are incorrect!\n");
        exit(EXIT_FAILURE);
    }
}

void read_update_xray_grid(int snap, int double_precision, xray_grid_t *thisXray_grid, confObj_t simParam)
{
    char xray_file[256];
    char snap_string[8];
    
    for(int i=0; i<256; i++)
    {
        xray_file[i] = '\0';
    }
    
    sprintf(snap_string,"%03d",snap); 
    strcat(xray_file, simParam->xray_file);
    strcat(xray_file, "_");
    strcat(xray_file, snap_string);
    
    if(file_exist(xray_file) == 1)
    {
        read_lum_xraygrid(thisXray_grid, xray_file, double_precision);
    }else{
        fprintf(stderr, "No xray file available, or names are incorrect!\n");
        exit(EXIT_FAILURE);
    }
}

void update_xray_spectrum(xray_spectrum_t *thisXray_spectrum, double lumX, double alphaX, double nuX_min)
{
    thisXray_spectrum->lumX = lumX;
    thisXray_spectrum->alphaX = alphaX;
    thisXray_spectrum->nu_min = nuX_min;
}

void read_update_lya_grid(int snap, int double_precision, lya_grid_t *thisLya_grid, confObj_t simParam)
{
    char lya_file[256];
    char snap_string[8];
    
    for(int i=0; i<256; i++)
    {
        lya_file[i] = '\0';
    }
    
    sprintf(snap_string,"%03d",snap); 
    strcat(lya_file, simParam->lya_file);
    strcat(lya_file, "_");
    strcat(lya_file, snap_string);
    
    if(file_exist(lya_file) == 1)
    {
        read_lum_lyagrid(thisLya_grid, lya_file, double_precision);
    }else{
        fprintf(stderr, "No lya file available, or names are incorrect!\n");
        exit(EXIT_FAILURE);
    }
}

void update_lya_spectrum(lya_spectrum_t *thisLya_spectrum, double lumLya, double alphaLya, double nuLya_min)
{
    thisLya_spectrum->lum = lumLya;
    thisLya_spectrum->alpha = alphaLya;
    thisLya_spectrum->nu_min = nuLya_min;
}

void update_21cmgrid(grid_21cm_t *this21cmGrid, cell_t *thisCell, int x, int y, int z)
{
    int nbins = this21cmGrid->nbins;
    
    this21cmGrid->temp[x*nbins*nbins+y*nbins+z] = thisCell->temp;
    this21cmGrid->XHI[x*nbins*nbins+y*nbins+z] = thisCell->XHI;
    this21cmGrid->XHeI[x*nbins*nbins+y*nbins+z] = thisCell->XHeI;
    this21cmGrid->XHeII[x*nbins*nbins+y*nbins+z] = thisCell->XHeII;
    
    this21cmGrid->Xe[x*nbins*nbins+y*nbins+z] = 1.-thisCell->XHI;
}

void BubbleSort(double a[], int array_size)
{
    double temp;
    for(int i = 0; i < (array_size - 1); ++i)
    {
        for(int j = 0; j < array_size - 1 - i; ++j )
        {
            if(a[j] > a[j+1])
            {
                temp = a[j+1];
                a[j+1] = a[j];
                a[j] = temp;
            }
        }
    }
}   

double *create_redshift_table(inputlist_t *thisInputlist, double zstart, double zend, double dz)
{
    int num_tmp = thisInputlist->num + (zstart-zend)/dz + 1;
    int num;
    double *redshift;
    double *tmp;
    double z; 
    double tmp_value;
    int snap;
    
    tmp  = malloc(num_tmp*sizeof(double));
    
    snap = 0;
    for(int i=0; i<thisInputlist->num; i++)
    {
        tmp[i] = thisInputlist->redshift[snap];
        snap++;
    }
    z = zstart;
    for(int i=thisInputlist->num; i<num_tmp; i++)
    {
        tmp[i] = z;
        z = z-dz;
    }

    BubbleSort(tmp, num_tmp);
        
    num = num_tmp;
    for(int i=0; i<num_tmp-1; i++)
    {
        if(tmp[i] == tmp[i+1] || tmp[i] < zend)
        {
            num--;
        }
    }

    redshift = malloc(num*sizeof(double));
    
    int counter = 0;
    for(int i=0; i<num_tmp; i++)
    {
        tmp_value = tmp[num_tmp-i-1];
        redshift[counter] = tmp_value;

        if(tmp[num_tmp-i-1] == tmp[num_tmp-i-2])
        {
            i++;
        }
        
        if(redshift[counter] <= zend) break;
        counter++;
    }

    free(tmp);

    debug_printf(DEBUG_REDSHIFTLIST, "++++\nList of redshift steps (%d):\n", num);
    for(int i=0; i<num; i++)
    {
        debug_printf(DEBUG_REDSHIFTLIST, " redshift[%d] = %e\n", i, redshift[i]);
    }
    debug_printf(DEBUG_REDSHIFTLIST, "++++\n\n");
        
    thisInputlist->num_total = num;
    
    return redshift;
}

void evolve(confObj_t simParam, int myRank)
{
/*-------------------------------------------------------------------------------------*/
/* OUTPUT */
/*-------------------------------------------------------------------------------------*/ 
    char Tb_filename[128];
    char Ts_filename[128];
    char Xe_filename[128];
    char temp_filename[128];
    
    char redshift_string[6];

/*-------------------------------------------------------------------------------------*/
/* COSMOLOGY */
/*-------------------------------------------------------------------------------------*/
    double h = simParam->h;//0.7;
    double omega_m = simParam->omega_m;//0.27;
    double omega_b = simParam->omega_b;//0.045;
    double omega_l = simParam->omega_l;//0.73;
    
    double Y = simParam->Y;//0.24;
    
/*-------------------------------------------------------------------------------------*/
/* GRID */
/*-------------------------------------------------------------------------------------*/
    int nbins = simParam->grid_size;
    double box_size = simParam->box_size;
    
    inputlist_t *thisInputlist;
    thisInputlist = read_inputlist(simParam->input_list);
    
    int snap = 0;
    
/*-------------------------------------------------------------------------------------*/
/* XRAY Spectrum */
/*-------------------------------------------------------------------------------------*/
    double alphaX = thisInputlist->alphaX[snap]-1.;//-1.5;
    double nuX_min = thisInputlist->nuX_min[snap];//0.1e3*ev_to_erg/planck_cgs;
    double lumX = thisInputlist->lumX[snap];//1.29e22/(1.e3*ev_to_erg)*pow(1.e3*ev_to_erg/planck_cgs,-alphaX)*1.e-5;
    if(myRank == 0) printf("X-ray sources:\t alphaX = %e\t\t nuX_min = %e\t\tlumX = %e\n \t\t emissivity = lumX/(h nuX_min) * (nu / nuX_min)^alphaX\n", alphaX, nuX_min, lumX);
    
/*-------------------------------------------------------------------------------------*/
/* LYA Spectrum */
/*-------------------------------------------------------------------------------------*/
    double alphaLya = thisInputlist->alphaLya[snap]-1.;//-1.5;
    double nuLya_min = thisInputlist->nuLya_min[snap];//10.3*ev_to_erg/planck_cgs;
    double lumLya = thisInputlist->lumLya[snap];//3.4e40/(1.e3*ev_to_erg);
    if(myRank == 0) printf("Lya sources:\t alphaLya = %e\t nuLya_min = %e\t lumLya = %e\n \t\t emissivity = lumLya/(h nuLya_min) * (nu / nuLya_min)^alphaLya\n", alphaLya, nuLya_min, lumLya);

/*-------------------------------------------------------------------------------------*/
/* Initialize for ionization and temperature evolution */
/*-------------------------------------------------------------------------------------*/
    double zstart = simParam->zstart;
    double zend = simParam->zend;
    double dz = simParam->dz_max;
    double z = 0., z_new = 0.;
    
    double Xe = 0.;
    double temp = 0.;
  
    
    if(myRank == 0) printf("Initializing ... ");
    cell_t *cell;
    cell = initCell_temp(simParam->igm_temperature_start, 1.);
    
    recomb_t * recomb_rates;
    recomb_rates = calcRecRate();
        
/*-------------------------------------------------------------------------------------*/
/* Initialize for 21cm computation */
/*-------------------------------------------------------------------------------------*/
    
    cosmology_t *thisCosmology;
    
    xray_grid_t *thisXray_grid;
    xray_spectrum_t *thisXray_spectrum;
    
    lya_grid_t *thisLya_grid;
    lya_spectrum_t *thisLya_spectrum;
    
    grid_21cm_t *this21cmGrid;
    grid_3cm_t *this3cmGrid = NULL;
    
    k10_t *k10_table;
    
    int double_precision = simParam->input_doubleprecision;
    
    thisCosmology = allocate_cosmology(h, omega_m, omega_l, omega_b, zstart, Y);
    
    thisXray_grid = allocate_xray_grid(nbins, box_size);
    thisXray_spectrum = allocate_xray_spectrum(thisInputlist->lumX[snap], thisInputlist->alphaX[snap], thisInputlist->nuX_min[snap]);
    
    thisLya_grid = allocate_lya_grid(nbins, box_size);
    thisLya_spectrum = allocate_lya_spectrum(thisInputlist->lumLya[snap], thisInputlist->alphaLya[snap], thisInputlist->nuLya_min[snap]);
    
    this21cmGrid = allocate_21cmgrid(nbins, box_size);
    if(simParam->solve_3He == 1) this3cmGrid = allocate_3cmgrid(nbins, box_size);

    read_update_density_grid(thisInputlist->dens_snap[snap], double_precision, this21cmGrid, simParam);
    read_update_xray_grid(thisInputlist->xray_snap[snap], double_precision, thisXray_grid, simParam);
    read_update_lya_grid(thisInputlist->lya_snap[snap], double_precision, thisLya_grid, simParam);

    set_temperature_21cmgrid(this21cmGrid, simParam->igm_temperature_start);
    set_temperature_3cmgrid(this3cmGrid, 1.e4);

    k10_table = create_k10_table();
    
    if(myRank == 0) printf("done\n");
    
/*-------------------------------------------------------------------------------------*/
/* REDSHIFT EVOLUTION */
/*-------------------------------------------------------------------------------------*/
    
    this21cmGrid->z = zstart;
    this21cmGrid->z_old = zstart;
    
    if(snap<thisInputlist->num-1) snap++;
    
    double *redshift;
    
    redshift = create_redshift_table(thisInputlist, zstart, zend, dz);
    
    if(myRank == 0) printf("Number of steps: %d\n", thisInputlist->num_total);
    
    for(int i=0; i<thisInputlist->num_total-1; i++)
    {            
        z = redshift[i];
        if(myRank == 0) printf("++++\nStep %d at redshift z = %e\n", i, z);

        if(fabs(thisInputlist->redshift[snap]-z)<1.e-10 && snap<thisInputlist->num)
        {
            read_update_density_grid(thisInputlist->dens_snap[snap], double_precision, this21cmGrid, simParam);
            read_update_xray_grid(thisInputlist->xray_snap[snap], double_precision, thisXray_grid, simParam);
            update_xray_spectrum(thisXray_spectrum, thisInputlist->lumX[snap], thisInputlist->alphaX[snap], thisInputlist->nuX_min[snap]);
            read_update_lya_grid(thisInputlist->lya_snap[snap], double_precision, thisLya_grid, simParam);
            update_lya_spectrum(thisLya_spectrum, thisInputlist->lumLya[snap], thisInputlist->alphaLya[snap], thisInputlist->nuLya_min[snap]);
            
            snap++;
        }
        
        overwrite_old_values_21cmgrid(this21cmGrid);
        cosmology_update_z(thisCosmology, z);
        
        Xe = get_mean_Xe_21cmgrid(this21cmGrid);
        temp = get_mean_temp_21cmgrid(this21cmGrid);
        if(myRank == 0) printf(" <Xe> = %e\t <T> = %e\n", Xe, temp);
        
        /* Compute HI 21cm emission at redshift z */
        if(myRank == 0) printf(" computing HI 21cm emission\n");
        do_step_21cm_emission(thisCosmology, thisXray_grid, thisXray_spectrum, thisLya_grid, thisLya_spectrum, this21cmGrid, k10_table, Xe, myRank);
        
        /* Compute HeII 3cm emission at redshift z */
        if(simParam->solve_3He == 1)
        {
            if(myRank == 0) printf(" computing HeII 3cm emission\n");
            do_step_3cm_emission(thisCosmology, thisXray_spectrum, this21cmGrid, this3cmGrid, myRank);
        }
        
        /* Compute time evolution of Xe and temp from z=z_old to z=z_new */
        z_new = redshift[i+1];
        if(myRank == 0) printf(" updating ionization fraction and temperature from z = %e to %e\n", z, z_new);
        this21cmGrid->z = z_new;
        if(simParam->solve_3He == 1) this3cmGrid->z = z_new;
        compute_temp_ion_grid(this21cmGrid, recomb_rates, thisCosmology, thisXray_grid);
        
        
        /* OUTPUT */
        
        struct stat st;
        
        if (stat("output", &st) == -1) 
        {
            printf("creating directory\n");
            mkdir("output", 0755);
        }
        sprintf(redshift_string, "%3.1f", z);
        for(int j=0; j<128; j++)
        {
            Tb_filename[j] = '\0';
            Ts_filename[j] = '\0';
            Xe_filename[j] = '\0';
            temp_filename[j] = '\0';
        }
        
        if (stat("output/Tb", &st) == -1) 
        {
            printf("creating directory\n");
            mkdir("output/Tb", 0755);
        }
        strcat(Tb_filename, "output/Tb/");
        strcat(Tb_filename, "Tb_z");
        strcat(Tb_filename, redshift_string);
        strcat(Tb_filename, ".dat");
        write_Tb_field_file(this21cmGrid, Tb_filename);
        
        if (stat("output/Ts", &st) == -1) 
        {
            printf("creating directory\n");
            mkdir("output/Ts", 0755);
        }
        strcat(Ts_filename, "output/Ts/");
        strcat(Ts_filename, "Ts_z");
        strcat(Ts_filename, redshift_string);
        strcat(Ts_filename, ".dat");
        write_Ts_field_file(this21cmGrid, Ts_filename);
        
        if (stat("output/Xe", &st) == -1) 
        {
            printf("creating directory\n");
            mkdir("output/Xe", 0755);
        }
        strcat(Xe_filename, "output/Xe/");
        strcat(Xe_filename, "Xe_z");
        strcat(Xe_filename, redshift_string);
        strcat(Xe_filename, ".dat");
        write_Xe_field_file(this21cmGrid, Xe_filename);
        
        if (stat("output/temp", &st) == -1) 
        {
            printf("creating directory\n");
            mkdir("output/temp", 0755);
        }
        strcat(temp_filename, "output/temp/");
        strcat(temp_filename, "temp_z");
        strcat(temp_filename, redshift_string);
        strcat(temp_filename, ".dat");
        write_temp_field_file(this21cmGrid, temp_filename);
        
        if(simParam->solve_3He == 1) 
        {
            for(int j=0; j<128; j++)
            {
                Tb_filename[j] = '\0';
                Ts_filename[j] = '\0';
            }
            
            if (stat("output/Tb_He", &st) == -1) 
            {
                printf("creating directory\n");
                mkdir("output/Tb_He", 0755);
            }
            strcat(Tb_filename, "output/Tb_He/");
            strcat(Tb_filename, "Tb_He_z");
            strcat(Tb_filename, redshift_string);
            strcat(Tb_filename, ".dat");
            write_Tb_3cm_field_file(this3cmGrid, Tb_filename);
            
            if (stat("output/Ts_He", &st) == -1) 
            {
                printf("creating directory\n");
                mkdir("output/Ts_He", 0755);
            }
            strcat(Ts_filename, "output/Ts_He/");
            strcat(Ts_filename, "Ts_He_z");
            strcat(Ts_filename, redshift_string);
            strcat(Ts_filename, ".dat");
            write_Ts_3cm_field_file(this3cmGrid, Ts_filename);
        }
        if(myRank == 0) printf("++++\n\n");
    }
    
    if(myRank == 0) printf("deallocating ... ");
    free(redshift);
      
/*-------------------------------------------------------------------------------------*/
/* Deallocate ionization and temperature evolution structs */
/*-------------------------------------------------------------------------------------*/
    deallocate_recomb(recomb_rates);
    deallocate_cell(cell);
    
/*-------------------------------------------------------------------------------------*/
/* Deallocate 21cm computation structs */
/*-------------------------------------------------------------------------------------*/
    deallocate_cosmology(thisCosmology);
    deallocate_xray_grid(thisXray_grid);
    deallocate_xray_spectrum(thisXray_spectrum);
    deallocate_lya_grid(thisLya_grid);
    deallocate_lya_spectrum(thisLya_spectrum);
    deallocate_21cmgrid(this21cmGrid);
    if(simParam->solve_3He == 1)
    {
        deallocate_3cmgrid(this3cmGrid);
    }
    deallocate_k10_table(k10_table);
    deallocateInputlist(thisInputlist);
    
    if(myRank == 0) printf("done\n");
}


void evolve_old(int myRank)
{
/*-------------------------------------------------------------------------------------*/
/* COSMOLOGY */
/*-------------------------------------------------------------------------------------*/
    double h = 0.7;
    double omega_m = 0.27;
    double omega_b = 0.045;
    double omega_l = 0.73;
    
    double Y = 0.24;
    
/*-------------------------------------------------------------------------------------*/
/* XRAY Spectrum */
/*-------------------------------------------------------------------------------------*/
    double alphaX = -1.5;
    double nuX_min = 0.1e3*ev_to_erg/planck_cgs;
    double lumX = 3.4e40/(1.e3*ev_to_erg);
    printf("alphaX = %e\t nuX_min = %e\t lumX = %e\n", alphaX, nuX_min, lumX);
/*-------------------------------------------------------------------------------------*/
/* LYA Spectrum */
/*-------------------------------------------------------------------------------------*/
    double alphaLya = -1.5;
    double nuLya_min = 10.3*ev_to_erg/planck_cgs;
    double lumLya = 3.4e40/(1.e3*ev_to_erg);

/*-------------------------------------------------------------------------------------*/
/* GRID */
/*-------------------------------------------------------------------------------------*/
    int nbins = 128;
    int local_n0 = nbins;
    double box_size = 80.;
    
    char densfile[128];
    char xray_sourcefile[128];
    char lya_sourcefile[128];
    
    for (int i=0; i<128; i++)
    {
        densfile[i] = '\0';
        xray_sourcefile[i] = '\0';
        lya_sourcefile[i] = '\0';
    }
    strcat(densfile, "/home/anne/PostdocSwinburne/grid_model_mpi/inputfiles/grid128/density11_ic.in");
    strcat(xray_sourcefile, "/home/anne/PostdocSwinburne/grid_model_mpi/inputfiles/grid128/density11_ic.in");
    strcat(lya_sourcefile, "/home/anne/PostdocSwinburne/grid_model_mpi/inputfiles/grid128/density11_ic.in");
    
/*-------------------------------------------------------------------------------------*/
/* Initialize for ionization and temperature evolution */
/*-------------------------------------------------------------------------------------*/
    double zstart = 25.;
    double zend = 10.;
    double dz = 0.5;
  
    cell_t *cell;
    cell = initCell_temp(14., 1.);
    
    recomb_t * recomb_rates;
    recomb_rates = calcRecRate();
    
    int dim_matrix = 2;
    
/*-------------------------------------------------------------------------------------*/
/* Initialize for 21cm computation */
/*-------------------------------------------------------------------------------------*/

    cosmology_t *thisCosmology;
    
    xray_grid_t *thisXray_grid;
    xray_spectrum_t *thisXray_spectrum;
    
    lya_grid_t *thisLya_grid;
    lya_spectrum_t *thisLya_spectrum;
    
    grid_21cm_t *this21cmGrid;
    
    k10_t *k10_table;
    
    int double_precision = 0;
    
    thisCosmology = allocate_cosmology(h, omega_m, omega_l, omega_b, zstart, Y);
    
    thisXray_grid = allocate_xray_grid(nbins, box_size);
    thisXray_spectrum = allocate_xray_spectrum(lumX, alphaX, nuX_min);
    
    thisLya_grid = allocate_lya_grid(nbins, box_size);
    thisLya_spectrum = allocate_lya_spectrum(lumLya, alphaLya, nuLya_min);
    
    this21cmGrid = allocate_21cmgrid(nbins, box_size);
    
    read_density_21cmgrid(this21cmGrid, densfile, double_precision);
    read_lum_xraygrid(thisXray_grid, xray_sourcefile, double_precision);
    read_lum_lyagrid(thisLya_grid, lya_sourcefile, double_precision);
    
    set_temperature_21cmgrid(this21cmGrid, 14.);
    
    k10_table = create_k10_table();
    
/*-------------------------------------------------------------------------------------*/
/* REDSHIFT EVOLUTION */
/*-------------------------------------------------------------------------------------*/

    const double tmp_HI = 1./(boltzman_cgs*nu_HI);
    const double tmp_HeI = 1./(boltzman_cgs*nu_HeI);
    const double tmp_HeII = 1./(boltzman_cgs*nu_HeII);
    
    double fion_HI = 1.;
    double fion_HeI = 1.;
    double fion_HeII = 1.;
    
    double fheat = 1.;
    
    double Xe = 0.;
    
    double dz_tmp = dz;
    
    for(double z = zstart; z>zend; z = z-dz_tmp)
    {
        printf("z = %e\t %e\n", z, fmod(z,dz));
        cosmology_update_z(thisCosmology, z);
        if(fmod(z,dz) == 0.)
        {
            Xe = get_mean_Xe_21cmgrid(this21cmGrid);
            printf("Xe = %e\n", Xe);
            do_step_21cm_emission(thisCosmology, thisXray_grid, thisXray_spectrum, thisLya_grid, thisLya_spectrum, this21cmGrid, k10_table, Xe, myRank);
        }
        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    if(fmod(z,dz) == 0.)
                    {
                        double XHII = 1.-creal(this21cmGrid->XHI[i*nbins*nbins+j*nbins+k]);
                        double XHeII = creal(this21cmGrid->XHeII[i*nbins*nbins+j*nbins+k]);
                        double XHeIII = 1.-XHeII-creal(this21cmGrid->XHeI[i*nbins*nbins+j*nbins+k]);
                        update_X_cell(cell, XHII, XHeII, XHeIII);
                        
                        double dens = creal(this21cmGrid->dens[i*nbins*nbins+j*nbins+k]);
                        update_dens_cell(cell, dens);
                        
                        double temp = creal(this21cmGrid->temp[i*nbins*nbins+j*nbins+k]);
                        update_temp_cell(cell, temp);
                        
                        double heatHI = fheat*creal(thisXray_grid->xray_heating_HI[i*nbins*nbins+j*nbins+k]);
                        double heatHeI = fheat*creal(thisXray_grid->xray_heating_HeI[i*nbins*nbins+j*nbins+k]);
                        double heatHeII = fheat*creal(thisXray_grid->xray_heating_HeII[i*nbins*nbins+j*nbins+k]);
                        update_heat_cell(cell, heatHI, heatHeI, heatHeII);
                        
                        double photIonHI = creal(thisXray_grid->xray_ionization_HI[i*nbins*nbins+j*nbins+k]) + heatHI*fion_HI*tmp_HI;
                        double photIonHeI = creal(thisXray_grid->xray_ionization_HeI[i*nbins*nbins+j*nbins+k]) + heatHeI*fion_HeI*tmp_HeI;
                        double photIonHeII = creal(thisXray_grid->xray_ionization_HeII[i*nbins*nbins+j*nbins+k]) + heatHeII*fion_HeII*tmp_HeII;
                        update_photIon_cell(cell, photIonHI, photIonHeI, photIonHeII);
                        
//                         printf("photHI = %e\n", cell->photIonHI);
                    }

                    calc_step(recomb_rates, cell, dim_matrix, z, dz_tmp);
                    printf("%d %d %d: z = %e:\tT = %e\t T = T0/(1+z)^2 = %e\t XHII = %e\t XHeII = %e\t XHeIII = %e\n", i,j,k,z, cell->temp, T_CMB(zstart)/((1.+zstart)*(1.+zstart))*((1.+z-dz_tmp)*(1.+z-dz_tmp)), cell->XHII, cell->XHeII, cell->XHeIII);
                    
                    if(fmod(z,dz) == 0.)
                    {
                        update_21cmgrid(this21cmGrid, cell, i, j, k);
                    }
                }
            }
        }
    }
      
      
/*-------------------------------------------------------------------------------------*/
/* Deallocate ionization and temperature evolution structs */
/*-------------------------------------------------------------------------------------*/
    deallocate_recomb(recomb_rates);
    deallocate_cell(cell);
    
/*-------------------------------------------------------------------------------------*/
/* Deallocate 21cm computation structs */
/*-------------------------------------------------------------------------------------*/
    deallocate_cosmology(thisCosmology);
    
    deallocate_xray_grid(thisXray_grid);
    deallocate_xray_spectrum(thisXray_spectrum);
    
    deallocate_lya_grid(thisLya_grid);
    deallocate_lya_spectrum(thisLya_spectrum);
    
    deallocate_21cmgrid(this21cmGrid);

    deallocate_k10_table(k10_table);
}

