//
//  Post_processing.c
//  Astrophysical_Turbulence
//
//  Created by Mrinal Jetti on 09/09/20.
//  Copyright © 2020 Mrinal Jetti. All rights reserved.
//

#include "post_processing.h"
#include "parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
    double *all_data_array = (double *)malloc(sizeof(double)*((nv)*(nx*ny*nz)));
    double *rho = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *vx1 = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *vx2 = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *vx3 = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *prs = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *temp = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *mach = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *sb = (double *)malloc(sizeof(double)*(nx*ny));
    double *vlos = (double *)malloc(sizeof(double)*(nx*ny));
    double *sigvlos = (double *)malloc(sizeof(double)*(nx*ny));
    
    double cs;
    double num_density;
    double radiation_rate;
    double sb_rms;
    double sb_avg;
    double temp_avg;
    double mach_rms;

    int i,j,k,i1,j1,k1;

    FILE *fp;
    FILE *fpsb;
    FILE *ftemp;
    FILE *fmrms;
    char filenum[5];
    char filename[100];
    
    fpsb = fopen(sb_data_dir,"w");
    fprintf(fpsb,"del_sb_rms/sb_avg\n");
    ftemp = fopen(temp_data_dir,"w");
    fprintf(ftemp,"temp\n");
    fmrms = fopen(mrms_data_dir,"w");
    fprintf(fmrms,"mrms\n");

    for(i=f1;i<f2;i++)
    {
        read_dbl(i, all_data_array);

        sprintf(filenum,"%04d",i);
        strcpy(filename, datdir );
        strcat(filename,data_prefix);
        strcat(filename,filenum);
        strcat(filename,data_suffix);

        fp = fopen( filename ,"w");

        fprintf(fp,"rho\tvx1\tvx2\tvx3\tprs\t\n");

        if(fp==NULL)
        {
            printf("Error!!");
            exit(1);
        }
        
        temp_avg = 0.0;
        

        for(j=0;j<(nx*ny*nz);j++)
        {
            rho[j] = all_data_array[j];
            vx1[j] = all_data_array[j+nx*ny*nz];
            vx2[j] = all_data_array[j+2*nx*ny*nz];
            vx3[j] = all_data_array[j+3*nx*ny*nz];
            prs[j] = all_data_array[j+4*nx*ny*nz];
            
            temp[j] = (prs[j]/rho[j])*((CONST_mu*UNIT_VELOCITY*UNIT_VELOCITY*CONST_mp)/CONST_kB);
            temp_avg += temp[j];
            cs = sqrt(gamma*(prs[j])/rho[j]);
            mach[j] = sqrt(vx1[j]*vx1[j]+vx2[j]*vx2[j]+vx3[j]*vx3[j])/cs;
            mach_rms += mach[j]*mach[j];
            
            if(temp[j]<0.)
            {
                printf("error here %lf %lf %d\n", rho[j],prs[j],j);
                exit(0);
            }
            
            fprintf(fp,"%f\t%f\t%f\t%f\t%f\t",rho[j],vx1[j],vx2[j],vx3[j],prs[j]);
            fprintf(fp,"\n");
        }
        mach_rms = mach_rms/(nx*ny*nz);
        temp_avg = temp_avg/(nx*ny*nz);
        fprintf(ftemp,"%f\n",temp_avg);
        fprintf(fmrms,"%f\n",mach_rms);
        
        fclose(fp);
        
        sb_rms = 0.0;
        sb_avg = 0.0;
        for(i1=0;i1<nx;i1++)
        {
            for(j1=0;j1<ny;j1++)
            {
                sb[i1*ny+j1] = 0.0;
                vlos[i1*ny+j1] = 0.0;
                for(k1=0;k1<nz;k1++)
                {
                    num_density = rho[i1*(ny*nz)+(j1*nz+k1)]*UNIT_DENSITY/CONST_mu/CONST_mp;
                    radiation_rate = num_density*num_density*lambda(temp[i1*(ny*nz)+(j1*nz+k1)]);
                    sb[i1*ny+j1] += radiation_rate*(1.0/(double)nz)*UNIT_LENGTH;
                    vlos[i1*ny+j1] += vx3[i1*(ny*nz)+(j1*nz+k1)]*UNIT_VELOCITY*radiation_rate*(1.0/(double)nz)*UNIT_LENGTH;
                }
                vlos[i1*ny+j1] = vlos[i1*ny+j1]/sb[i1*ny+j1];
                sb_rms += sb[i1*ny+j1]*sb[i1*ny+j1];
                sb_avg += sb[i1*ny+j1];
                
                for(k1=0;k1<nz;k1++)
                {
                    num_density = rho[i1*(ny*nz)+(j1*nz+k1)]*UNIT_DENSITY/CONST_mu/CONST_mp;
                    radiation_rate = num_density*num_density*lambda(temp[i1*(ny*nz)+(j1*nz+k1)]);
                    sigvlos[i1*ny+j1] += (vx3[i1*(ny*nz)+(j1*nz+k1)]-vlos[i1*ny+j1])*(vx3[i1*(ny*nz)+(j1*nz+k1)]-vlos[i1*ny+j1])*UNIT_VELOCITY*radiation_rate*(1.0/(double)nz)*UNIT_LENGTH;
                    
                }
                sigvlos[i1*ny+j1] = sigvlos[i1*ny+j1]/sb[i1*ny+j1];
            }
        }
        sb_avg = sb_avg/(nx*ny);
        sb_rms = (sb_rms/(nx*ny)) - pow(sb_avg,2);
        
        fprintf(fpsb,"%lf\n",sb_rms/sb_avg);
        
    }
    
    fclose(fpsb);
    fclose(fmrms);
    fclose(ftemp);
    return 1;
    

}

