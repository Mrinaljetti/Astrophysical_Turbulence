//
//  Post_processing.c
//  Astrophysical_Turbulence
//
//  Created by Mrinal Jetti on 09/09/20.
//  Copyright © 2020 Mrinal Jetti. All rights reserved.
//
#include "parameters.h"
#include "post_processing.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{

    //Variables, data stored at each point of the box
    double *all_data_array = (double *)malloc(sizeof(double)*((nv)*(nx*ny*nz)));
    double *rho = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *vx1 = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *vx2 = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *vx3 = (double *)malloc(sizeof(double)*(nx*ny*nz));
    double *prs = (double *)malloc(sizeof(double)*(nx*ny*nz));
    
    //Derived variables
    double *temp = (double *)malloc(sizeof(double)*(nx*ny*nz));//Temperature=(P/rho)*(mu*amu/Kb) multiply by UNIT_LENGTH^2 in units
    double *mach = (double *)malloc(sizeof(double)*(nx*ny*nz));//Mach = |v|/cs
    double *sb = (double *)malloc(sizeof(double)*(nx*ny));     //Surface brigthness = integral(n^2*lamda(T)*dz)
    double *vlos = (double *)malloc(sizeof(double)*(nx*ny));   //vlos = integral(vz*n^2*lamda(T)*dz)/sb
    double *sigvlos = (double *)malloc(sizeof(double)*(nx*ny));//sigvlos = integral((vz-vlos)^2*n^2*lamda(T)*dz)/sb
    
    //Derived variables
    double cs;
    double cs_avg;
    double num_density;
    double radiation_rate;  //(rho)/(mu*amu)
    double sb_rms;          //sb_rms = sqrt(<sb^2>-<sb>^2) where <> denotes average
    double sb_avg;          //sb average
    double temp_avg;        //temp average
    double v_rms;           //sqrt(<|v|^2>)
    double mach_rms;        //v_rms/cs_avg
    
    //NOTE sb_rms/sb_avg is stored in the file fpsb
    
    //Looping variables
    int i,j,i1,j1,k1;

    FILE *fp;               // This file stores rho,vx,vy,vz,prs
    FILE *fpsb;             // This file stores sb_rms/sb_avg,temp_avg, mach_rms : file name sb_data.txt
    char filenum[5];
    char filename[100];
    
    fpsb = fopen(sb_data_dir,"w");
    fprintf(fpsb,"del_sb_rms/sb_avg\ttemp_avg\tmrms\n");

    
    //This loop is over the .dbl files(snapshots of our box)
    for(i=f1;i<=f2;i++)
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
        
        
        //This loop is over all the points of our simulation box
        //rho,vx,vy,vz,prs variables are written to fp
        //temp_avg,mach_rms,v_rms are evaluated
        temp_avg = 0.0;
        mach_rms = 0.0;
        v_rms = 0.0;
        cs_avg = 0.0;
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
            cs_avg += cs;
            mach[j] = sqrt(vx1[j]*vx1[j]+vx2[j]*vx2[j]+vx3[j]*vx3[j])/cs;
            v_rms += vx1[j]*vx1[j]+vx2[j]*vx2[j]+vx3[j]*vx3[j];
            if(temp[j]<0.)
            {
                printf("error here %lf %lf %d\n", rho[j],prs[j],j);
                exit(0);
            }
            
            fprintf(fp,"%f\t%f\t%f\t%f\t%f\t",rho[j],vx1[j],vx2[j],vx3[j],prs[j]);
            fprintf(fp,"\n");

        }
        fclose(fp);
        cs_avg = cs_avg/(nx*ny*nz);
        v_rms = sqrt(v_rms/(nx*ny*nz));
        temp_avg = temp_avg/(nx*ny*nz);
        mach_rms = v_rms/cs_avg;
        
        //This loop is to calcualte sb,vlos by keeping in focus the integral to be evaluted along the z-direction
        sb_rms = 0.0;
        sb_avg = 0.0;
        for(i1=0;i1<nx;i1++)
        {
            for(j1=0;j1<ny;j1++)
            {
                sb[i1*ny+j1] = 0.0;
                vlos[i1*ny+j1] = 0.0;
                
                //integral along z-direction
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
                
                sigvlos[i1*ny+j1] = 0.0;
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
        sb_rms = (sb_rms/(nx*ny));    // - pow(sb_avg,2);
        sb_rms = sqrt(sb_rms - pow(sb_avg,2));
        
        //Write the variables to sb_data.txt
        fprintf(fpsb,"%lf\t%lf\t%lf\n",sb_rms/sb_avg,temp_avg,mach_rms);
        
    }
    fclose(fpsb);
    return 0;

}


