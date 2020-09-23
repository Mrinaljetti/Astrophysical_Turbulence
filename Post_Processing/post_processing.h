//
//  Post_processing.h
//  Astrophysical_Turbulence
//
//  Created by Mrinal Jetti on 09/09/20.
//  Copyright © 2020 Mrinal Jetti. All rights reserved.
//

#ifndef post_processing_h
#define post_processing_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parameters.h"

void read_dbl(int file_num, double *data_array)
{
    FILE *fp;
    int i;
    char filenumb[5];
    char filename[100];
    sprintf(filenumb,"%04d",file_num);
    strcpy(filename,datdir);
    strcat(filename,prefix);
    strcat(filename,filenumb);
    strcat(filename,suffix);
    
    printf("%s\n",filename);
    fp = fopen(filename,"rb");
    
    if(fp == NULL)
    {
      printf("fopen error (1)");
      exit(0);
    }
    
    for(i=0;i<(nv)*nx*ny*nz;i++){
      fread(&data_array[i],sizeof(double),1,fp);
    }
    fclose(fp);
    return;
}




#endif /* Post_processing_h */
