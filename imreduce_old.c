/*
  This file was downloaded from the CFITSIO utilities web page:
    http://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html

  That page contains this text:
    You may freely modify, reuse, and redistribute these programs as you wish.

  We assume it was originally written by the CFITSIO authors (primarily William
  D. Pence).

  We (the Astrometry.net team) have modified it slightly.
  # Licensed under a 3-clause BSD style license - see LICENSE


*/

#include <string.h>
#include <stdio.h>
#include "fitsio.h"

int main(int argc, char *argv[])
{
    fitsfile *afptr, *bfptr, *cfptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, bnaxis, check = 1, ii, op, bitpix=-32;
    long npixels = 1, firstpix[2] = {1,1};
    long anaxes[2] = {1,1}, bnaxes[2]={1,1}, cnaxes[2]={1,1};
    double *apix, *bpix, *cpix;
    float flat_min=0.;

    if (argc != 5) { 
      printf("Usage: imreduce datafile biasfile flatfile\n");
      printf("\n");
      printf("Just reduce an image\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */
    fits_open_file(&bfptr, argv[2], READONLY, &status);
    fits_open_file(&cfptr, argv[3], READONLY, &status);

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_dim(bfptr, &bnaxis, &status);
    fits_get_img_size(afptr, 2, anaxes, &status);
    fits_get_img_size(bfptr, 2, bnaxes, &status);
    fits_get_img_size(cfptr, 2, cnaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }
         /* check that the input 2 images have the same size */
    else if ( anaxes[0] != bnaxes[0] || 
              anaxes[0] != cnaxes[0] || 
              anaxes[1] != bnaxes[1] || 
              anaxes[1] != cnaxes[1] ) {
       printf("Error: input images don't have same size\n");
       check = 0;
    }

    /* create the new empty output file if the above checks are OK */
    if (check && !fits_create_file(&outfptr, argv[4], &status) )
    {
      /* copy all the header keywords from first image to new output file */
      fits_copy_header(afptr, outfptr, &status);
      fits_update_key(outfptr, TINT, "BITPIX", &bitpix, NULL, &status);
      fits_delete_key(outfptr, "BZERO", &status);
      fits_delete_key(outfptr, "BSCALE", &status);

      fits_read_key(cfptr, TFLOAT, "FLATMIN", &flat_min, NULL, &status);
      if (status!=0) {
          status=0;
          flat_min=0.;
      }

      npixels = anaxes[0];  /* no. of pixels to read in each row */

      apix = (double *) malloc(npixels * sizeof(double)); /* mem for 1 row */
      bpix = (double *) malloc(npixels * sizeof(double));
      cpix = (double *) malloc(npixels * sizeof(double));

      if (apix == NULL || bpix == NULL || cpix == NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      /* loop over all rows of the plane */
      for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++)
      {
        /* Read both images as doubles, regardless of actual datatype.  */
        /* Give starting pixel coordinate and no. of pixels to read.    */
        /* This version does not support undefined pixels in the image. */

        if (fits_read_pix(afptr, TDOUBLE, firstpix, npixels, NULL, apix,
                          NULL, &status)  ||         
            fits_read_pix(bfptr, TDOUBLE, firstpix, npixels,  NULL, bpix,
                          NULL, &status)  ||        
            fits_read_pix(cfptr, TDOUBLE, firstpix, npixels,  NULL, cpix,
                          NULL, &status)  )        
            break;   /* jump out of loop on error */

        for(ii=0; ii< npixels; ii++) {
          if (cpix[ii]>flat_min) cpix[ii] = (apix[ii]-bpix[ii])/cpix[ii];
          else cpix[ii]=0;
        }

        fits_write_pix(outfptr, TDOUBLE, firstpix, npixels,
                     cpix, &status); /* write new values to output image */
      }

      //fits_update_key(outfptr, TINT, "BITPIX", &bitpix, NULL, &status);
      fits_close_file(outfptr, &status);
      free(apix);
      free(bpix);
      free(cpix);
    }

    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);
    fits_close_file(cfptr, &status);
 
    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}

