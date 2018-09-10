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
    fitsfile *afptr, *bfptr, *cfptr, *dfptr, *efptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, bnaxis, cnaxis, dnaxis, enaxis, check = 1, ii, op;
    long npixels = 1, firstpix[2] = {1,1}, lastpix[2] = {1,1};
    long anaxes[2] = {1,1}, bnaxes[2]={1,1}, cnaxes[2]={1,1}, dnaxes[2]={1,1}, enaxes[2]={1,1};
    double *apix, *bpix, *cpix, *dpix, *epix;
    float bscale,ibscale,srcfac=10.;
    double nsig=3.,tmp,s0=50.,s1, s12,nsig2,v00,v11,diff,vdiff,sky,sky0,sat,sat0;

    if (argc < 6) {
      printf("Usage: weight_clip datafile weightfile maskfile stack_file stack_weight [nsig] [srcfac] [s0]\n");
      printf("\n");
      printf("clip out bogux pixels, updates datafile and weightfile\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READWRITE, &status); /* open input images */
    fits_open_file(&bfptr, argv[2], READWRITE, &status);
    fits_open_file(&cfptr, argv[3], READONLY, &status);
    fits_open_file(&dfptr, argv[4], READONLY, &status);
    fits_open_file(&efptr, argv[5], READONLY, &status);

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_dim(bfptr, &bnaxis, &status);
    fits_get_img_dim(cfptr, &cnaxis, &status);
    fits_get_img_dim(dfptr, &dnaxis, &status);
    fits_get_img_dim(efptr, &enaxis, &status);
    fits_get_img_size(afptr, 2, anaxes, &status);
    fits_get_img_size(bfptr, 2, bnaxes, &status);
    fits_get_img_size(cfptr, 2, cnaxes, &status);
    fits_get_img_size(dfptr, 2, dnaxes, &status);
    fits_get_img_size(efptr, 2, enaxes, &status);

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
              anaxes[0] != dnaxes[0] ||
              anaxes[0] != enaxes[0] ||
              anaxes[1] != bnaxes[1] ||
              anaxes[1] != cnaxes[1] ||
              anaxes[1] != dnaxes[1] ||
              anaxes[1] != enaxes[1] ) {
       printf("Error: input images don't have same size\n");
       check = 0;
    }

    /* if the above checks are OK */
    if (check==1) {
      /* copy all the header keywords from first image to new output file */
      fits_read_key(afptr, TDOUBLE, "SATURATE", &sat, NULL, &status);
      if (status!=0) {
          status=0;
          sat=7000.;
      }
      fits_read_key(dfptr, TDOUBLE, "SATURATE", &sat0, NULL, &status);
      if (status!=0) {
          status=0;
          sat0=7000.;
      }

      fits_read_key(afptr, TDOUBLE, "SKYLEV", &sky, NULL, &status);
      if (status!=0) {
          status=0;
          sky=50.;
      }
      fits_read_key(dfptr, TDOUBLE, "SKYLEV", &sky0, NULL, &status);
      if (status!=0) {
          status=0;
          sky0=50.;
      }

      npixels = anaxes[0];  /* no. of pixels to read in each row */

      apix = (double *) malloc(npixels * sizeof(double)); /* mem for 1 row */
      bpix = (double *) malloc(npixels * sizeof(double));
      cpix = (double *) malloc(npixels * sizeof(double));
      dpix = (double *) malloc(npixels * sizeof(double));
      epix = (double *) malloc(npixels * sizeof(double));

      if (apix == NULL || bpix == NULL || cpix == NULL || dpix == NULL || epix == NULL) {
        printf("Memory allocation error\n");
        return(1);
      }
      if (argc>7) nsig=atof(argv[7]);
      if (argc>8) srcfac=atof(argv[8]);
      if (argc>9) s0=atof(argv[9]);

      s1 = sky0/sky;
      s12 = s1*s1;

      /* loop over all rows of the plane */
      for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
        lastpix[1] = anaxes[1] - firstpix[1] + 1;
        /* Read both images as doubles, regardless of actual datatype.  */
        /* Give starting pixel coordinate and no. of pixels to read.    */
        /* This version does not support undefined pixels in the image. */

        if (fits_read_pix(afptr, TDOUBLE, firstpix, npixels, NULL, apix,
                          NULL, &status)  ||
            fits_read_pix(bfptr, TDOUBLE, firstpix, npixels,  NULL, bpix,
                          NULL, &status)  ||
            fits_read_pix(cfptr, TDOUBLE, firstpix, npixels,  NULL, cpix,
                          NULL, &status)  ||
            fits_read_pix(dfptr, TDOUBLE, firstpix, npixels,  NULL, dpix,
                          NULL, &status)  ||
            fits_read_pix(efptr, TDOUBLE, firstpix, npixels,  NULL, epix,
                          NULL, &status)  )
            break;   /* jump out of loop on error */

        for(ii=0; ii< npixels; ii++) {
          if (bpix[ii]>0 && epix[ii]>0 && dpix[ii]<sat0 && apix[ii]<sat) {
              v00 = 1/bpix[ii];
              v11 = s12/epix[ii];
              vdiff = v00 + v11;
              if (apix[ii]>0) vdiff += v00*apix[ii]/s0;
              if (dpix[ii]>0) {
                  tmp = s1*dpix[ii];
                  vdiff += tmp*tmp + v11*dpix[ii]/s0;
              }
              diff = apix[ii]-s1*dpix[ii];
              nsig2=nsig*nsig;
              if (cpix[ii]==0) nsig2*=srcfac;
              if (diff*diff>nsig2*vdiff) apix[ii]=bpix[ii]=0;
          }
        }

        fits_write_pix(afptr, TDOUBLE, firstpix, npixels, apix, &status); /* write new values to output image */
        fits_write_pix(bfptr, TDOUBLE, firstpix, npixels, bpix, &status); /* write new values to output image */
      }

      free(apix);
      free(bpix);
      free(cpix);
      free(dpix);
      free(epix);
    }

    //bscale = sky/s0;
    //ibscale = 1/(bscale*bscale);
    //fits_update_key(afptr, TFLOAT, "BSCALE", &bscale, NULL, &status);
    //fits_update_key(bfptr, TFLOAT, "BSCALE", &ibscale, NULL, &status);
    //bscale *= sat;
    //fits_update_key(afptr, TFLOAT, "SATURATE", &bscale, NULL, &status);
    //fits_update_key(bfptr, TFLOAT, "SATURATE", &bscale, NULL, &status);

    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);
    fits_close_file(cfptr, &status);
    fits_close_file(dfptr, &status);
    fits_close_file(efptr, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
