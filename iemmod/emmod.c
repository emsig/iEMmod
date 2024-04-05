#include "emmod.h"

void emmod(complex *data, complex *ddata, REAL freq, int nx, REAL dx, int ny, REAL dy, int nrz, REAL *z, REAL *econdV, REAL *econdH, REAL *epermV, REAL *epermH, REAL *mpermV, REAL *mpermH, REAL zsrc, REAL zrcv, int *component, int nd, REAL kmax, REAL startlogx, REAL deltalogx, int nlogx, REAL c1, REAL c2, int maxpt, int xdirect, int dograd, int verbose)
{
        int     iz, fullspace, above, zrcv_layer, zsrc_layer;
        REAL    sumecondV, sumecondH, sumepermV, sumepermH, summpermV, summpermH;
        REAL    avecondV, avecondH, avepermV, avepermH, avmpermV, avmpermH;
        REAL    sumallpar;

        /* Check if medium is a homogeneous fullspace. If that is the case, 
           the EM-field is computed directly in the space-domain.
           Note: Also a stack of layers with the same material parameters is 
                 treated as a homogeneous fullspace. */
        sumecondV = 0.0;
        sumecondH = 0.0;
        sumepermV = 0.0;
        sumepermH = 0.0;
        summpermV = 0.0;
        summpermH = 0.0;
        for (iz=0; iz<nrz; iz++) {
                if (verbose) vmess("econdV[%d] = %30.20f",iz,econdV[iz]);
                if (verbose) vmess("econdH[%d] = %30.20f",iz,econdH[iz]);
                sumecondV = sumecondV + econdV[iz];
                sumecondH = sumecondH + econdH[iz];
                sumepermV = sumepermV + epermV[iz];
                sumepermH = sumepermH + epermH[iz];
                summpermV = summpermV + mpermV[iz];
                summpermH = summpermH + mpermH[iz];
        }
        avecondV = sumecondV/nrz;
        avecondH = sumecondH/nrz;
        avepermV = sumepermV/nrz;
        avepermH = sumepermH/nrz;
        avmpermV = summpermV/nrz;
        avmpermH = summpermH/nrz;
        sumecondV = 0.0;
        sumecondH = 0.0;
        sumepermV = 0.0;
        sumepermH = 0.0;
        summpermV = 0.0;
        summpermH = 0.0;
        for (iz=0; iz<nrz; iz++) {
                sumecondV = sumecondV+(econdV[iz]-avecondV)*(econdV[iz]-avecondV);
                sumecondH = sumecondH+(econdH[iz]-avecondH)*(econdH[iz]-avecondH);
                sumepermV = sumepermV+(epermV[iz]-avepermV)*(epermV[iz]-avepermV);
                sumepermH = sumepermH+(epermH[iz]-avepermH)*(epermH[iz]-avepermH);
                summpermV = summpermV+(mpermV[iz]-avmpermV)*(mpermV[iz]-avmpermV);
                summpermH = summpermH+(mpermH[iz]-avmpermH)*(mpermH[iz]-avmpermH);
        }
        sumallpar = sumecondV + sumecondH + sumepermV + sumepermH + summpermV + summpermH;
        if (sumallpar <= 10e-12) {
                fullspace = 1;
        } else {
                fullspace = 0;
        }

        /* Determine where the source level is in relation to the receiver level.
           Note: If zsrc or zrcv are on a layer interface, the layer above the interface
                 is chosen. */
        zrcv_layer = -1;
        zsrc_layer = -1;
        if (zrcv<=z[1]) zrcv_layer = 0;
        if (zsrc<=z[1]) zsrc_layer = 0;
        for (iz=2; iz<nrz; iz++) {
                if (zrcv<=z[iz] && zrcv>z[iz-1]) zrcv_layer = iz-1;
                if (zsrc<=z[iz] && zsrc>z[iz-1]) zsrc_layer = iz-1;
        }
        /* If zrcv_layer or zsrc_layer remain -1 after the loop, 
           the receivers or the source must be below all interfaces. */
        if (zrcv_layer == -1) zrcv_layer = nrz-1;
        if (zsrc_layer == -1) zsrc_layer = nrz-1;

        /* Determine if the receivers are above the source (above=1), 
           below the source (above=-1) or in the same layer as the source (above=0). */
        if (zrcv_layer == zsrc_layer) {
                above=0;
                // Receivers are in the same layer as the source
        } else if (zrcv_layer < zsrc_layer) {
                above=1;
                // Receivers are in a layer above the source
        } else {
                above=-1;
                // Receivers are in a layer below the source
        }

    if(above!=0) verr("The source and the receivers must be located in the same layer.");

        /* Call kxwmod with all the loaded parameters to do the actual modeling. */
        if (verbose) vmess("fullspace: %d",fullspace);
        kxwmod(data, ddata, freq, nx, dx, ny, dy, nrz, z, econdV, econdH, epermV, epermH, mpermV, mpermH, zsrc, zrcv, zrcv_layer, zsrc_layer, above, component, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, fullspace, xdirect, dograd, verbose);

        return;
}
