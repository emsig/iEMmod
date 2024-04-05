#include "emmod.h"

void kxwmod(complex *cdata, complex *dcdata, REAL freq, int nx, REAL dx, int ny, REAL dy, int nz, REAL *z, REAL *econdV, REAL *econdH, REAL *epermV, REAL *epermH, REAL *mpermV, REAL *mpermH, REAL zsrc, REAL zrcv, int zrcv_layer, int zsrc_layer, int above, int *component, int nd, REAL kmax, REAL startlogx, REAL deltalogx, int nlogx, REAL c1, REAL c2, int maxpt, int fullspace, int xdirect, int dograd, int verbose)
{    // This function computes the electromagnetic field for one source-receiver configuration in the wavenumber domain followed by a Hankel-Transformation to the space domain.
    // All source-receiver combinations are similar. Component 11 is commented with more details than the other components.
    int     *marker, *dmarker, *markernew, *dmarkernew, iz, ikx, iky, nlayers, temp, ind, ipos, corel, dnlogx, nlogxnew, dnlogxnew, dnlogxrdy, nlogxdata, dnlogxdata, newit, ixs, itcount, component1, component2, componentold, nocomp, nxh, nyh, directfield;
    REAL    om, eperm0, mperm0;
    REAL    kintel, a, b;
    REAL    *cor, *tempcor, *xpos, *dxpos, *xposnew, *dxposnew, *dxposrdy;
    complex *Gamma, *dGamma, *GammaB, *dGammaB, comptemp;
    complex *etaV, *etaH, *zetaV, *zetaH, *gamA, *gamB;
    complex *Rp, *dRp, *Rpbar, *dRpbar, *Rm, *dRm, *Rmbar, *dRmbar, *Pupplus, *dPupplus, *Pupplusbar, *dPupplusbar, *Pupmin, *dPupmin, *Pupminbar, *dPupminbar, *Pdownplus, *dPdownplus, *Pdownplusbar, *dPdownplusbar, *Pdownmin, *dPdownmin, *Pdownminbar, *dPdownminbar, *Wup, *dWup, *Wupbar, *dWupbar, *Wdown, *dWdown, *Wdownbar, *dWdownbar, *temptot, *dtemptot, *temptot1, *temptot2, *dtemptot1, *dtemptot2, *xtotrad1, *dxtotrad1, *xtotrad2, *dxtotrad2, *xtotrad1new, *dxtotrad1new, *xtotrad2new, *dxtotrad2new, *Ptot, *dPtot, *Gin;
    
    /* First only one quadrant is computed and subsequently copied to the other quadrants.
       The variables nxh and nyh contain half plus one element of the size of the complete setup. */
    nxh = nx/2+1;
    nyh = ny/2+1;
    
    /* The variable nlayers contains the amount of layers between the receiver and the source
       including the receiver and the source layer. */
    nlayers = abs(zsrc_layer - zrcv_layer) + 1;
    if (verbose) vmess("Number of layers between ind ncluding source and receiver layers (nlayers) = %d", nlayers);
    
    /* Initialize the variable nocomp.
       It will be set to 1 if the entered component does not exist. */
    nocomp = 0;
    
    /* Calculate eta, zeta and gamma for each layer. */
    etaV = (complex *)calloc(nz,sizeof(complex));
    etaH = (complex *)calloc(nz,sizeof(complex));
    zetaV = (complex *)calloc(nz,sizeof(complex));
    zetaH = (complex *)calloc(nz,sizeof(complex));
    gamA = (complex *)calloc(nz,sizeof(complex));
    gamB = (complex *)calloc(nz,sizeof(complex));
    
    om = freq*2.0*M_PI; // angular frequency
    eperm0 = 8.854187817e-12; // electric permittivity in free space
    mperm0 = 4.0e-7*M_PI; // magnetic permeability in free space
    
    /* Compute the complex material parameters eta and zeta in the horizontal and vertical direction for each layer. */
    for (iz=0; iz<nz; iz++) {
        etaH[iz].r = econdH[iz];
        etaH[iz].i = om*epermH[iz]*eperm0;
        etaV[iz].r = econdV[iz];
        etaV[iz].i = om*epermV[iz]*eperm0;
        zetaH[iz].r = 0.0;
        zetaH[iz].i = om*mpermH[iz]*mperm0;
        zetaV[iz].r = 0.0;
        zetaV[iz].i = om*mpermV[iz]*mperm0;
    }
    
    /* Determine on which points in the wavenumber domain the EM-field needs to be computed
     in order to integrate it with the 61-point Gauss-Kronrod rule (zqk61n.f). */
    corel = 61*nd;
    cor = (REAL *)calloc(corel,sizeof(REAL));
    tempcor = (REAL *)calloc(61,sizeof(REAL));
    kintel = kmax/nd;
    /* 61 points are necessary for each of the nd subdomains. */
    for (ind=0; ind<nd; ind++) {
        a = ind*kintel; // lower integration boundary
        b = (ind+1)*kintel; // upper integration boundary
        zqk61_setup_grid_(tempcor, &a, &b);
        for (ipos=0; ipos<61; ipos++) {
            cor[ind*61+ipos] = tempcor[ipos];
        }
    }
    
    /* The computation of the different components begins. */
    if (component[0] == 11) {
        /* This component has slightly more comments than the other components.
         It serves as an example, because all other components work similarly. */
        if (verbose) vmess("Computing now component %d", component[0]);
        
        /* Allocate memory */
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wdownbar = (complex *)calloc(corel,sizeof(complex));
        dWdownbar = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Wupbar = (complex *)calloc(corel,sizeof(complex));
        dWupbar = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rpbar = (complex *)calloc(corel*nz,sizeof(complex));
        dRpbar = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rmbar = (complex *)calloc(corel*nz,sizeof(complex));
        dRmbar = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupmin = (complex *)calloc(corel,sizeof(complex));
        dPupmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownmin = (complex *)calloc(corel,sizeof(complex));
        dPdownmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pupplusbar = (complex *)calloc(corel,sizeof(complex));
        dPupplusbar = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
        dPdownplusbar = (complex *)calloc(corel*nz,sizeof(complex));
        temptot1 = (complex *)calloc(corel,sizeof(complex));
        dtemptot1 = (complex *)calloc(corel*nz,sizeof(complex));
        temptot2 = (complex *)calloc(corel,sizeof(complex));
        dtemptot2 = (complex *)calloc(corel*nz,sizeof(complex));
        // The following six arrays are allocated for 100'000 elements, because they change size during running.
        // Deallocation followed by reallocation is not efficient and, therefore, avoided.
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        xtotrad2 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad2 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        GammaB = (complex *)calloc(nz*corel,sizeof(complex));
        dGammaB = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
            dgammarad_(&GammaB[iz*corel], &dGammaB[iz*corel], &gamB[iz], &corel, cor, &zetaV[iz], &zetaH[iz], &etaH[iz], &zetaH[iz]);
        }
        
        /* Wavenumber domain calculations */
        /* compute field propagaters and its derivatives in the layer where the source and the receivers are located */
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        dwprop_(Wupbar, Wdownbar, dWupbar, dWdownbar, GammaB, dGammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        /* compute the reflection coming from below the receiver R+ => Pu and its derivative */
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        rplusderiv_(Rpbar, dRpbar, &corel, cor, zetaH, GammaB, dGammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        /* compute the reflection coming from above the receiver R- => Pd and its derivative */
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        rminderiv_(Rmbar, dRmbar, &corel, cor, zetaH, GammaB, dGammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        /* compute the field at receiver level coming from below the receiver Pu and its derivative */
        pupminderiv_(Pupmin, dPupmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        /* compute the field at receiver level coming from above the receiver Pd and its derivative */
        pdownminderiv_(Pdownmin, dPdownmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        /* compute the field at receiver level coming from below the receiver Pu and its derivative */
        pupplusderiv_(Pupplusbar, dPupplusbar, Rpbar, Rmbar, dRpbar, dRmbar, GammaB, dGammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        /* compute the field at receiver level coming from above the receiver Pd and its derivative */
        pdownplusderiv_(Pdownplusbar, dPdownplusbar, Rpbar, Rmbar, dRpbar, dRmbar, GammaB, dGammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        /* compute total field and its derivative */
        dptotalxx_(temptot1, temptot2, dtemptot1, dtemptot2, Pdownmin, Pupmin, Pdownplusbar, Pupplusbar, dPdownmin, dPupmin, dPdownplusbar, dPupplusbar, Wup, Wdown, Wupbar, Wdownbar, dWup, dWdown, dWupbar, dWdownbar, Gamma, GammaB, dGamma, dGammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, gamB, &nlayers, &above, &xdirect);
        
        /* Transform from the wavenumber to the space domain */
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine
         no more points are added to the coordinate vector
         and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        
        /* In the while-loop the points in space for which the Hankel-transformation has to be carried out
         are determined. The loop runs until no more points need to be added for the specified precision
         (newit=0) or until the maximum amount of points (maxpt) has been reached. */
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        xtotrad2new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform for points that the marker is 2 */
            hankeltransxx_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
            
            /* The variable nlogxnew has the size of the maximum amount of points after this iteration
             is completed. That is nlogx points from the previous iteration plus nlogx-1 points that
             could be added (In the worst case, between every point a new datapoint is required.). */
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data */
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed */
            evalpoints_lin_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    xtotrad2[ixs].r = xtotrad2new[ikx].r;
                    xtotrad2[ixs].i = xtotrad2new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                // More points than the maximum amount of points are necessary to achieve the required precision.
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                // The Hankel Transform is still computed for the points selected in this iteration, before the while-loop is exited.
                hankeltransxx_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        free(xtotrad2new);
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // The field is now known on a irregularly sampled line in the space domain. The call to the gridit-function uses symmetry to compute the field on a 2D-grid.
        gridit_xx_lin_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
        
        free(marker);
        free(xpos);
        free(xtotrad1);
        free(xtotrad2);
        
        // The same is now also done for the gradient
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { // The gradient of the field is computed
            /* To get the field with respect to all the layers, the Hankel Transform needs to be done nz times. */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            dxtotrad2new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity with respect to layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                /* In the while-loop the points in space for which the Hankel-transformation has to be carried out
                 are determined. The loop runs until no more points need to be added for the specified precision
                 (newit=0) or until the maximum amount of points (maxpt) has been reached. */
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltransxx_(dxtotrad1,dxtotrad2,dmarker,&dtemptot1[iz*corel],&dtemptot2[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    
                    /* The variable nlogxnew has the size of the maximum amount of points after this iteration
                     is completed. That is nlogx points from the previous iteration plus nlogx-1 points that
                     could be added (In the worst case, between every point a new datapoint is required.). */
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data */
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed */
                    evalpoints_lin_(dxposnew,dxtotrad1new,dxtotrad2new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,dxtotrad2,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            dxtotrad2[ixs].r = dxtotrad2new[ikx].r;
                            dxtotrad2[ixs].i = dxtotrad2new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltransxx_(dxtotrad1,dxtotrad2,dmarker,&dtemptot1[iz*corel],&dtemptot2[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_xx_lin_(&dPtot[iz*nxh*nyh], dxtotrad1, dxtotrad2, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect, &directfield);
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
            free(dxtotrad2new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        free(dxtotrad2);
        
        /* Uncomment this piece of code to write intermediate results to the harddisk. */
        /*      {
         int nwrite;
         FILE *fp_out = fopen("dtemptot2_.bin","w");
         nwrite = fwrite( &dtemptot2[0].r, sizeof(complex), corel, fp_out);
         assert(nwrite == corel);
         fclose(fp_out);
         }
         {
         int nwrite;
         FILE *fp_out = fopen("dtemptot2_.bin","w");
         nwrite = fwrite( &dtemptot2[0].r, sizeof(complex), corel, fp_out);
         assert(nwrite == corel);
         fclose(fp_out);
         }*/
        
        // The field has been only computed for one quadrant.
        // Here we copy quadrants to other parts to build up the whole field.
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        /* Free memory */
        if (fullspace == 0) {
            free(Wdown);
            free(Wdownbar);
            free(Wup);
            free(Wupbar);
            free(Rp);
            free(Rpbar);
            free(Rm);
            free(Rmbar);
            free(Pupmin);
            free(Pdownmin);
            free(Pupplusbar);
            free(Pdownplusbar);
            free(temptot1);
            free(temptot2);
            free(dWdown);
            free(dWdownbar);
            free(dWup);
            free(dWupbar);
            free(dRp);
            free(dRpbar);
            free(dRm);
            free(dRmbar);
            free(dPupmin);
            free(dPdownmin);
            free(dPupplusbar);
            free(dPdownplusbar);
            free(dtemptot1);
            free(dtemptot2);
        }
    }else if (component[0] == 12 || component[0] == 21) {
        if (verbose) vmess("Computing now component %d", component[0]);
        
        // Allocate memory
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wdownbar = (complex *)calloc(corel,sizeof(complex));
        dWdownbar = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Wupbar = (complex *)calloc(corel,sizeof(complex));
        dWupbar = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rpbar = (complex *)calloc(corel*nz,sizeof(complex));
        dRpbar = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rmbar = (complex *)calloc(corel*nz,sizeof(complex));
        dRmbar = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupmin = (complex *)calloc(corel,sizeof(complex));
        dPupmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownmin = (complex *)calloc(corel,sizeof(complex));
        dPdownmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pupplusbar = (complex *)calloc(corel,sizeof(complex));
        dPupplusbar = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
        dPdownplusbar = (complex *)calloc(corel*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        GammaB = (complex *)calloc(nz*corel,sizeof(complex));
        dGammaB = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
            dgammarad_(&GammaB[iz*corel], &dGammaB[iz*corel], &gamB[iz], &corel, cor, &zetaV[iz], &zetaH[iz], &etaH[iz], &zetaH[iz]);
        }
        
        /* compute field propagaters in the layer where the source and the receivers are located */
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        dwprop_(Wupbar, Wdownbar, dWupbar, dWdownbar, GammaB, dGammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        /* compute the reflection coming from below the receiver R+ => Pu */
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        rplusderiv_(Rpbar, dRpbar, &corel, cor, zetaH, GammaB, dGammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        /* compute the reflection coming from above the receiver R- => Pd */
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        rminderiv_(Rmbar, dRmbar, &corel, cor, zetaH, GammaB, dGammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        /* compute the field at receiver level coming below the receiver Pu */
        pupminderiv_(Pupmin, dPupmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        /* compute the field at receiver level coming above the receiver Pd */
        pdownminderiv_(Pdownmin, dPdownmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        /* compute the field at receiver level coming below the receiver Pu */
        pupplusderiv_(Pupplusbar, dPupplusbar, Rpbar, Rmbar, dRpbar, dRmbar, GammaB, dGammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        /* compute the field at receiver level coming above the receiver Pd */
        pdownplusderiv_(Pdownplusbar, dPdownplusbar, Rpbar, Rmbar, dRpbar, dRmbar, GammaB, dGammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        
        // compute total field
        dptotalxy_(temptot, dtemptot, Pdownmin, Pupmin, Pdownplusbar, Pupplusbar, dPdownmin, dPupmin, dPdownplusbar, dPupplusbar, Wup, Wdown, Wupbar, Wdownbar, dWup, dWdown, dWupbar, dWdownbar, Gamma, GammaB, dGamma, dGammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, gamB, &nlayers, &above, &xdirect);
        
        /* Set pre-determined points equal for the derivatives*/
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        // getcoords_(dxposrdy,&startlogx,&deltalogx,&dnlogxrdy);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltransxy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltransxy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_xy_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { //The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                /* In the while-loop the points in space for which the Hankel-transformation has to be carried out
                 are determined. The loop runs until no more points need to be added for the specified precision
                 (newit=0) or until the maximum amount of points (maxpt) has been reached. */
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltransxy_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    
                    /* The variable nlogxnew has the size of the maximum amount of points after this iteration
                     is completed. That is nlogx points from the previous iteration plus nlogx-1 points that
                     could be added (In the worst case, between every point a new datapoint is required.). */
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltransxy_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_xy_lin_(&dPtot[iz*nxh*nyh], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect, &directfield);
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = -dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = -dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = -dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = -dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        if (fullspace == 0) {
            free(Wdown);
            free(Wdownbar);
            free(Wup);
            free(Wupbar);
            free(Rp);
            free(Rpbar);
            free(Rm);
            free(Rmbar);
            free(Pupmin);
            free(Pdownmin);
            free(Pupplusbar);
            free(Pdownplusbar);
            free(temptot);
            free(dWdown);
            free(dWdownbar);
            free(dWup);
            free(dWupbar);
            free(dRp);
            free(dRpbar);
            free(dRm);
            free(dRmbar);
            free(dPupmin);
            free(dPdownmin);
            free(dPupplusbar);
            free(dPdownplusbar);
            free(dtemptot);
        }
    }else if (component[0] == 13) {
        if (verbose) vmess("Computing now component %d", component[0]);
        
        // Allocate memory
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupplus = (complex *)calloc(corel,sizeof(complex));
        dPupplus = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownplus = (complex *)calloc(corel,sizeof(complex));
        dPdownplus = (complex *)calloc(corel*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the reflection coming from above the receiver R- => Pd
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the field at receiver level coming below the receiver Pu
        pupplusderiv_(Pupplus, dPupplus, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming above the receiver Pd
        pdownplusderiv_(Pdownplus, dPdownplus, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute total field: incident + Pdownmin + Pupmin
        dptotalxz_(temptot, dtemptot, Pdownplus, Pupplus, dPdownplus, dPupplus, Wup, Wdown, dWup, dWdown, Gamma, dGamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);
        
        /* Set pre-determined points equal for the derivatives*/
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltransxz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltransxz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_xz_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { //The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                /* In the while-loop the points in space for which the Hankel-transformation has to be carried out
                 are determined. The loop runs until no more points need to be added for the specified precision
                 (newit=0) or until the maximum amount of points (maxpt) has been reached. */
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltransxz_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    
                    /* The variable nlogxnew has the size of the maximum amount of points after this iteration
                     is completed. That is nlogx points from the previous iteration plus nlogx-1 points that
                     could be added (In the worst case, between every point a new datapoint is required.). */
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltransxz_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_xz_lin_(&dPtot[iz*nxh*nyh], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect, &directfield);
            }// end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        }// end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = -dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = -dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        if (fullspace == 0) {
            free(Wdown);
            free(Wup);
            free(Rp);
            free(Rm);
            free(Pupplus);
            free(Pdownplus);
            free(temptot);
            free(dWdown);
            free(dWup);
            free(dRp);
            free(dRm);
            free(dPupplus);
            free(dPdownplus);
            free(dtemptot);
        }
    }else if (component[0] == 22) {
        if (verbose) vmess("Computing now component %d", component[0]);
        
        // Allocate memory
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wdownbar = (complex *)calloc(corel,sizeof(complex));
        dWdownbar = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Wupbar = (complex *)calloc(corel,sizeof(complex));
        dWupbar = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rpbar = (complex *)calloc(corel*nz,sizeof(complex));
        dRpbar = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rmbar = (complex *)calloc(corel*nz,sizeof(complex));
        dRmbar = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupmin = (complex *)calloc(corel,sizeof(complex));
        dPupmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownmin = (complex *)calloc(corel,sizeof(complex));
        dPdownmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pupplusbar = (complex *)calloc(corel,sizeof(complex));
        dPupplusbar = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
        dPdownplusbar = (complex *)calloc(corel*nz,sizeof(complex));
        temptot1 = (complex *)calloc(corel,sizeof(complex));
        dtemptot1 = (complex *)calloc(corel*nz,sizeof(complex));
        temptot2 = (complex *)calloc(corel,sizeof(complex));
        dtemptot2 = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        xtotrad2 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad2 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        GammaB = (complex *)calloc(nz*corel,sizeof(complex));
        dGammaB = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
            dgammarad_(&GammaB[iz*corel], &dGammaB[iz*corel], &gamB[iz], &corel, cor, &zetaV[iz], &zetaH[iz], &etaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        dwprop_(Wupbar, Wdownbar, dWupbar, dWdownbar, GammaB, dGammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        rplusderiv_(Rpbar, dRpbar, &corel, cor, zetaH, GammaB, dGammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the reflection coming from above the receiver R- => Pd
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        rminderiv_(Rmbar, dRmbar, &corel, cor, zetaH, GammaB, dGammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the field at receiver level coming below the receiver Pu
        pupminderiv_(Pupmin, dPupmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming above the receiver Pd
        pdownminderiv_(Pdownmin, dPdownmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming below the receiver Pu
        pupplusderiv_(Pupplusbar, dPupplusbar, Rpbar, Rmbar, dRpbar, dRmbar, GammaB, dGammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming above the receiver Pd
        pdownplusderiv_(Pdownplusbar, dPdownplusbar, Rpbar, Rmbar, dRpbar, dRmbar, GammaB, dGammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute total field
        dptotalyy_(temptot1, temptot2, dtemptot1, dtemptot2, Pdownmin, Pupmin, Pdownplusbar, Pupplusbar, dPdownmin, dPupmin, dPdownplusbar, dPupplusbar, Wup, Wdown, Wupbar, Wdownbar, dWup, dWdown, dWupbar, dWdownbar, Gamma, GammaB, dGamma, dGammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, gamB, &nlayers, &above, &xdirect);
        
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        xtotrad2new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltransyy_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    xtotrad2[ixs].r = xtotrad2new[ikx].r;
                    xtotrad2[ixs].i = xtotrad2new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltransyy_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        free(xtotrad2new);
        
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_yy_lin_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
        
        free(marker);
        free(xpos);
        free(xtotrad1);
        free(xtotrad2);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { //The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            dxtotrad2new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                /* In the while-loop the points in space for which the Hankel-transformation has to be carried out
                 are determined. The loop runs until no more points need to be added for the specified precision
                 (newit=0) or until the maximum amount of points (maxpt) has been reached. */
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltransyy_(dxtotrad1,dxtotrad2,dmarker,&dtemptot1[iz*corel],&dtemptot2[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    
                    /* The variable nlogxnew has the size of the maximum amount of points after this iteration
                     is completed. That is nlogx points from the previous iteration plus nlogx-1 points that
                     could be added (In the worst case, between every point a new datapoint is required.). */
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_(dxposnew,dxtotrad1new,dxtotrad2new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,dxtotrad2,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            dxtotrad2[ixs].r = dxtotrad2new[ikx].r;
                            dxtotrad2[ixs].i = dxtotrad2new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltransyy_(dxtotrad1,dxtotrad2,dmarker,&dtemptot1[iz*corel],&dtemptot2[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_yy_lin_(&dPtot[iz*nxh*nyh], dxtotrad1, dxtotrad2, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect, &directfield);
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
            free(dxtotrad2new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        free(dxtotrad2);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        if (fullspace == 0) {
            free(Wdown);
            free(Wdownbar);
            free(Wup);
            free(Wupbar);
            free(Rp);
            free(Rpbar);
            free(Rm);
            free(Rmbar);
            free(Pupmin);
            free(Pdownmin);
            free(Pupplusbar);
            free(Pdownplusbar);
            free(temptot1);
            free(temptot2);
            free(dWdown);
            free(dWdownbar);
            free(dWup);
            free(dWupbar);
            free(dRp);
            free(dRpbar);
            free(dRm);
            free(dRmbar);
            free(dPupmin);
            free(dPdownmin);
            free(dPupplusbar);
            free(dPdownplusbar);
            free(dtemptot1);
            free(dtemptot2);
        }
    }else if (component[0] == 23) {
        if (verbose) vmess("Computing now component %d", component[0]);
        
        // Allocate memory
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupplus = (complex *)calloc(corel,sizeof(complex));
        dPupplus = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownplus = (complex *)calloc(corel,sizeof(complex));
        dPdownplus = (complex *)calloc(corel*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the reflection coming from above the receiver R- => Pd
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the field at receiver level coming below the receiver Pu
        pupplusderiv_(Pupplus, dPupplus, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming above the receiver Pd
        pdownplusderiv_(Pdownplus, dPdownplus, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute total field: incident + Pdownmin + Pupmin
        dptotalyz_(temptot, dtemptot, Pdownplus, Pupplus, dPdownplus, dPupplus, Wup, Wdown, dWup, dWdown, Gamma, dGamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);
        
        /* Set pre-determined points equal for the derivatives*/
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltransyz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltransyz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_yz_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
        
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { //The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                /* In the while-loop the points in space for which the Hankel-transformation has to be carried out
                 are determined. The loop runs until no more points need to be added for the specified precision
                 (newit=0) or until the maximum amount of points (maxpt) has been reached. */
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltransyz_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    
                    /* The variable nlogxnew has the size of the maximum amount of points after this iteration
                     is completed. That is nlogx points from the previous iteration plus nlogx-1 points that
                     could be added (In the worst case, between every point a new datapoint is required.). */
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltransyz_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_yz_lin_(&dPtot[iz*nxh*nyh], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect, &directfield);
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = -dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = -dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        if (fullspace == 0) {
            free(Wdown);
            free(Wup);
            free(Rp);
            free(Rm);
            free(Pupplus);
            free(Pdownplus);
            free(temptot);
            free(dWdown);
            free(dWup);
            free(dRp);
            free(dRm);
            free(dPupplus);
            free(dPdownplus);
            free(dtemptot);
        }
    }else if (component[0] == 31) {
        if (verbose) vmess("Computing now component %d", component[0]);
        
        // Allocate memory
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupmin = (complex *)calloc(corel,sizeof(complex));
        dPupmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownmin = (complex *)calloc(corel,sizeof(complex));
        dPdownmin = (complex *)calloc(corel*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the reflection coming from above the receiver R- => Pd
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the field at receiver level coming below the receiver Pu
        pupminderiv_(Pupmin, dPupmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming above the receiver Pd
        pdownminderiv_(Pdownmin, dPdownmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute total field: incident + Pdownmin + Pupmin
        dptotalzx_(temptot, dtemptot, Pdownmin, dPdownmin, Pupmin, dPupmin, Wup, dWup, Wdown, dWdown, Gamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);
        
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltranszx_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltranszx_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        gridit_zx_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
        
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { //The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltranszx_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    /* In the worst case, between every point a new datapoint is required.*/
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltranszx_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_zx_lin_(&dPtot[iz*(nxh*nyh)], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect, &directfield);
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = -dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = -dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        if (fullspace == 0) {
            free(Wdown);
            free(dWdown);
            free(Wup);
            free(dWup);
            free(Rp);
            free(dRp);
            free(Rm);
            free(dRm);
            free(Pupmin);
            free(dPupmin);
            free(Pdownmin);
            free(dPdownmin);
            free(temptot);
            free(dtemptot);
        }
    }else if (component[0] == 32) {
        if (verbose) vmess("Computing now component %d", component[0]);
        
        // Allocate memory
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupmin = (complex *)calloc(corel,sizeof(complex));
        dPupmin = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownmin = (complex *)calloc(corel,sizeof(complex));
        dPdownmin = (complex *)calloc(corel*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the reflection coming from above the receiver R- => Pd
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the field at receiver level coming below the receiver Pu
        pupminderiv_(Pupmin, dPupmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming above the receiver Pd
        pdownminderiv_(Pdownmin, dPdownmin, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute total field: incident + Pdownmin + Pupmin
        dptotalzy_(temptot, dtemptot, Pdownmin, dPdownmin, Pupmin, dPupmin, Wup, dWup, Wdown, dWdown, Gamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);
        
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltranszy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltranszy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_zy_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
        
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { //The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltranszy_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    /* In the worst case, between every point a new datapoint is required.*/
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltranszy_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_zy_lin_(&dPtot[iz*(nxh*nyh)], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect, &directfield);
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = -dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = -dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        if (fullspace == 0) {
            free(Wdown);
            free(dWdown);
            free(Wup);
            free(dWup);
            free(Rp);
            free(dRp);
            free(Rm);
            free(dRm);
            free(Pupmin);
            free(dPupmin);
            free(Pdownmin);
            free(dPdownmin);
            free(temptot);
            free(dtemptot);
        }
    }else if (component[0] == 33) {
        if (verbose) vmess("Computing now component %d", component[0]);
        
        // Allocate memory
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Rm = (complex *)calloc(corel*nz,sizeof(complex));
        dRm = (complex *)calloc(corel*nz*nz,sizeof(complex));
        Pupplus = (complex *)calloc(corel,sizeof(complex));
        dPupplus = (complex *)calloc(corel*nz,sizeof(complex));
        Pdownplus = (complex *)calloc(corel,sizeof(complex));
        dPdownplus = (complex *)calloc(corel*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the reflection coming from above the receiver R- => Pd
        rminderiv_(Rm, dRm, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above, zetaH);
        // compute the field at receiver level coming below the receiver Pu
        pupplusderiv_(Pupplus, dPupplus, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute the field at receiver level coming above the receiver Pd
        pdownplusderiv_(Pdownplus, dPdownplus, Rp, Rm, dRp, dRm, Gamma, dGamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
        // compute total field: incident + Pdown + Pup
        dptotalzz_(temptot, dtemptot, Pdownplus, dPdownplus, Pupplus, dPupplus, Wup, dWup, Wdown, dWdown, Gamma, dGamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);
        
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltranszz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltranszz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            }
            
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_zz_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { //The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0;
                }
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltranszz_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    /* In the worst case, between every point a new datapoint is required.*/
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltranszz_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                dgridit_zz_lin_(&dPtot[iz*(nxh*nyh)], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect, &directfield);
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        if (dograd==1) { //The gradient of the field is computed
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        if (fullspace == 0) {
            free(Wdown);
            free(Wup);
            free(Rp);
            free(Rm);
            free(Pupplus);
            free(Pdownplus);
            free(temptot);
            free(dWdown);
            free(dWup);
            free(dRp);
            free(dRm);
            free(dPupplus);
            free(dPdownplus);
            free(dtemptot);
        }
    }else if (component[0] == 77) {
        if (verbose) vmess("Computing now component %d: TM reflection response", component[0]);
        // Allocate memory
        Wup = (complex *)calloc(corel,sizeof(complex));
        dWup = (complex *)calloc(corel,sizeof(complex));
        Wdown = (complex *)calloc(corel,sizeof(complex));
        dWdown = (complex *)calloc(corel,sizeof(complex));
        Rp = (complex *)calloc(corel*nz,sizeof(complex));
        dRp = (complex *)calloc(corel*nz*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        Gamma = (complex *)calloc(nz*corel,sizeof(complex));
        dGamma = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&Gamma[iz*corel], &dGamma[iz*corel], &gamA[iz], &corel, cor, &etaV[iz], &etaH[iz], &zetaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wup, Wdown, dWup, dWdown, Gamma, dGamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        temp = 0;
        rplusderiv_(Rp, dRp, &corel, cor, etaH, Gamma, dGamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &temp, zetaH);
        // compute total field
        dptotalref_(temptot, dtemptot, Wup, dWup, Rp, dRp, &corel, &nz, &zrcv_layer);
        
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            } else {
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_ref_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx);
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { // The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0; }
                //free(dxpos);
                //free(dxtotrad1);
                //free(dxtotrad2);
                //if (iz>0) {
                //        free(dmarker);
                //}
                
                //dxpos = (REAL *)calloc(dnlogxrdy,sizeof(REAL));
                //dxtotrad1 = (complex *)calloc(dnlogxrdy,sizeof(complex));
                //dxtotrad2 = (complex *)calloc(dnlogxrdy,sizeof(complex));
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltransref_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    /* In the worst case, between every point a new datapoint is required.*/
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltransref_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    } else {
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                gridit_ref_lin_(&dPtot[iz*nxh*nyh], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx);
                
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        
        // copy quadrants to other parts
        if (dograd==1) {
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        free(Wup);
        free(Wdown);
        free(Rp);
        free(dWup);
        free(dWdown);
        free(dRp);
        free(temptot);
        free(dtemptot);
    }else if (component[0] == 88) {
        if (verbose) vmess("Computing now component %d: TE reflection response", component[0]);
        // Allocate memory
        Wupbar = (complex *)calloc(corel,sizeof(complex));
        dWupbar = (complex *)calloc(corel,sizeof(complex));
        Wdownbar = (complex *)calloc(corel,sizeof(complex));
        dWdownbar = (complex *)calloc(corel,sizeof(complex));
        Rpbar = (complex *)calloc(corel*nz,sizeof(complex));
        dRpbar = (complex *)calloc(corel*nz*nz,sizeof(complex));
        temptot = (complex *)calloc(corel,sizeof(complex));
        dtemptot = (complex *)calloc(corel*nz,sizeof(complex));
        xpos = (REAL *)calloc(100000,sizeof(REAL));
        dxpos = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1 = (complex *)calloc(100000,sizeof(complex));
        dxtotrad1 = (complex *)calloc(100000,sizeof(complex));
        Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
        dPtot = (complex *)calloc(nxh*nyh*nz,sizeof(complex));
        
        GammaB = (complex *)calloc(nz*corel,sizeof(complex));
        dGammaB = (complex *)calloc(nz*corel,sizeof(complex));
        
        for (iz=0; iz<nz; iz++) {
            dgammarad_(&GammaB[iz*corel], &dGammaB[iz*corel], &gamB[iz], &corel, cor, &zetaV[iz], &zetaH[iz], &etaH[iz], &zetaH[iz]);
        }
        
        // compute field propagaters in the layer where the source and the receivers are located
        dwprop_(Wupbar, Wdownbar, dWupbar, dWdownbar, GammaB, dGammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
        // compute the reflection coming from below the receiver R+ => Pu
        temp = 0;
        rplusderiv_(Rpbar, dRpbar, &corel, cor, zetaH, GammaB, dGammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &temp, zetaH);
        // compute total field
        dptotalref_(temptot, dtemptot, Wupbar, dWupbar, Rpbar, dRpbar, &corel, &nz, &zrcv_layer);
        
        dnlogxrdy = nlogx;
        
        /* Set up coordinate vector*/
        getcoords_(xpos,&startlogx,&deltalogx,&nlogx);
        
        /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
        newit = 1;
        /* The Hankel transformation will be evaluated where the marker is 2.*/
        marker = (int *)calloc(100000,sizeof(int));
        for (ikx=0; ikx<nlogx; ikx++) {
            marker[ikx] = 2;
        }
        itcount = 0;
        markernew = (int *)calloc(100000,sizeof(int));
        xposnew = (REAL *)calloc(100000,sizeof(REAL));
        xtotrad1new = (complex *)calloc(100000,sizeof(complex));
        while ((newit == 1) && (nlogx < maxpt)) {
            if (verbose) vmess("Iteration %d", itcount);
            /* Compute the Hankel Transform */
            hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            /* In the worst case, between every point a new datapoint is required.*/
            nlogxnew = 2*nlogx-1;
            /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
            nlogxdata = 0;
            /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
            evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
            if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
            /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
             The marker is set to 2 at these locations.*/
            nlogx = nlogxdata;
            ixs = 0;
            for (ikx=0; ikx<nlogxnew; ikx++) {
                if (markernew[ikx]!=0) {
                    marker[ixs] = markernew[ikx];
                    xpos[ixs] = xposnew[ikx];
                    xtotrad1[ixs].r = xtotrad1new[ikx].r;
                    xtotrad1[ixs].i = xtotrad1new[ikx].i;
                    ixs++;
                }
            }
            if (nlogx>maxpt) {
                if (verbose) vwarn("More points are needed to achieve required precision,");
                if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
            } else {
            }
            itcount = itcount + 1;
        }
        free(markernew);
        free(xposnew);
        free(xtotrad1new);
        
        if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);
        
        // Compute the field for one quadrant based on the radial data
        gridit_ref_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx);
        free(marker);
        free(xpos);
        free(xtotrad1);
        
        dmarker = (int *)calloc(100000,sizeof(int));
        if (dograd==1) { // The gradient of the field is computed
            /* To get the field wrt to all the layers the hankeltransform needs to be done nz times  */
            dmarkernew = (int *)calloc(100000,sizeof(int));
            dxposnew = (REAL *)calloc(100000,sizeof(REAL));
            dxtotrad1new = (complex *)calloc(100000,sizeof(complex));
            for (iz=0; iz<nz; iz++) {
                if (verbose) vmess("Conductivity wrt layer %d of %d layers", iz+1, nz);
                if (iz==zsrc_layer) {
                    directfield = 1;
                } else {
                    directfield = 0; }
                //free(dxpos);
                //free(dxtotrad1);
                //free(dxtotrad2);
                //if (iz>0) {
                //        free(dmarker);
                //}
                
                //dxpos = (REAL *)calloc(dnlogxrdy,sizeof(REAL));
                //dxtotrad1 = (complex *)calloc(dnlogxrdy,sizeof(complex));
                //dxtotrad2 = (complex *)calloc(dnlogxrdy,sizeof(complex));
                
                getcoords_(dxpos,&startlogx,&deltalogx,&dnlogxrdy);
                dnlogx = dnlogxrdy;
                
                /* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
                newit = 1;
                /* The Hankel transformation will be evaluated where the marker is 2.*/
                for (ikx=0; ikx<dnlogx; ikx++) {
                    dmarker[ikx] = 2;
                }
                itcount = 0;
                while ((newit == 1) && (dnlogx < maxpt)) {
                    if (verbose) vmess("Derivative iteration %d", itcount);
                    /* Compute the Hankel Transform */
                    hankeltransref_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    /* In the worst case, between every point a new datapoint is required.*/
                    dnlogxnew = 2*dnlogx-1;
                    /* nlogxdata will show how many of the nlogxnew elements actually contain data*/
                    dnlogxdata = 0;
                    /* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
                    evalpoints_lin_mono_(dxposnew,dxtotrad1new,dmarkernew,&dnlogxnew,&dnlogxdata,&newit,dxpos,dxtotrad1,&dnlogx,&c1,&c2);
                    if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",dnlogxdata-dnlogx,dnlogxnew-dnlogx);
                    /* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
                     The marker is set to 2 at these locations.*/
                    dnlogx = dnlogxdata;
                    ixs = 0;
                    for (ikx=0; ikx<dnlogxnew; ikx++) {
                        if (dmarkernew[ikx]!=0) {
                            dmarker[ixs] = dmarkernew[ikx];
                            dxpos[ixs] = dxposnew[ikx];
                            dxtotrad1[ixs].r = dxtotrad1new[ikx].r;
                            dxtotrad1[ixs].i = dxtotrad1new[ikx].i;
                            ixs++;
                        }
                    }
                    if (dnlogx>maxpt) {
                        if (verbose) vwarn("More points are needed to achieve required precision,");
                        if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
                        hankeltransref_(dxtotrad1,dmarker,&dtemptot[iz*corel],&corel,cor,dxpos,&dnlogx,&nd,&kmax);
                    } else {
                    }
                    itcount = itcount + 1;
                }
                if (verbose) vmess("Using %d datapoints in the space domain.",dnlogx);
                
                // Compute the field for one quadrant based on the radial data
                gridit_ref_lin_(&dPtot[iz*nxh*nyh], dxtotrad1, &dx, &nxh, &dy, &nyh, dxpos, &dnlogx);
                
            } // end for
            free(dmarkernew);
            free(dxposnew);
            free(dxtotrad1new);
        } // end dograd
        free(dmarker);
        free(dxpos);
        free(dxtotrad1);
        
        // copy quadrants to other parts
        for (iky = 0; iky<nyh; iky++) {
            // 1st quadrant
            for (ikx = 0; ikx<nxh; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
            }
            // 2nd quadrant
            for (ikx = nxh; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
                cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
            }
        }
        for (iky = nyh; iky<ny; iky++) {
            // 3rd and 4th quadrant
            for (ikx = 0; ikx<nx; ikx++) {
                cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
                cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
            }
        }
        
        // copy quadrants to other parts
        if (dograd==1) {
            for (iz = 0; iz<nz; iz++) {
                for (iky = 0; iky<nyh; iky++) {
                    // 1st quadrant
                    for (ikx = 0; ikx<nxh; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((nxh-ikx-1)+(nyh-iky-1)*nxh)].i;
                    }
                    // 2nd quadrant
                    for (ikx = nxh; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dPtot[(iz*(nxh*nyh))+((ikx-nxh+1)+(nyh-iky-1)*nxh)].i;
                    }
                }
                for (iky = nyh; iky<ny; iky++) {
                    // 3rd and 4th quadrant
                    for (ikx = 0; ikx<nx; ikx++) {
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].r = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].r;
                        dcdata[(iz*(nx*ny))+ikx+iky*nx].i = dcdata[(iz*(nx*ny))+ikx+(ny-iky)*nx].i;
                    }
                }
            }
        }
        
        // Free memory
        free(Wupbar);
        free(dWupbar);
        free(Wdownbar);
        free(dWdownbar);
        free(Rpbar);
        free(dRpbar);
        free(temptot);
        free(dtemptot);
    }else { // in case the user entered a component that does not exist
        if (verbose) vmess("The component %d is not implemented in the current version of emmod.",component[0]);
        nocomp = 1;
    }
    
    /* Free memory */
    if (nocomp != 1) { // if the component entered exists
        if (fullspace == 0) { // if the model is layered
            if (component[0] == 31 || component[0] == 13 || component[0] == 32 || component[0] == 33 || component[0] == 23 || component[0] == 43 || component[0] == 53 || component[0] == 77) { // components containing only the TM-mode
                free(Gamma);
                free(dGamma);
                free(Ptot);
                free(dPtot);
            }else if (component[0] == 61 || component[0] == 62 || component[0] == 88) { // components containing onlhy the TE-mode
                free(GammaB);
                free(dGammaB);
                free(Ptot);
                free(dPtot);
            }else if (component[0] == 11 || component[0] == 12 || component[0] == 21 || component[0] == 22 || component[0] == 41 || component[0] == 42 || component[0] == 51 || component[0] == 52) { // components containing TE- as well as TM-modes
                free(Gamma);
                free(dGamma);
                free(GammaB);
                free(dGammaB);
                free(Ptot);
                free(dPtot);
            }
        } else { // if the model is a homogeneous fullspace
            free(Ptot);
        }
    }
    free(etaV);
    free(etaH);
    free(zetaV);
    free(zetaH);
    free(gamA);
    free(gamB);
    free(cor);
    free(tempcor);
    
    return;
}
