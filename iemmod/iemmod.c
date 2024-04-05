#include "emmod.h"
#include "float.h"
#include <errno.h>

/************************ self documentation ***********************/
char *sdoc[] = {
    "                                                                 ",
    " iemmod - non-linear conjugate-gradient inversion algorithm      ",
    "          to invert for electric conductivity in a layered earth.",
    "                                                                 ",
    " Required parameters:                                            ",
    "                                                                 ",
    " nx ................. number of x-samples                        ",
    " ny ................. number of y-samples                        ",
    " dx ................. x-sampling                                 ",
    " dy ................. y-sampling                                 ",
    " zsrc ............... depth of source                            ",
    " zrcv ............... depth of receiver                          ",
    " file_weight ........ location of weighting file                 ",
    " file_in ............ location of input data file                ",
    " path_out ........... folder for storing output                  ",
    " z .................. depth of layer interfaces                  ",
    " econd .............. starting electric conductivity             ",
    " eperm .............. starting rel. electric permittivity        ",
    " mperm .............. starting rel. magnetic permeability        ",
    "                                                                 ",
    " Optional parameters (default in brackets):                      ",
    "                                                                 ",
    " freqvec ............ frequency of the input data (0.5)          ",
    " known .............. 1: conductivity of that layer is known     ",
    "                      0: conductivity of that layer is unknown   ",
    "                      (0)                                        ",
    " component .......... receiver and source orientation (11)       ",
    " nd ................. number of integration domains (1000)       ",
    " startlogx .......... first integration point in space (-6)      ",
    " deltalogx .......... logarithmic sampling rate of integration   ",
    "                      points in space at first iteration (0.025) ",
    " nlogx .............. amount of integration points in space      ",
    "                      at first iteration (512)                   ",
    " kmax ............... largest wavenumber to be integrated        ",
    "                      (0.628625)                                 ",
    " c1 ................. first precision parameter (0.0)            ",
    " c2 ................. second precision parameter (0.01)          ",
    " maxpt .............. maximum amount of integration points (500) ",
    " xdirect ............ 1: direct field in space domain            ",
    "                      0: direct field in the wavenumber domain   ",
    "                      (1)                                        ",
    " whichstep .......... 1: step length determined by parabola      ",
    "                         fitting                                 ",
    "                      2: step length determined by grid search   ",
    "                      (2)                                        ",
    " whichdir ........... direction of next update using             ",
    "                      1: the Fletcher-Reeves formula             ",
    "                      2: the Polak-Ribière formula               ",
    "                      (1)                                        ",
    " eps ................ critical model update size (1e-5)          ",
    " maxit .............. maximum amount of iterations (500)         ",
    " restartgrad ........ amount of iterations after which the       ",
    "                      direction of the update is set to the      ",
    "                      negative gradient for one iteration (5)    ",
    " ngradrestart ....... amount of extraordinary gradient restarts  ",
    "                      (4)                                        ",
    " epsa ............... threshold to stop the line search for      ",
    "                      determining the step length (1e-5)         ",
    " maxita ............. maximum amount of iterations for           ",
    "                      determining the step length (20)           ",
    " maxnogood .......... maximum amount of line search trials       ",
    "                      without improvement (5)                    ",
    NULL};
/***************** end self doc ************************************/

/*
 * By using this software, you agree to read and agree to the
 * disclaimer at:
 * http://software.seg.org/disclaimer2.txt
 *
 * The SEG static snapshot of this code may be found at:
 * http://software.seg.org/2016/0005
 */

REAL compmisfit(complex *datatrue, complex *data, complex *weight, int nel);
void compgrad(REAL *gradient, complex *datatrue, complex *data, complex *ddata, complex *weight, int nrz, int nel);
void inv3x3(REAL *invA, REAL *A);
REAL scalprod(REAL *v1, REAL *v2, int nel);
REAL vecsum(REAL *vec, int nel);

int main(int argc, char **argv)
{
    FILE    *fp_in;
    int     nrz, xdirect, maxpt, nread, in, im, il, ia, m;
    int     *component, *onecomponent, *known, nx, ny, verbose, nd, nlogx, donumgrad, nfreq, ncomp, nweight, nfiles, ifile;
    int     maxit, maxita, whichstep, whichdir, lastdir, switcher;
    int     restartgrad, ngradrestart, maxnogood, newgrad, stopcounter, inc, dostop, firsttry, nogood;
    REAL    kmax, dx, dy, freq, *freqvec, *z, *d, *r, *rold, *tempvec, eps, epsa, alpha, alphaprev, dalpha, beta;
    REAL    realtemp, *kernel, *invkernel, *coef, userfac, deltacond;
    REAL    *econd, *econd1, *econd2, *econdold, *econdtemp, *eperm, *mperm, *grad, *gradold, *gradtemp;
    REAL    zrcv, zsrc, startlogx, deltalogx, c1, c2, lengthstep, normfac, invnormfac;
    REAL    cosgam, gam, misfit, misfittempprev, misfit1, misfit2, misfitold, misfitup, misfitdo;
    complex *datatrue, *data, *datatemp, *data1, *data2, *ddata, *ddatatemp, *weight;
    char    **file_in, **file_weight, *path_out, *filenameout;
    double  t0, t1;
    
    initargs(argc, argv);
    requestdoc(1);
    
    /* Time when program was started. */
    t0=wallclock_time();
    
    vmess("Welcome to iEMmod");
    vmess("Inverting electromagnetic fields with a conjugate gradient scheme using EMmod as forward solver.");
    
    /* Mandatory Parameters */
    if(!getparfloat("zsrc", &zsrc)) verr("zsrc (source depth) must be specified.");
    if(!getparfloat("zrcv", &zrcv)) verr("zrcv (receiver depth) must be specified.");
    if(!getparfloat("dx", &dx)) verr("dx (x-sampling) must be specified.");
    if(!getparfloat("dy", &dy)) verr("dy (y-sampling) must be specified.");
    if(!getparint("nx", &nx)) verr("nx (number of x-samples) must be specified.");
    if(!getparint("ny", &ny)) verr("ny (number of y-samples) must be specified.");
    /* If an uneven number is entered for nx or ny, it is made even. */
    if (fabs(nx/2.0-ceil(nx/2.0))>=0.000001) {nx = nx+1;}
    if (fabs(ny/2.0-ceil(ny/2.0))>=0.000001) {ny = ny+1;}
    
    /* Optional Parameters */
    nfiles = countparval("file_in");
    nweight = countparval("file_weight");
    nfreq  = countparval("freqvec");
    ncomp  = countparval("component");
    if (nfiles != nweight) verr("The amount of weighting files must be equal to the amount of data files specified.");
    if (nfiles != nfreq) verr("The length of the frequency vector must be equal to the amount of data files specified.");
    if (nfiles != ncomp) verr("The length of the component vector must be equal to the amount of data files specified.");
    freqvec = (REAL *)calloc(nfiles,sizeof(REAL));
    if(!getparfloat("freqvec", freqvec)) freqvec[0] = 0.5;
    //if(!getparfloat("freq", &freq)) freq = 0.5;
    if(!getparint("nd", &nd)) nd = 1000;
    if(!getparfloat("kmax", &kmax)) kmax = 0.628625;
    if(!getparfloat("startlogx", &startlogx)) startlogx = -6;
    if(!getparfloat("deltalogx", &deltalogx)) deltalogx = 0.025;
    if(!getparint("nlogx", &nlogx)) nlogx = 512;
    if(!getparfloat("c1", &c1)) c1 = 0.0;
    if(!getparfloat("c2", &c2)) c2 = 0.01;
    if(!getparint("maxpt", &maxpt)) maxpt = 500;
    if(!getparint("xdirect", &xdirect)) xdirect = 1;
    
    /* Inversion Parameters */
    // Determine the step length by (1) fitting of a parabola or (2) by a line search.
    if(!getparint("whichstep", &whichstep)) whichstep = 2;
    // Use (1) the Fletcher-Reeves formula to determine the search direction or (2) the Polak-Ribière formula.
    if(!getparint("whichdir", &whichdir)) whichdir = 1;
    // Tolerance to stop search, i.e., the distance a new solution is away from the current solution |alpha*d|
    if(!getparfloat("eps", &eps)) eps = 1e-5;
    // Maximum amount of iterations
    if(!getparint("maxit", &maxit)) maxit = 500;
    // Amount of iterations after which the gradient is restarted.
    if(!getparint("restartgrad", &restartgrad)) restartgrad = 5;
    // How many times an extraordinary restart of the gradient is allowed before the program stops.
    if(!getparint("ngradrestart", &ngradrestart)) ngradrestart = 4;
    // Tolerance to stop search for the step length alpha
    if(!getparfloat("epsa", &epsa)) epsa = 1e-5;
    // Maximum amount of iterations to find the step length alpha
    if(!getparint("maxita", &maxita)) maxita = 20;
    // If the misfit during the line search for the step length alpha has not improved
    // for the following amount of steps in a row, the line search is stopped and
    // the previously best alpha is selected.
    if(!getparint("maxnogood", &maxnogood)) maxnogood = 5;
    // Source and receiver orientation
    component = (int *)calloc(nfiles,sizeof(int));
    onecomponent = (int *)calloc(1,sizeof(int));
    if(!getparint("component", component)) component[0]=11;
    
    /* Read path to write results. */
    if(!getparstring("path_out", &path_out)) verr("path to write results must be specified.");
    vmess("Write files to: %s\n",path_out);
    
    // Read name of datafiles and weighting files
    file_in=(char **)malloc(nfiles*sizeof(char*));
    file_weight=(char **)malloc(nweight*sizeof(char*));
    for (ifile=0; ifile<nfiles; ifile++) {
        file_in[ifile]=(char *)malloc(512*sizeof(char));
        file_weight[ifile]=(char *)malloc(512*sizeof(char));
    }
    getparstringarray("file_in", file_in);
    getparstringarray("file_weight", file_weight);
    for (ifile=0; ifile<nfiles; ifile++) {
        vmess("Data file %d: %s\n",ifile,file_in[ifile]);
        vmess("Weighting file %d: %s\n",ifile,file_weight[ifile]);
    }
    
    /* Allocate matrix for measured data. */
    datatrue = (complex *)calloc(nx*ny*nfiles,sizeof(complex));
    
    /* Read measured data. */
    for (ifile=0; ifile<nfiles; ifile++) {
        fp_in = fopen(file_in[ifile],"r");
        if (fp_in == NULL) verr("error in reading measured data (file_in)");
        nread = fread( &datatrue[ifile*(nx*ny)].r, sizeof(complex), nx*ny, fp_in);
        assert(nread == nx*ny);
        fclose(fp_in);
    }
    
    /* Allocate matrix for weights. */
    // Also weight is complex, only the real part is used.
    // If a complex weight is used, the functions to compute the misfit and the gradient need to be adjusted.
    weight = (complex *)calloc(nx*ny*nfiles,sizeof(complex));
    
    /* Read weights data. */
    vmess("nweight=%d",nweight);
    for (ifile=0; ifile<nweight; ifile++) {
        fp_in = fopen(file_weight[ifile],"r");
        if (fp_in == NULL) verr("error in reading weights (file_weight)");
        nread = fread( &weight[ifile*(nx*ny)].r, sizeof(complex), nx*ny, fp_in);
        assert(nread == nx*ny);
        fclose(fp_in);
    }
    
    /* Read medium parameters of starting model. */
    nrz  = countparval("z");
    z  = (REAL *)calloc(nrz,sizeof(REAL));
    econd = (REAL *)calloc(nrz,sizeof(REAL));
    econd1 = (REAL *)calloc(nrz,sizeof(REAL));
    econd2 = (REAL *)calloc(nrz,sizeof(REAL));
    econdold = (REAL *)calloc(nrz,sizeof(REAL));
    econdtemp = (REAL *)calloc(nrz,sizeof(REAL));
    eperm = (REAL *)calloc(nrz,sizeof(REAL));
    mperm = (REAL *)calloc(nrz,sizeof(REAL));
    known = (int *)calloc(nrz,sizeof(int));
    getparfloat("z", z);
    getparfloat("econd", econd);
    getparfloat("eperm", eperm);
    getparfloat("mperm", mperm);
    if(!getparint("known", known)) {
        for (im=0; im<nrz; im++) {
            known[im]=0;
        }
    }
    
    /* In the forward simulation no information is printed to the screen. */
    verbose = 0;
    
    /* Allocate matrices for forward simulated data, gradient and CG algorithm. */
    data = (complex *)calloc(nx*ny,sizeof(complex));
    data1 = (complex *)calloc(nx*ny,sizeof(complex));
    data2 = (complex *)calloc(nx*ny,sizeof(complex));
    ddata = (complex *)calloc(nx*ny*nrz,sizeof(complex));
    ddatatemp = (complex *)calloc(nx*ny*nrz,sizeof(complex));
    grad = (REAL *)calloc(nrz,sizeof(REAL));
    gradold = (REAL *)calloc(nrz,sizeof(REAL));
    gradtemp = (REAL *)calloc(nrz,sizeof(REAL));
    d = (REAL *)calloc(nrz,sizeof(REAL));
    r = (REAL *)calloc(nrz,sizeof(REAL));
    rold = (REAL *)calloc(nrz,sizeof(REAL));
    tempvec = (REAL *)calloc(nrz,sizeof(REAL));
    kernel = (REAL *)calloc(9,sizeof(REAL));
    invkernel = (REAL *)calloc(9,sizeof(REAL));
    coef = (REAL *)calloc(3,sizeof(REAL));
    
    // User determined scaling factor for alpha used in the parabola fitting
    userfac = 1.0/1000.0;
    
    // Compute the gradient numerically (1) or analytically (0)
    donumgrad = 0;
    
    // The value added and subtracted from the conductivity to compute the gradient numerically
    deltacond = 5e-4;
    
    // Print information to the screen
    vmess("Parameter starting values:");
    for (im=0; im<nrz; im++) {
        vmess("econd[%d] = %f, known=%d",im,econd[im],known[im]);
        //vmess("econd[%d] = %f, eperm[%d] = %f, mperm[%d] = %f, known=%d",im,econd[im],im,eperm[im],im,mperm[im],known[im]);
    }
    
    /* Compute forward model and gradient for starting model */
    if (donumgrad == 0) {// compute the gradient analytically
        vmess("Computing forward model and gradient for starting model...");
        // Initialize misfit and gradient
        misfit = 0.0;
        for (im=0; im<nrz; im++) {
            grad[im] = 0.0;
        }
        for (ifile=0; ifile<nfiles; ifile++) {// Loop over amount of input data files
            onecomponent[0] = component[ifile];
            freq = freqvec[ifile];
            vmess("The component %d is %d",ifile,onecomponent[0]);
            vmess("The frequency %d is %f",ifile,freq);
            // call the forward routine to compute semi-analytically the forward problem and the derivatives of it
            emmod(data, ddata, freq, nx, dx, ny, dy, nrz, z, econd, econd, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 1, verbose);
            // compute the misfit
            realtemp = compmisfit(&datatrue[(nx*ny)*ifile], data, &weight[(nx*ny)*ifile], nx*ny);
            // the misfit of all data files are summed up
            misfit = misfit + realtemp;
            // compute the gradient of the misfit function
            compgrad(gradtemp, &datatrue[(nx*ny)*ifile], data, ddata, &weight[(nx*ny)*ifile], nrz, nx*ny);
            for (im=0; im<nrz; im++) {
                grad[im] = grad[im] + gradtemp[im];
            }
        }
    } else {// compute the gradient numerically
        vmess("Computing forward model for starting model...");
        misfit = 0.0;
        for (ifile=0; ifile<nfiles; ifile++) {
            onecomponent[0] = component[ifile];
            freq = freqvec[ifile];
            // call the forward routine to compute semi-analytically the forward problem
            emmod(data, ddata, freq, nx, dx, ny, dy, nrz, z, econd, econd, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
            // compute the misfit
            realtemp = compmisfit(&datatrue[(nx*ny)*ifile], data, &weight[(nx*ny)*ifile], nx*ny);
            // the misfit of all data files are summed up
            misfit = misfit + realtemp;
        }
        vmess("Computing gradient for starting model numerically...");
        // to compute the gradient numerically, for each layer the misfit is computed for a slightly too high and a slightly too low conductivity
        for (im=0; im<nrz; im++) {// Loop over layers
            // copy the conductivity model into another array
            for (il=0; il<nrz; il++) {
                econdtemp[il] = econd[il];
            }
            vmess("Layer %d of %d",im+1,nrz);
            vmess("Upper Misfit");
            // the conductivity of layer im is slightly increased
            econdtemp[im] = econd[im]+deltacond;
            datatemp = (complex *)calloc(nx*ny,sizeof(complex));
            misfitup = 0.0;
            for (ifile=0; ifile<nfiles; ifile++) {// Loop over amount of input data files
                onecomponent[0] = component[ifile];
                freq = freqvec[ifile];
                // Forward problem computing the misfit with the conductivity of layer im slightly too high
                emmod(datatemp, ddata, freq, nx, dx, ny, dy, nrz, z, econdtemp, econdtemp, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                realtemp = compmisfit(&datatrue[(nx*ny)*ifile], datatemp, &weight[(nx*ny)*ifile], nx*ny);
                misfitup = misfitup + realtemp;
            }
            free(datatemp);
            vmess("Lower Misfit");
            // the conductivity of layer im is slightly decrased
            econdtemp[im] = econd[im]-deltacond;
            datatemp = (complex *)calloc(nx*ny,sizeof(complex));
            misfitdo = 0.0;
            for (ifile=0; ifile<nfiles; ifile++) {
                onecomponent[0] = component[ifile];
                freq = freqvec[ifile];
                // Forward problem computing the misfit with the conductivity of layer im slightly too low
                emmod(datatemp, ddata, freq, nx, dx, ny, dy, nrz, z, econdtemp, econdtemp, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                realtemp = compmisfit(&datatrue[(nx*ny)*ifile], datatemp, &weight[(nx*ny)*ifile], nx*ny);
                misfitdo = misfitdo + realtemp;
            }
            free(datatemp);
            // the gradient of the misfit function for layer im is computed using the two values of the misfit function previously determined
            grad[im] = (misfitup-misfitdo)/(2.0*deltacond);
        }
    }
    vmess("misfit = %e",misfit);
    for (im=0; im<nrz; im++) {
        vmess("gradient[%d] = %e",im,grad[im]);
    }
    
    vmess("-----------------------------------------------------------------------\n");
    
    // The nonlinear conjugate-gradient scheme starts here. The implementation follows:
    // Shewchuk, J. R., 1994, An Introduction to the Conjugate Gradient Method Without
    // the Agonizing Pain, School of Computer Science Carnegie Mellon University Pittsburgh
    
    // The negative gradient of the starting model are the initial values of the direction d
    // and the residual r
    for (im=0; im<nrz; im++) {
        d[im] = -grad[im];
        r[im] = -grad[im];
    }
    
    switcher = 0;
    newgrad = 0;
    in = 1;
    lengthstep = 1.0;
    stopcounter = 0;
    // loop until convergence is achieved or the maximum amount of iterations reached
    while (stopcounter<ngradrestart && in<=maxit) {
        if (whichstep == 1) {// Determining the steplength by fitting a parabola
            // The misfit function in the direction d is approximated by a parabola.
            // If three points of the misfit function along d are known, a parabola
            // of the form y = a*x^2 + b*x + c can be fitted with a simple linear
            // inverse problem with the unknowns a, b, c (in this case the elements
            // of the vector coef).
            // The best step-length will bring the model to the minimum of the parabola.
            vmess("Entering Parabola fitting.");
            normfac = vecsum(d,nrz);
            invnormfac = 1.0/normfac;
            // The kernel is the matrix A of the linear inverse problem A*x = b
            kernel[0] = 0.0;
            kernel[1] = (1.0*userfac*invnormfac)*(1.0*userfac*invnormfac);
            kernel[2] = (2.0*userfac*invnormfac)*(2.0*userfac*invnormfac);
            kernel[3] = 0.0;
            kernel[4] = 1.0*userfac*invnormfac;
            kernel[5] = 2.0*userfac*invnormfac;
            kernel[6] = 1.0;
            kernel[7] = 1.0;
            kernel[8] = 1.0;
            // Determine two more models that are used to fit the parabola
            // Note, that this choice is connected how the kernel-variable is set up.
            for (im=0; im<nrz; im++) {
                if (known[im]==0){ // 0 means that we want to estimate that conductivity value.
                    econd1[im] = econd[im]+1.0*userfac*invnormfac*d[im];
                    econd2[im] = econd[im]+2.0*userfac*invnormfac*d[im];
                    if (econd1[im]<0) {
                        econd1[im] = 0.0;
                    }
                    if (econd2[im]<0) {
                        econd2[im] = 0.0;
                    }
                } else {
                    econd1[im] = econd[im];
                    econd2[im] = econd[im];
                }
            }
            // Three points are needed to fit a parabola to them.
            // The first point is the misfit of the previous model already computed. Thus, the misfit needs
            // to be computed only for two more models.
            // In the following two calls to emmod, only the data, but not the gradient are required
            vmess("Determine second parabola point...");
            misfit1 = 0.0;
            for (ifile=0; ifile<nfiles; ifile++) {// Loop over amount of input data files
                onecomponent[0] = component[ifile];
                freq = freqvec[ifile];
                // Forward problem computing the misfit using model econd1
                emmod(data1, ddatatemp, freq, nx, dx, ny, dy, nrz, z, econd1, econd1, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                realtemp = compmisfit(&datatrue[(nx*ny)*ifile], data1, &weight[(nx*ny)*ifile], nx*ny);
                misfit1 = misfit1 + realtemp;
            }
            vmess("Determine third parabola point...");
            misfit2 = 0.0;
            for (ifile=0; ifile<nfiles; ifile++) {// Loop over amount of input data files
                onecomponent[0] = component[ifile];
                freq = freqvec[ifile];
                // Forward problem computing the misfit using model econd2
                emmod(data2, ddatatemp, freq, nx, dx, ny, dy, nrz, z, econd2, econd2, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                realtemp = compmisfit(&datatrue[(nx*ny)*ifile], data2, &weight[(nx*ny)*ifile], nx*ny);
                misfit2 = misfit2 + realtemp;
            }
            // Compute the inverse of the kernel
            inv3x3(invkernel,kernel);
            // Determine the three parameters describing the fitted parabola by a matrix multiplication
            coef[0] = misfit*invkernel[0]+misfit1*invkernel[3]+misfit2*invkernel[6];
            coef[1] = misfit*invkernel[1]+misfit1*invkernel[4]+misfit2*invkernel[7];
            coef[2] = misfit*invkernel[2]+misfit1*invkernel[5]+misfit2*invkernel[8];
            // alpha is the step-length that will move the model to the minimum of the parabola
            alpha = -coef[1]/(2.0*coef[0]);
            if (alpha<0.0) { // A negative alpha means that we are walking 180 degrees in the wrong direction
                alpha = invnormfac;
            }
            
        } else if (whichstep == 2) {// Determining the steplegnth by a line search
            // The line search computes various misfits along the direction d to check, where there is a minimum.
            vmess("Entering line search.");
            misfittempprev = misfit;
            normfac = vecsum(d,nrz);
            invnormfac = 1.0/normfac;
            dalpha = 100.0;
            alpha = dalpha;
            alphaprev = alpha;
            inc = 0;
            dostop = 0;
            firsttry = 1;
            nogood = 0;
            ia = 1;
            while (dostop==0) {// stops when the line search has found a step-length
                // compute the model for a first guess of alpha, i.e., the step-length
                for (im=0; im<nrz; im++) {
                    if (known[im]==0){ // 0 means that we want to estimate that conductivity value.
                        econd1[im] = econd[im]+alpha*invnormfac*d[im];
                        if (econd1[im]<0) {
                            econd1[im] = 0.0;
                        }
                    } else {
                        econd1[im] = econd[im];
                    }
                }
                // In the following call to emmod, only the data, but not the gradient is required
                misfit1 = 0.0;
                for (ifile=0; ifile<nfiles; ifile++) {// Loop over amount of input data files
                    onecomponent[0] = component[ifile];
                    freq = freqvec[ifile];
                    // Forward problem computing the misfit using model econd1
                    emmod(data1, ddatatemp, freq, nx, dx, ny, dy, nrz, z, econd1, econd1, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                    realtemp = compmisfit(&datatrue[(nx*ny)*ifile], data1, &weight[(nx*ny)*ifile], nx*ny);
                    misfit1 = misfit1 + realtemp;
                }
                vmess("Iteration %d: alpha = %e; J = %e; J0 = %e",ia,alpha,misfit1,misfit);
                if (misfit1<misfittempprev) {// The misfit of the new model is smaller than the previous one
                    misfittempprev = misfit1;
                    alphaprev = alpha; // The step-length is accepted
                    while (dalpha>alpha) {
                        dalpha = dalpha/2.0;
                    }
                    alpha = alpha+dalpha;// A new step-length is computed
                    firsttry=0;
                    nogood=0;
                } else {// The misfit of the new model is larger than the previous one
                    if (firsttry==1) {// No alpha so far has been an improvement
                        alphaprev=alphaprev/10.0; // The step-lenght might have been too large and is therefore reduced
                        alpha=alphaprev;
                    } else {// The best estimate of alpha so far is used to determine the next guess of alpha
                        alpha = alphaprev;
                        dalpha = dalpha/2.0;
                        alpha = alpha+dalpha;
                        inc = 1;
                        nogood = nogood+1;
                    }
                }
                // The various rules that stop the grid search:
                // 1.) The maximum amount of iterations is reached and a good alpha has been found.
                // 2.) The maximum amount of iterations is reached, but no good alpha has been found.
                // 3.) The change of alpha (dalpha) has become smaller than the critical value
                // 4.) To many times in a row no improvement has been achieved
                if ((inc==1 && ia>maxita) || (firsttry==1 && ia>maxita) || dalpha<epsa || nogood>maxnogood) {
                    dostop = 1;
                    alpha = alphaprev;
                }
                ia = ia + 1; // counts the amount of iterations of the grid-search
            }
            vmess("Final alpha = %e\n",alpha);
            alpha = alpha*invnormfac;
            if (switcher==1) {
                whichstep = 1;
                switcher = 0;
            }
            
        }
        
        // The model is updated according to the step-length found with either method above
        for (im=0; im<nrz; im++) {
            econdold[im] = econd[im];
            if (known[im]==0){ // 0 means that we want to estimate that conductivity value.
                tempvec[im] = alpha*d[im];
            } else {
                tempvec[im] = 0.0;
            }
            econd[im] = econd[im]+tempvec[im];
            if (econd[im]<0) {
                econd[im] = 0.0;
            }
        }
        lengthstep = vecsum(tempvec,nrz);
        misfitold = misfit;
        for (im=0; im<nrz; im++) {
            gradold[im] = grad[im];
        }
        
        // What was already done for the starting model is repeated here:
        // The misfit and the gradient are determined. The latter either analytically or numerically.
        if (donumgrad == 0) {
            vmess("Computing forward model and gradient for current model...");
            misfit = 0.0;
            for (im=0; im<nrz; im++) {
                grad[im] = 0.0;
            }
            for (ifile=0; ifile<nfiles; ifile++) {
                onecomponent[0] = component[ifile];
                freq = freqvec[ifile];
                // In the following evaluation of the forward problem, both the field and the gradient are required
                emmod(data, ddata, freq, nx, dx, ny, dy, nrz, z, econd, econd, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 1, verbose);
                realtemp = compmisfit(&datatrue[(nx*ny)*ifile], data, &weight[(nx*ny)*ifile], nx*ny);
                misfit = misfit + realtemp;
                compgrad(gradtemp, &datatrue[(nx*ny)*ifile], data, ddata, &weight[(nx*ny)*ifile], nrz, nx*ny);
                for (im=0; im<nrz; im++) {
                    grad[im] = grad[im] + gradtemp[im];
                }
            }
        } else {
            vmess("Computing forward model for current model...");
            misfit = 0.0;
            for (ifile=0; ifile<nfiles; ifile++) {
                onecomponent[0] = component[ifile];
                freq = freqvec[ifile];
                emmod(data, ddata, freq, nx, dx, ny, dy, nrz, z, econd, econd, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                realtemp = compmisfit(&datatrue[(nx*ny)*ifile], data, &weight[(nx*ny)*ifile], nx*ny);
                misfit = misfit + realtemp;
            }
            vmess("Computing gradient for current model numerically...");
            for (im=0; im<nrz; im++) {
                for (il=0; il<nrz; il++) {
                    econdtemp[il] = econd[il];
                }
                vmess("Layer %d of %d",im+1,nrz);
                vmess("Upper Misfit");
                econdtemp[im] = econd[im]+deltacond;
                datatemp = (complex *)calloc(nx*ny,sizeof(complex));
                misfitup = 0.0;
                for (ifile=0; ifile<nfiles; ifile++) {
                    onecomponent[0] = component[ifile];
                    freq = freqvec[ifile];
                    emmod(datatemp, ddata, freq, nx, dx, ny, dy, nrz, z, econdtemp, econdtemp, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                    realtemp = compmisfit(&datatrue[(nx*ny)*ifile], datatemp, &weight[(nx*ny)*ifile], nx*ny);
                    misfitup = misfitup + realtemp;
                }
                free(datatemp);
                vmess("Lower Misfit");
                econdtemp[im] = econd[im]-deltacond;
                datatemp = (complex *)calloc(nx*ny,sizeof(complex));
                misfitdo = 0.0;
                for (ifile=0; ifile<nfiles; ifile++) {
                    onecomponent[0] = component[ifile];
                    freq = freqvec[ifile];
                    emmod(datatemp, ddata, freq, nx, dx, ny, dy, nrz, z, econdtemp, econdtemp, eperm, eperm, mperm, mperm, zsrc, zrcv, onecomponent, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, xdirect, 0, verbose);
                    realtemp = compmisfit(&datatrue[(nx*ny)*ifile], datatemp, &weight[(nx*ny)*ifile], nx*ny);
                    misfitdo = misfitdo + realtemp;
                }
                free(datatemp);
                grad[im] = (misfitup-misfitdo)/(2.0*deltacond);
            }
        }
        
        // Print information to the screen
        vmess("Parameter update after iteration %d:",in);
        for (im=0; im<nrz; im++) {
            vmess("econd[%d] = %f",im,econd[im]);
        }
        vmess("misfit = %e",misfit);
        for (im=0; im<nrz; im++) {
            vmess("gradient[%d] = %e",im,grad[im]);
        }
        if (misfit>misfitold) {// The misfit for this iteration has increased compared to the previous iteration
            vmess("Function value for next solution has increased.\n");
            // The previous model is reloaded
            for (im=0; im<nrz; im++) {
                econd[im] = econdold[im];
                grad[im] = gradold[im];
            }
            misfit = misfitold;
            if (whichstep == 1) {
                // If parabola fitting has been used to find the step-length, next time we do a line-search to find the step-length.
                // This is switched back after the next iteration to parabola fitting.
                vmess("Find a new alpha with a line search.\n");
                whichstep = 2;
                switcher = 1;
            }
        }
        // keeping the old residual and computing the new one
        for (im=0; im<nrz; im++) {
            rold[im] = r[im];
            r[im] = -1.0*grad[im];
        }
        // The step-length is smaller than the critical value (lengthstep<eps), but actually the algorithm could not find a suitable step-length (firsttry==1).
        // The direction is exceptionally restarted using the negative gradient (newgrad = 1).
        // To ensure, that the algorithm is not thinking that convergence has been reached, lengthstep is changed to 10 times the critical value.
        if (lengthstep<eps && firsttry==1) {
            lengthstep=10*eps;
            newgrad = 1;
        }
        // The step-length is smaller than the critical value indicates that convergence might be reached. Therefore the stopcounter is increased.
        // Still, the direction is exceptionally restarted using the negative gradient (newgrad = 1).
        if (lengthstep<eps) {
            stopcounter = stopcounter+1;
            newgrad = 1;
        }
        // The direction is restarted using the negative gradient. This is done regularly with a period of restartgrad iterations or exceptionally if newgrad=1.
        if (fabs(fmod(in,restartgrad))<=1e-6 || newgrad==1) {
            for (im=0; im<nrz; im++) {
                d[im] = -1.0*grad[im];
            }
        } else {// If the direction is not restarted one of the two following formulas is used to compute the direction of the next iteration.
            if (whichdir == 1) {// Fletcher-Reeves
                beta = scalprod(r,r,nrz)/scalprod(rold,rold,nrz);
            } else {// Polak-Ribière
                for (im=0; im<nrz; im++) {
                    tempvec[im] = r[im]-rold[im];
                }
                beta = scalprod(r,tempvec,nrz)/scalprod(rold,rold,nrz);
            }
            for (im=0; im<nrz; im++) {
                d[im] = r[im] + beta*d[im];
            }
        }
        vmess("-----------------------------------------------------------------------\n");
        
        // The following commented lines can be used to write intermediate data to a binary file for debugging purposes.
        /*{
         int nwrite;
         FILE *fp_out = fopen("/data_a/hunziker/testing/deriv/intermediate.bin","w");
         nwrite = fwrite( &data[0].r, sizeof(complex), nx*ny, fp_out);
         assert(nwrite == nx*ny);
         fclose(fp_out);
         }
         
         {
         int nwrite;
         FILE *fp_out = fopen("/data_a/hunziker/testing/deriv/intermediate_deriv.bin","w");
         nwrite = fwrite( &ddata[0].r, sizeof(complex), nx*ny*nrz, fp_out);
         assert(nwrite == nx*ny*nrz);
         fclose(fp_out);
         }*/
        
        
        /* Writing data to harddisk */
        filenameout = (char *)calloc(200,sizeof(char));
        sprintf(filenameout,"%s/econd_%d.bin\0",path_out,in);
        {
            int nwrite;
            FILE *fp_out = fopen(filenameout,"w");
            nwrite = fwrite( &econd[0], sizeof(REAL), nrz, fp_out);
            assert(nwrite == nrz);
            nwrite = fwrite( &misfit, sizeof(REAL), 1, fp_out);
            assert(nwrite == 1);
            fclose(fp_out);
        }
        free(filenameout);
        in = in + 1;
    } // end while
    
    /* Time when program is completed. */
    t1=wallclock_time();
    
    /* Free memory */
    free(component);
    free(onecomponent);
    free(datatrue);
    free(weight);
    free(z);
    free(econd);
    free(econd1);
    free(econd2);
    free(econdold);
    free(econdtemp);
    free(eperm);
    free(mperm);
    free(known);
    free(data);
    free(data1);
    free(data2);
    free(ddata);
    free(ddatatemp);
    free(grad);
    free(gradold);
    free(gradtemp);
    free(d);
    free(r);
    free(rold);
    free(tempvec);
    free(kernel);
    free(invkernel);
    free(coef);
    
    /* Show runtimes. */
    vmess("iemmod runtime: ");
    vmess(" - time=%.3f", t1-t0);
    
    exit ( 0 );
}

REAL compmisfit(complex *datatrue, complex *data, complex *weight, int nel)
{// Compute the misfit function
    
    int     im;
    REAL    misfit, realtemp;
    complex    tempmisfit;
    
    misfit = 0.0;
    for (im=0; im<nel; im++) {
        tempmisfit.r = (datatrue[im].r - data[im].r)*(datatrue[im].r - data[im].r);
        tempmisfit.i = (datatrue[im].i - data[im].i)*(datatrue[im].i - data[im].i);
        realtemp = weight[im].r*weight[im].r*(tempmisfit.r + tempmisfit.i);
        misfit = misfit + realtemp;
    }
    //vmess("misfit = %e",misfit);
    return misfit;
}

void compgrad(REAL *gradient, complex *datatrue, complex *data, complex *ddata, complex *weight, int nrz, int nel)
{// compute the gradient
    
    int     in, im;
    REAL    gradfac;
    complex    tempmisfit, signfun;
    
    for (in=0; in<(nrz); in++) {
        gradient[in] = 0.0;
        for (im=0; im<nel; im++) {
            gradfac = -2.0*weight[im].r*weight[im].r;
            tempmisfit.r = (datatrue[im].r - data[im].r)*ddata[im+in*nel].r;
            tempmisfit.i = (datatrue[im].i - data[im].i)*ddata[im+in*nel].i;
            gradient[in] = gradient[in] + gradfac*(tempmisfit.r+tempmisfit.i);
        }
        //vmess("gradient[%d] = %e",in,gradient[in]);
    }
}

void inv3x3(REAL *invA, REAL *A)
{// Inverts a 3 by 3 matrix
    REAL    a, b, c, d, e, f, g, h, i;
    REAL    detA, invdetA;
    
    a =  (A[4]*A[8]-A[5]*A[7]);
    b = -(A[3]*A[8]-A[5]*A[6]);
    c =  (A[3]*A[7]-A[4]*A[6]);
    d = -(A[1]*A[8]-A[2]*A[7]);
    e =  (A[0]*A[8]-A[2]*A[6]);
    f = -(A[0]*A[7]-A[1]*A[6]);
    g =  (A[1]*A[5]-A[2]*A[4]);
    h = -(A[0]*A[5]-A[2]*A[3]);
    i =  (A[0]*A[4]-A[1]*A[3]);
    
    detA = A[0]*(A[4]*A[8]-A[5]*A[7]);
    detA = detA -A[1]*(A[8]*A[3]-A[5]*A[6]);
    detA = detA -A[2]*(A[3]*A[7]-A[4]*A[6]);
    invdetA = 1.0/detA;
    
    invA[0] = invdetA*a;
    invA[1] = invdetA*d;
    invA[2] = invdetA*g;
    invA[3] = invdetA*b;
    invA[4] = invdetA*e;
    invA[5] = invdetA*h;
    invA[6] = invdetA*c;
    invA[7] = invdetA*f;
    invA[8] = invdetA*i;
    
}

REAL scalprod(REAL *v1, REAL *v2, int nel)
{// Computes the scalar product of two vectors
    
    int     im;
    REAL    floattemp;
    
    floattemp = 0.0;
    for (im=0; im<nel; im++) {
        floattemp = floattemp + v1[im]*v2[im];
    }
    return floattemp;
}

REAL vecsum(REAL *vec, int nel)
{// Computes the sum of a vector
    
    int     im;
    REAL     floattemp;
    
    floattemp = 0.0;
    for (im=0; im<nel; im++) {
        floattemp = floattemp + vec[im]*vec[im];
    }
    return sqrt(floattemp);
}
