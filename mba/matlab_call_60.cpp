#include "mba.hpp"

#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t mrows0, mrows1, mrows2;  
    size_t ncols0, ncols1, ncols2;  
    double *in_bev, *in_voxel, *out_bev; 
    
    int i = 0;
    double w = 0, t1 = 0, t2 = 0, t3 = 0;
    double *ret;
    
    mrows0 = mxGetM(prhs[0]); 
    ncols0 = mxGetN(prhs[0]); 

    mrows1 = mxGetM(prhs[1]); 
    ncols1 = mxGetN(prhs[1]); 
    
    mrows2 = mxGetM(prhs[2]); 
    ncols2 = mxGetN(prhs[2]); 
    
    in_bev = mxGetPr(prhs[0]);
    in_voxel = mxGetPr(prhs[1]);
    out_bev = mxGetPr(prhs[2]);

    plhs[0]  = mxCreateDoubleMatrix(100,1,mxREAL);
    ret = mxGetPr(plhs[0]);
    
    std::vector<mba::point<3>> coo(60);
    for (i=0; i<60; i++)
    {
        t1 = in_bev[i];
        t2 = in_bev[i+60];
        t3 = in_bev[i+60*2];
        
        coo[i] = {t1, t2, t3};
    }
   

	// Data values.
	std::vector<double> val(60);
    for (i=0; i<60; i++)
    {
        t1 = in_voxel[i];
        
        val[i] = t1;
    }
    
	// Bounding box containing the data points.
	mba::point<3> lo = { -1, -1, -1};
	mba::point<3> hi = { 1, 1, 1};

	// Initial grid size.
	mba::index<3> grid = { 3, 3, 3 };

	// Algorithm setup.
	mba::MBA<3> interp(lo, hi, grid, coo, val);

	// Get interpolated value at arbitrary location.
	    
    for (i=0; i<100; i++)
    {
        t1 = out_bev[i];
        t2 = out_bev[i+100];
        t3 = out_bev[i+100*2];
        w = interp(mba::point<3>{t1, t2, t3});
        ret[i] = w;
    }
}

