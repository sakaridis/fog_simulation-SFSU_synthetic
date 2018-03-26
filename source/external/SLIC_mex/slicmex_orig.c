//=================================================================================
//  slicmex.c
//
//
//  AUTORIGHTS
//  Copyright (C) 2015 Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland.
//
//  Created by Radhakrishna Achanta on 12/01/15.
//=================================================================================
/*Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met
 
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of EPFL nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

void rgbtolab(int* rin, int* gin, int* bin, int sz, double* lvec, double* avec, double* bvec)
{
    int i; int sR, sG, sB;
    double R,G,B;
    double X,Y,Z;
    double r, g, b;
    const double epsilon = 0.008856;	//actual CIE standard
    const double kappa   = 903.3;		//actual CIE standard
    
    const double Xr = 0.950456;	//reference white
    const double Yr = 1.0;		//reference white
    const double Zr = 1.088754;	//reference white
    double xr,yr,zr;
    double fx, fy, fz;
    double lval,aval,bval;
    
    for(i = 0; i < sz; i++)
    {
        sR = rin[i]; sG = gin[i]; sB = bin[i];
        R = sR/255.0;
        G = sG/255.0;
        B = sB/255.0;
        
        if(R <= 0.04045)	r = R/12.92;
        else				r = pow((R+0.055)/1.055,2.4);
        if(G <= 0.04045)	g = G/12.92;
        else				g = pow((G+0.055)/1.055,2.4);
        if(B <= 0.04045)	b = B/12.92;
        else				b = pow((B+0.055)/1.055,2.4);
        
        X = r*0.4124564 + g*0.3575761 + b*0.1804375;
        Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
        Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
        
        //------------------------
        // XYZ to LAB conversion
        //------------------------
        xr = X/Xr;
        yr = Y/Yr;
        zr = Z/Zr;
        
        if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
        else				fx = (kappa*xr + 16.0)/116.0;
        if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
        else				fy = (kappa*yr + 16.0)/116.0;
        if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
        else				fz = (kappa*zr + 16.0)/116.0;
        
        lval = 116.0*fy-16.0;
        aval = 500.0*(fx-fy);
        bval = 200.0*(fy-fz);
        
        lvec[i] = lval; avec[i] = aval; bvec[i] = bval;
    }
}

void getLABXYSeeds(int STEP, int width, int height, int* seedIndices, int* numseeds)
{
    const bool hexgrid = false;
	int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip,yerrperstrip;
    int xoff,yoff;
    int x,y;
    int xe,ye;
    int seedx,seedy;
    int i;

	xstrips = (0.5+(double)(width)/(double)(STEP));
	ystrips = (0.5+(double)(height)/(double)(STEP));
    
    xerr = width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = width - STEP*xstrips;}
    yerr = height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = height- STEP*ystrips;}
    
	xerrperstrip = (double)(xerr)/(double)(xstrips);
	yerrperstrip = (double)(yerr)/(double)(ystrips);
    
	xoff = STEP/2;
	yoff = STEP/2;
    
    n = 0;
	for( y = 0; y < ystrips; y++ )
	{
		ye = y*yerrperstrip;
		for( x = 0; x < xstrips; x++ )
		{
			xe = x*xerrperstrip;
            seedx = (x*STEP+xoff+xe);
            if(hexgrid){ seedx = x*STEP+(xoff<<(y&0x1))+xe; if(seedx >= width)seedx = width-1; }//for hex grid sampling
            seedy = (y*STEP+yoff+ye);
            i = seedy*width + seedx;
			seedIndices[n] = i;
			n++;
		}
	}
    *numseeds = n;
}

void PerformSuperpixelSLIC(double* lvec, double* avec, double* bvec, double* kseedsl, double* kseedsa, double* kseedsb, double* kseedsx, double* kseedsy, int width, int height, int numseeds, int* klabels, int STEP, double compactness)
{
    int x1, y1, x2, y2;
	double l, a, b;
	double dist;
	double distxy;
    int itr;
    int n;
    int x,y;
    int i;
    int ind;
    int r,c;
    int k;
    int sz = width*height;
	const int numk = numseeds;
	int offset = STEP;
    
    double* clustersize = mxMalloc(sizeof(double)*numk);
    double* inv         = mxMalloc(sizeof(double)*numk);
    double* sigmal      = mxMalloc(sizeof(double)*numk);
    double* sigmaa      = mxMalloc(sizeof(double)*numk);
    double* sigmab      = mxMalloc(sizeof(double)*numk);
    double* sigmax      = mxMalloc(sizeof(double)*numk);
    double* sigmay      = mxMalloc(sizeof(double)*numk);
    double* distvec     = mxMalloc(sizeof(double)*sz);
	double invwt = 1.0/((STEP/compactness)*(STEP/compactness));
    
	for( itr = 0; itr < 10; itr++ )
	{
		for(i = 0; i < sz; i++){distvec[i] = DBL_MAX;}
     
		for( n = 0; n < numk; n++ )
		{
            x1 = kseedsx[n]-offset; if(x1 < 0) x1 = 0;
            y1 = kseedsy[n]-offset; if(y1 < 0) y1 = 0;
            x2 = kseedsx[n]+offset; if(x2 > width)  x2 = width;
            y2 = kseedsy[n]+offset; if(y2 > height) y2 = height;
            
			for( y = y1; y < y2; y++ )
			{
				for( x = x1; x < x2; x++ )
				{
					i = y*width + x;
                    
					l = lvec[i];
					a = avec[i];
					b = bvec[i];
                    
					dist =			(l - kseedsl[n])*(l - kseedsl[n]) +
                                    (a - kseedsa[n])*(a - kseedsa[n]) +
                                    (b - kseedsb[n])*(b - kseedsb[n]);
                    
					distxy =		(x - kseedsx[n])*(x - kseedsx[n]) + (y - kseedsy[n])*(y - kseedsy[n]);
					
					dist += distxy*invwt;
                    
					if(dist < distvec[i])
					{
						distvec[i] = dist;
						klabels[i]  = n;
					}
				}
			}
		}
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
        for(k = 0; k < numk; k++)
        {
            sigmal[k] = 0;
            sigmaa[k] = 0;
            sigmab[k] = 0;
            sigmax[k] = 0;
            sigmay[k] = 0;
            clustersize[k] = 0;
        }
        
		ind = 0;
        for( r = 0; r < height; r++ )
        {
            for( c = 0; c < width; c++ )
            {
                if(klabels[ind] >= 0)
                {
                    sigmal[klabels[ind]] += lvec[ind];
                    sigmaa[klabels[ind]] += avec[ind];
                    sigmab[klabels[ind]] += bvec[ind];
                    sigmax[klabels[ind]] += c;
                    sigmay[klabels[ind]] += r;
                    clustersize[klabels[ind]] += 1.0;
                }
                ind++;
            }
        }
        
		{for( k = 0; k < numk; k++ )
		{
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
		}}
		
		{for( k = 0; k < numk; k++ )
		{
			kseedsl[k] = sigmal[k]*inv[k];
			kseedsa[k] = sigmaa[k]*inv[k];
			kseedsb[k] = sigmab[k]*inv[k];
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
		}}
	}
    mxFree(sigmal);
    mxFree(sigmaa);
    mxFree(sigmab);
    mxFree(sigmax);
    mxFree(sigmay);
    mxFree(clustersize);
    mxFree(inv);
    mxFree(distvec);
}

void EnforceSuperpixelConnectivity(int* labels, int width, int height, int numSuperpixels,int* nlabels, int* finalNumberOfLabels)
{
    int i,j,k;
    int n,c,count;
    int x,y;
    int ind;
    int oindex, adjlabel;
    int label;
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    const int sz = width*height;
    const int SUPSZ = sz/numSuperpixels;
    int* xvec = mxMalloc(sizeof(int)*SUPSZ*10);
	int* yvec = mxMalloc(sizeof(int)*SUPSZ*10);

	for( i = 0; i < sz; i++ ) nlabels[i] = -1;
    oindex = 0;
    adjlabel = 0;//adjacent label
    label = 0;
	for( j = 0; j < height; j++ )
	{
		for( k = 0; k < width; k++ )
		{
			if( 0 > nlabels[oindex] )
			{
				nlabels[oindex] = label;
				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				{for( n = 0; n < 4; n++ )
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					if( (x >= 0 && x < width) && (y >= 0 && y < height) )
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}}
                
				count = 1;
				for( c = 0; c < count; c++ )
				{
					for( n = 0; n < 4; n++ )
					{
						x = xvec[c] + dx4[n];
						y = yvec[c] + dy4[n];
                        
						if( (x >= 0 && x < width) && (y >= 0 && y < height) )
						{
							int nindex = y*width + x;
                            
							if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}
                        
					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= SUPSZ >> 2)
				{
					for( c = 0; c < count; c++ )
					{
                        ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}
			oindex++;
		}
	}
	*finalNumberOfLabels = label;
    
	mxFree(xvec);
	mxFree(yvec);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int width;
    int height;
    int sz;
    int i, ii;
    int x, y;
    int* rin; int* gin; int* bin;
    int* klabels;
    int* clabels;
    double* lvec; double* avec; double* bvec;
    int step;
    int* seedIndices;
    int numseeds;
    double* kseedsx;double* kseedsy;
    double* kseedsl;double* kseedsa;double* kseedsb;
    int k;
    const mwSize* dims;//int* dims;
    int* outputNumSuperpixels;
    int* outlabels;
    int finalNumberOfLabels;
    unsigned char* imgbytes;
    int numelements;
    int numSuperpixels = 200;//default value
    double compactness = 10;//default value

    if (nrhs < 1) {
        mexErrMsgTxt("At least one argument is required.") ;
    } else if(nrhs > 3) {
        mexErrMsgTxt("Too many input arguments.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("SLIC:nlhs","Two outputs required, a labels and the number of labels, i.e superpixels.");
    }
    //---------------------------
    numelements   = mxGetNumberOfElements(prhs[0]) ;
    mwSize numdims = mxGetNumberOfDimensions(prhs[0]) ;
    dims  = mxGetDimensions(prhs[0]) ;
    imgbytes  = (unsigned char*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it
    width = dims[1]; height = dims[0];//Note: first dimension provided is height and second is width
    sz = width*height;
    //---------------------------
    numSuperpixels  = mxGetScalar(prhs[1]);
    compactness     = mxGetScalar(prhs[2]);
    
    //---------------------------
    // Allocate memory
    //---------------------------
    rin    = mxMalloc( sizeof(int)      * sz ) ;
    gin    = mxMalloc( sizeof(int)      * sz ) ;
    bin    = mxMalloc( sizeof(int)      * sz ) ;
    lvec    = mxMalloc( sizeof(double)      * sz ) ;
    avec    = mxMalloc( sizeof(double)      * sz ) ;
    bvec    = mxMalloc( sizeof(double)      * sz ) ;
    klabels = mxMalloc( sizeof(int)         * sz );//original k-means labels
    clabels = mxMalloc( sizeof(int)         * sz );//corrected labels after enforcing connectivity
    seedIndices = mxMalloc( sizeof(int)     * sz );
    
    //---------------------------
    // Perform color conversion
    //---------------------------
    //if(2 == numdims)
    if(numelements/sz == 1)//if it is a grayscale image, copy the values directly into the lab vectors
    {
        for(x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(y = 0; y < height; y++)
            {
                i = y*width+x;
                lvec[i] = imgbytes[ii];
                avec[i] = imgbytes[ii];
                bvec[i] = imgbytes[ii];
                ii++;
            }
        }
    }
    else//else covert from rgb to lab
    {
        for(x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(y = 0; y < height; y++)
            {
                i = y*width+x;
                rin[i] = imgbytes[ii];
                gin[i] = imgbytes[ii+sz];
                bin[i] = imgbytes[ii+sz+sz];
                ii++;
            }
        }
        rgbtolab(rin,gin,bin,sz,lvec,avec,bvec);
    }
    //---------------------------
    // Find seeds
    //---------------------------
    step = sqrt((double)(sz)/(double)(numSuperpixels))+0.5;
    getLABXYSeeds(step,width,height,seedIndices,&numseeds);
    
    kseedsx    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsy    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsl    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsa    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsb    = mxMalloc( sizeof(double)      * numseeds ) ;
    for(k = 0; k < numseeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width;
        kseedsy[k] = seedIndices[k]/width;
        kseedsl[k] = lvec[seedIndices[k]];
        kseedsa[k] = avec[seedIndices[k]];
        kseedsb[k] = bvec[seedIndices[k]];
    }
    //---------------------------
    // Compute superpixels
    //---------------------------
    PerformSuperpixelSLIC(lvec, avec, bvec, kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,width,height,numseeds,klabels,step,compactness);
    //---------------------------
    // Enforce connectivity
    //---------------------------
    EnforceSuperpixelConnectivity(klabels,width,height,numSuperpixels,clabels,&finalNumberOfLabels);
    //---------------------------
    // Assign output labels
    //---------------------------
    plhs[0] = mxCreateNumericMatrix(height,width,mxINT32_CLASS,mxREAL);
    outlabels = mxGetData(plhs[0]);
    for(x = 0, ii = 0; x < width; x++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(y = 0; y < height; y++)
        {
            i = y*width+x;
            outlabels[ii] = clabels[i];
            ii++;
        }
    }
    //---------------------------
    // Assign number of labels/seeds
    //---------------------------
    plhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    outputNumSuperpixels = (int*)mxGetData(plhs[1]);//gives a void*, cast it to int*
    *outputNumSuperpixels = finalNumberOfLabels;
    //---------------------------
    // Deallocate memory
    //---------------------------
    mxFree(rin);
    mxFree(gin);
    mxFree(bin);
    mxFree(lvec);
    mxFree(avec);
    mxFree(bvec);
    mxFree(klabels);
    mxFree(clabels);
    mxFree(seedIndices);
    mxFree(kseedsx);
    mxFree(kseedsy);
    mxFree(kseedsl);
    mxFree(kseedsa);
    mxFree(kseedsb);
}
