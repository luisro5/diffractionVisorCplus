/**
 *
 *
 * This source code is licensed under the MIT license found in the
 * MIT-LICENSE file in the root directory of this source tree.
 *
 *
 */
#include <iostream>
#include <vector>
#include <iomanip>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cmath>
using namespace std;
const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
//********************************PROTOTYPES****************/
void fft(CArray&);
void ifft(CArray&);
void initial(int[],int);
void transposeShift(Complex[],Complex[],int);
void fftshift2d(Complex[],Complex[],int);
void x0Filter(int[],int,int,int);
void generatePlate(double[],double[],int,double[]);
double erf(int,int,double,int,int);
void identityMatrix(int,Complex[]);
double getSizePlate(double[],double[],int);
void atomicPotential80(int[],int[],int,double,double,double,std::string,
	double,Complex[]);
Complex normalDivide(Complex,Complex);
void d2FFT(Complex[],Complex[],int,int);
void fresnelnu(Complex[],Complex[],int,int,int,double,double);
void convolve2D(Complex[],Complex[],int,int,Complex[],int,int);
//********************************PROTOTYPES****************/
// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
// void fft(CArray& x)
// {
//     const size_t N = x.size();
//     if (N <= 1) return;

//     // divide
//     CArray even = x[std::slice(0, N/2, 2)];
//     CArray  odd = x[std::slice(1, N/2, 2)];

//     // conquer
//     fft(even);
//     fft(odd);

//     // combine
//     for (size_t k = 0; k < N/2; ++k)
//     {
//         Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
//         x[k    ] = even[k] + t;
//         x[k+N/2] = even[k] - t;
//     }
// }

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
void fft(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0;void atomicPotential80(int l[],int m[],int side,double currentZ,double secondLastZ,
	double lastZ,string atomicNumber,double lam,Complex flatChem[]) l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
	//// Normalize (This section make it not working correctly)
	//Complex f = 1.0 / sqrt(N);
	//for (unsigned int i = 0; i < N; i++)
	//	x[i] *= f;
}

// inverse fft (in-place)
void ifft(CArray &x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);

    // forward fft
    fft( x );

    // conjugate the complex numbers again
    x = x.apply(std::conj);

    // scale the numbers
    x /= x.size();
}
void initial(int base[],int totalS)
{
	int k=0;
	while (k < totalS)
	{
		base[k]=0;
		k++;
	}
}
//this function must be applied after the FFT
//also asume that n*n matrix as inputs
void transposeShift(Complex original[],Complex shifted[],int r)
{
  int xshift=r/2;
  int yshift=r/2;
	int i,j;
  for ( i = 0; i < r; i++)
  {
    for (j = 0; j < r; j++)
    {
      shifted[i * r + j] = original[j * r + i];
    }
  }
}
//Asume that is n*n square
void fftshift2d(Complex original[],Complex shifted[],int r)
{
	int xdim=r,ydim=r,xshift,yshift,x,y,outx,outy;
	xshift=xdim/2;
	yshift=ydim/2;
	for (x = 0; x < xdim; x++)
	{
		outx=(x+xshift)%xdim;
		for (y = 0; y < ydim; y++)
		{
			outy=(y+yshift)%ydim;
			shifted[outx+xdim*outy]=original[x+xdim*y];
		}
	}
}
// n and m must be the same !!
void x0Filter(int out[],int m,int n,int dx0)
{
    int i=-n/2,da,counter,counter2;
		counter2=0;
    while(i<(n/2))
    {
      da=i*dx0;
      counter=0;
      while(counter<m)
      {
        out[counter2]=da;
				counter2++;
        counter++;
      }
      i++;
    }
}
//this function process each of the section of the array 1 PLANE AT TIME
void generatePlate(double theX[],double theY[],int side,double flat[])
{
	int square2,indexIm,ix,iy,position;
	double center;
  square2=pow(side,2);
  center=floor(side/2);
  for(indexIm=0;indexIm<square2;indexIm++)
  {
    ix=center+theX[indexIm];
    ky=center+theY[indexIm];
    position=(side*ix)+ky;
    flat[position]=1;
  }
}
// dt and dx must be 1 (discrete values)
double erf(int ll,int ul,double t,int dt, int dx)
{
	int i;
	double valErf,errorAcum;
	valErf = 0;
	for(i = ll; i <= ul; i+=dx)
	{
		valErf+=exp(-pow(t,2))*dt;
	}
	errorAcum=(2/sqrt(PI))*valErf;
	return errorAcum;
}
//identity Matrix
void identityMatrix(int side,Complex baseArray[])
{
	Complex theVoid=Complex(0,0);
	Complex identity=Complex(1,0);
	int x,y,counterIdentity;
	counterIdentity=0;
	for(y=0;y<side;y++)
	{
		for(x=0;x<side;x++)
		{
			if(x==y)
			{
				baseArray[counterIdentity]=identity;
			}
			counterIdentity++;
		}
	}
}
//We asume that x and y array are same size with the same Z
double getSizePlate(double theX[],double theY[],int all1D)
{
	double* per1,per2,biggy;
	per1=std::max_element(theX, theX + all1D);
	per2=std::max_element(theY, theY + all1D);
	if(*per1>*per2)
	{
		biggy=per1;
	}
	else
	{
		biggy=per2;
	}
	//by 4 it is only for presentation
  return biggy*4;
}
//****************************************************************
//IMPORTANT: remember to add the complex library to return the
//array of complex as index of Chemical bright to each atom
//in 2 DIMENSION
//****************************************************************
//atomicNumber,element,ya1,yb1,ya2,yb2,ya3,yb3,ya4,yb4,ya5,yb5
//this function is in adition with the initial set,
//to get the electro Chemical Potential,before the visual process
void atomicPotential80(int l[],int m[],int side,double currentZ,double secondLastZ,
	double lastZ,std::string atomicNumber,double lam,Complex flatChem[])
	{
		int centerAtomicPot;
		centerAtomicPot=floor(side/2);
		double ai[5] = {0,0,0,0,0};
		double bi[5] = {0,0,0,0,0};
		/**********************TERMINAR!!!!!!!!!**************************/
		//$atomSpecs=getAtomScatt($atomicNumber,'scatteringTables40.txt');
		/****************************************************************/
		//Ti
		ai[0]=2.6865;
	  ai[1]=1.2748;
	  ai[2]=0.3853;
	  ai[3]=1.0297;
	  ai[4]=1.9382;
	  bi[0]=16.7493;
	  bi[1]=4.7274;
	  bi[2]=0.2125;
	  bi[3]=1.3530;
	  bi[4]=50.8067;
		/*
		ai[0]=atomSpecs[2];
		ai[1]=atomSpecs[3];
		ai[2]=atomSpecs[4];
		ai[3]=atomSpecs[5];
		ai[4]=atomSpecs[6];
		bi[0]=atomSpecs[7];
		bi[1]=atomSpecs[8];
		bi[2]=atomSpecs[9];
		bi[3]=atomSpecs[10];
		bi[4]=atomSpecs[11];
		*/
	  double planck2=4.390477714e-67,massElectron=9.10938188e-31;
	  double preb=9.869604401,first,second,third,fourth,fifth,sixth,pseudoPos;
	  double particleBand=(-lam)/((2*PI*massElectron*lam)/planck2);
		int indexIm1,indexIm2,alpha,coeff,position;
		double electronicPotential;
	  for(indexIm1=0;indexIm1<side;indexIm1++)
	  {
	    for(indexIm2=0;indexIm2<side;indexIm2++)
	    {
	        alpha=pow(l[indexIm1],2)+pow(m[indexIm2],2);
	        for(coeff=0;coeff<5;coeff++)
	        {
	            first=(ai[coeff]/bi[coeff]);
	            second=(-preb*alpha)/bi[coeff];
	            third=sqrt((2*preb)/bi[coeff]);
	            fourth=erf(currentZ,lastZ,third,1,1);
	            fifth=erf(currentZ,secondLastZ,third,1,1);
	            sixth=fourth-fifth;
	            pseudoPos+=first*exp(second)*sixth;
	        }
	        electronicPotential=particleBand*pseudoPos;
	        position=(centerAtomicPot*l[indexIm1])+m[indexIm2];
	        flatChem[position]=Complex(electronicPotential,0);
	    }
	  }
	}
	//in the form a/b
	Complex normalDivide(Complex a,Complex b)
	{
		double acom,bcom,ccom,dcom,denom,anumer,bnumer,semA,semB;
		Complex cutter;
		acom=real(a);
		bcom=imag(a);
		ccom=real(b);
		dcom=imag(b);
		denom=pow(ccom,2)+pow(dcom,2);
		anumer=(acom*ccom)+(bcom*dcom);
		bnumer=(bcom*ccom)-(acom*dcom);
		semA = anumer/denom;
		semB = bnumer/denom;
		cutter = Complex(semA,semB);
		return cutter;
	}
	void d2FFT(Complex allsquare[],Complex squarePhase[],int rows,int columns)
	{
		int moveRow,moveColumns;
		int allunits=rows*columns;
		Complex miniA[]={};
		Complex miniB[]={};
		Complex miniC[]={};
	  //first we process the $rows
	  for(moveRow = 0; moveRow < allunits; moveRow+=rows)
	  {
	    miniA=std::slice(allsquare,moveRow,rows);
	    fft(miniA);
			copy(miniA, miniA + moveRow,allsquare+moveRow);
	  }
	  transposeShift(miniA,miniB,rows);
	  for(moveColumns = 0; moveColumns < allunits; moveColumns+=columns)
	  {
	    miniC=std::slice(allsquare,moveColumns,columns);
	    fft(miniC);
			copy(miniC, miniC + moveColumns,allsquare+moveColumns);
	  }
		int normal;
		double acom,bcom,ccom,dcom,denom,anumer,bnumer,semA,semB;
	  Complex normfft=Complex((double)allunits,0),pointfftA,pointfftB;
	  for(normal=0;normal<allunits;normal++)
	  {
			pointfftA=allsquare[normal];
			squarePhase[normal]=normalDivide(pointfftA,normfft);
	  }
	}
	void fresnelnu(Complex pref0[],Complex lastKernel[],int M,int N,int dx0=1,double z,double lamb)
	{
	  // NOTE: fo and m*n must be the same size
	  if(z<0.1)
	  {
	      z=1;
	  }
	  double k=(2*PI)/lamb;
		double X0adjusted1[]={},X0adjusted2[]={};
	  //the x0 optical filter
	  x0Filter(X0adjusted1,M,N,dx0);
	  Complex normalize=Complex(lamb*z,0),t,t2;
		Complex secondKernel[]={},thirdKernel[]={};
	  //the y0 optical filter
	  transposeShift(X0adjusted1,X0adjusted2,N);
	  //g=f0.*exp(j*0.5*k*(x0.^2+y0.^2)/z); this expression we can handle as complex
	  int preCharge,allprf0,charge;
		double firstKernel,kernel_0,kernel_1;
		allprf0=M*N;
	  for(preCharge=0;preCharge<allprf0;preCharge++)
	  {
	    firstKernel=0.5*k*(pow(X0adjusted1[preCharge],2)+pow(X0adjusted2[preCharge],2))/z;
	    kernel_0=cos(firstKernel);
	    kernel_1=-sin(firstKernel);
	    t = Complex(kernel_0, kernel_1);
	    t2 = pref0[preCharge]*t;
	    secondKernel[preCharge]=t2;
	  }
	  //G=fftshift(fft2(g));
	  d2FFT(secondKernel,thirdKernel,M,N);
		Complex fourthKernel[]={};
	  fftshift2d(thirdKernel,fourthKernel,N);
	  // we get the adjust of dx1
	  double du=1/(M*dx0);
	  double dx1=lamb*z*du;
		double X0adjusted3[]={},X0adjusted4[]={};
		double fifthKernel,kernel_2,kernel_3;
	  //the x1 optical filter
	  x0Filter(X0adjusted3,M,N,dx1);
	  //the y1 optical filter
		transposeShift(X0adjusted3,X0adjusted4,N);
	  //f1=G.*exp(i*0.5*k*(x1.^2+y1.^2)/z);
		Complex tk,t3,t4;
		Complex normfresnel=Complex((double)allprf0,0);
	  for(charge=0;charge<allprf0;charge++)
	  {
	    fifthKernel=0.5*k*(pow(X0adjusted3[charge],2)+pow(X0adjusted4[charge],2))/z;
	    kernel_2=cos(fifthKernel);
	    kernel_3=-sin(fifthKernel);
	    tk = Complex(kernel_2, kernel_3);
	    t3=fourthKernel[charge]*tk;
	    t4=normalDivide(t4,normfresnel);
	    lastKernel[charge]=t4;
	  }
	}
	// utitlities for matrix************
	//******NOTICE:must be same size both arguments!! padd with zeros and change panels as fftShift
	//after
	/*
	function convolution2($thex,$they,$side)
	{
	  $signalx=FFT($thex);
	  $signaly=FFT($they);
	  $both=multiplyArrComplex($signaly,$signalx);
	  $smooth=IFFT($both);
	  //need normalize vector
	  $smooth2=normalizeComplex($smooth);
	  return $smooth2;
	}
	*/
	void convolve2D( Complex in[], Complex out[], int dataSizeX,
	  int dataSizeY,Complex kernel[], int kernelSizeX, int kernelSizeY)
	{
	    int i, j, m, n, mm, nn;
	    int kCenterX, kCenterY;                         // center index of kernel
	    Complex sum=Complex(0,0);                                      // temp accumulation buffer
	    Complex temporal8=Complex(0,0);
	    Complex temporal9=Complex(0,0);
	    int rowIndex, colIndex;

	    // find center position of kernel (half of kernel size)
	    kCenterX = kernelSizeX / 2;
	    kCenterY = kernelSizeY / 2;

	    for(i=0; i < dataSizeY; ++i)                // rows
	    {
	        for(j=0; j < dataSizeX; ++j)            // columns
	        {
	            sum = Complex(0,0);                            // init to 0 before sum
	            for(m=0; m < kernelSizeY; ++m)      // kernel rows
	            {
	                mm = kernelSizeY - 1 - m;       // row index of flipped kernel

	                for(n=0; n < kernelSizeX; ++n)  // kernel columns
	                {
	                    nn = kernelSizeX - 1 - n;   // column index of flipped kernel

	                    // index of input signal, used for checking boundary
	                    rowIndex = i + (kCenterY-mm);
	                    colIndex = j + (kCenterX-nn);

	                    // ignore input samples which are out of bound
	                    if(rowIndex >= 0 && rowIndex < dataSizeY && colIndex >= 0 && colIndex < dataSizeX)
	                    {
	                      temporal8=in[dataSizeX * rowIndex + colIndex];
	                      temporal9=kernel[kernelSizeX * mm + nn];
	                      sum+=temporal8+temporal9;
	                    }
	                }
	            }
	            out[dataSizeX * i + j] = sum;
	        }
	    }
	}
