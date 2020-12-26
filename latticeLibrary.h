/**
 *
 *
 * This source code is licensed under the MIT license found in the
 * MIT-LICENSE file in the root directory of this source tree.
 *
 *
 */
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <valarray>
#include <cmath>
//********************************PROTOTYPES****************/
function next_pow(int toX);
double convertToRadians(double);
double findLambda(double);
void rotativeMatrix(double,double,double,float[]);
void lattice4(float,float,float,float,float,float,int,int,int,
  float,float,float,double[],double[],double[]);
void bases(float,float,float,int,int,int,float,float,float,
  double[],double[],double[]);
void inner(float,float,float,int,int,int,float,float,
  float,double[],double[],double[]);
void face(float,float,float,int,int,int,float,float,
  float,double[],double[],double[]);
void layersInZ(double[],double[],int);
//********************************PROTOTYPES****************/
//next power of 2:
function next_pow(int toX)
{
  int power = 1;
  while(power < toX)
  {
    power*=2;
  }
  return power;
}
// Function for conversion
double convertToRadians(double degree)
{
    double pi = 3.14159265359;
    return (degree * (pi / 180));
}
//this function find lambda based on the voltage
double findLambda(double volts)
{
  double subterm,term1,term2,newlambda,h,emass,squarelight,elementalCharge;
  double newlambda;
  h=6.62606957e-34;
  emass= 9.10938291e-31;
  squarelight=8.987551787e+16;
  elementalCharge=1.60217657e-19;
  subterm=(elementalCharge*volts)/(2*emass*squarelight);
  term1=h/sqrt(2*emass*elementalCharge*$volts);
  term2=1/sqrt(1+subterm);
  newlambda=term1*term2;
  return newlambda;
}
//rotate the matrix using rotation angles
void rotativeMatrix(double phi,double tetha,double psi,float results[])
{
  double phi2,tetha2,psi2,r11,r12,r13,r21,r22,r23,r31,r32,r33;
  phi2=convertToRadians(phi);
  tetha2=convertToRadians(tetha);
  psi2=convertToRadians(psi);
  r11=(cos(psi2)*cos(tetha2)*cos(phi2))-(sin(psi2)*sin(phi2));
  r12=(-cos(psi2)*cos(tetha2)*sin(phi2))-(sin(psi2)*sin(phi2));
  r13=cos(psi2)*sin(tetha2);
  r21=(sin(psi2)*cos(tetha2)*cos(phi2))+(cos(psi2)*sin(phi2));
  r22=(-sin(psi2)*cos(tetha2)*sin(phi2))+(cos(psi2)*cos(phi2));
  r23=sin(psi2)*sin(tetha2);
  r31=-sin(tetha2)*cos(phi2);
  r32=sin(tetha2)*sin(phi2);
  r33=cos(tetha2);
  results[0]=r11;
  results[1]=r12;
  results[2]=r13;
  results[3]=r21;
  results[4]=r22;
  results[5]=r23;
  results[6]=r31;
  results[7]=r32;
  results[8]=r33;
}
// generates the lattice points to generate the matrix of points
function lattice4(float aa,float bb,float cc,float alfie,float betie,float gammie,
  int l,int m,int n,float psi4,float phi4,float tetha4,double ax[],
  double ay[],double az[])
{
  double alti5,beti5,gati5,doublealfa,doublebeta,doublegamma;
  double prevol1,prevol2,vol,i,j,k,bas1,bas2,bas3;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double b11,b12,b13,b21,b22,b23,b31,b32,b33;
  double c11,c12,c13,c21,c22,c23,c31,c32,c33;
  int precounter=0;
  float rotates[9] = {};
  alti5=convertToRadians(alfie);
  beti5=convertToRadians(betie);
  gati5=convertToRadians(gammie);
  doublealfa=pow(cos(alti5),2);
  doublebeta=pow(cos(beti5),2);
  doublegamma=pow(cos(gati5),2);
  prevol1=aa*bb*cc;
  prevol2=2*cos(alti5)*cos(beti5)*cos(gati5);
  vol=prevol1*sqrt(1-doublealfa-doublebeta-doublegamma+prevol2);
  a11=aa;
  a12=bb*cos(gati5);
  a13=cc*cos(beti5);
  a21=0;
  a22=bb*sin(gati5);
  a23=cc*((cos(alti5)-cos(beti5)*cos(gati5))/sin(gati5));
  a31=0;
  a32=0;
  a33=vol/(aa*bb*sin(gati5));
  rotativeMatrix(phi4,tetha4,psi4,rotates);
  b11=rotates[0];
  b12=rotates[1];
  b13=rotates[2];
  b21=rotates[3];
  b22=rotates[4];
  b23=rotates[5];
  b31=rotates[6];
  b32=rotates[7];
  b33=rotates[8];
  c11=(b11*a11)+(b12*a21)+(b13*a31);
  c12=(b11*a12)+(b12*a22)+(b13*a32);
  c13=(b11*a13)+(b12*a23)+(b13*a33);
  c21=(b21*a11)+(b22*a21)+(b23*a31);
  c22=(b21*a12)+(b22*a22)+(b23*a32);
  c23=(b21*a13)+(b22*a23)+(b23*a33);
  c31=(b31*a11)+(b32*a21)+(b33*a31);
  c32=(b31*a12)+(b32*a22)+(b33*a32);
  c33=(b31*a13)+(b32*a23)+(b33*a33);
  for(i=-l;i<l;i++)
  {
    for(j=-m;j<m;j++)
    {
      for(k=-n;k<n;k++)
      {
        bas1=(j*c11)+(i*c12)+(k*c13);
        bas2=(j*c21)+(i*c22)+(k*c23);
        bas3=(j*c31)+(i*c32)+(k*c33);
        ax[precounter]=bas1;
        ay[precounter]=bas2;
        az[precounter]=bas3;
        precounter++;
      }
    }
  }
}
//base crystallography structure
void bases(float aa,float bb,float cc,int l,int m,int n,float psi4,float phi4,
float tetha4,double ax[],double ay[],double az[])
{
  float rotates[9] = {};
  double a11,a12,a13,a22,a23,a31,a32,a33;
  double b11,b12,b13,b21,b22,b23,b31,b32,b33;
  double c11,c12,c13,c21,c22,c23,c31,c32,c33;
  int i,j,k;
  a11=aa;
  a12=0;
  a13=0;
  a21=0;
  a22=bb;
  a23=0;
  a31=0;
  a32=0;
  a33=cc;
  rotativeMatrix(phi4,tetha4,psi4,rotates);
  b11=rotates[0];
  b12=rotates[1];
  b13=rotates[2];
  b21=rotates[3];
  b22=rotates[4];
  b23=rotates[5];
  b31=rotates[6];
  b32=rotates[7];
  b33=rotates[8];
  c11=b11*a11;
  c12=b12*a22;
  c13=b13*a33;
  c21=b21*a11;
  c22=b22*a22;
  c23=b23*a33;
  c31=b31*a11;
  c32=b32*a22;
  c33=b33*a33;
  double aaa,bbb,ccc,mon1,mon2;
  int precounter2=0;
  for(i=-l;i<l;i++)
  {
    for(j=-m;j<m;j++)
    {
      for(k=-n;k<n;k++)
      {
        aaa=(j*c11)+(i*c12)+(k*c13);
        bbb=(j*c21)+(i*c22)+(k*c23);
        ccc=(j*c31)+(i*c32)+(k*c33);
        mon2=(i+0.5)*$bbb;
        mon1=(j+0.5)*$aaa;
        ax[precounter2]=j*aaa;
        ay[precounter2]=i*bbb;
        az[precounter2]=k*ccc;
        ax[precounter2+(m)]=mon1;
        ay[precounter2+(l)]=mon2;
        az[precounter2+(n)]=k*ccc;
        precounter2++;
      }
    }
  }
}
//inner crystallography structure
void inner(float aa,float bb,float cc,int l,int m,int n,float psi4,float phi4,
float tetha4,double ax[],double ay[],double az[])
{
  float rotates[9] = {};
  double a11,a12,a13,a22,a23,a31,a32,a33;
  double b11,b12,b13,b21,b22,b23,b31,b32,b33;
  double c11,c12,c13,c21,c22,c23,c31,c32,c33;
  int i,j,k;
  a11=aa;
  a12=0;
  a13=0;
  a21=0;
  a22=bb;
  a23=0;
  a31=0;
  a32=0;
  a33=cc;
  rotativeMatrix(phi4,tetha4,psi4,rotates);
  b11=rotates[0];
  b12=rotates[1];
  b13=rotates[2];
  b21=rotates[3];
  b22=rotates[4];
  b23=rotates[5];
  b31=rotates[6];
  b32=rotates[7];
  b33=rotates[8];
  c11=b11*a11;
  c12=b12*a22;
  c13=b13*a33;
  c21=b21*a11;
  c22=b22*a22;
  c23=b23*a33;
  c31=b31*a11;
  c32=b32*a22;
  c33=b33*a33;
  double aaa,bbb,ccc,mon1,mon2,mon3;
  int precounter3=0;
  for(i=-l;i<l;i++)
  {
    for(j=-m;j<m;j++)
    {
      for(k=-n;k<n;k++)
      {
        aaa=(j*c11)+(i*c12)+(k*c13);
        bbb=(j*c21)+(i*c22)+(k*c23);
        ccc=(j*c31)+(i*c32)+(k*c33);
        mon2=(i+0.5)*bbb;
        mon1=(j+0.5)*aaa;
        mon3=(k+0.5)*ccc;
        ax[precounter3]=j*aaa;
        ay[precounter3]=i*bbb;
        az[precounter3]=k*ccc;
        ax[precounter3+(m)]=mon1;
        ay[precounter3+(l)]=mon2;
        az[precounter3+(n)]=mon3;
        precounter3++;
      }
    }
  }
}
//face crystallography structure
void face(float aa,float bb,float cc,int l,int m,int n,float psi4,float phi4,
float tetha4,double ax[],double ay[],double az[])
{
  float rotates[9] = {};
  double a11,a12,a13,a22,a23,a31,a32,a33;
  double b11,b12,b13,b21,b22,b23,b31,b32,b33;
  double c11,c12,c13,c21,c22,c23,c31,c32,c33;
  int i,j,k;
  a11=aa;
  a12=0;
  a13=0;
  a21=0;
  a22=bb;
  a23=0;
  a31=0;
  a32=0;
  a33=cc;
  rotativeMatrix(phi4,tetha4,psi4,rotates);
  b11=rotates[0];
  b12=rotates[1];
  b13=rotates[2];
  b21=rotates[3];
  b22=rotates[4];
  b23=rotates[5];
  b31=rotates[6];
  b32=rotates[7];
  b33=rotates[8];
  c11=b11*a11;
  c12=b12*a22;
  c13=b13*a33;
  c21=b21*a11;
  c22=b22*a22;
  c23=b23*a33;
  c31=b31*a11;
  c32=b32*a22;
  c33=b33*a33;
  double aaa,bbb,ccc,mon1,mon2,mon3;
  int precounter4=0;
  for(i=-l;i<l;i++)
  {
    for(j=-m;j<m;j++)
    {
      for(k=-n;k<n;k++)
      {
        aaa=(j*c11)+(i*c12)+(k*c13);
        bbb=(j*c21)+(i*c22)+(k*c23);
        ccc=(j*c31)+(i*c32)+(k*c33);
        mon2=(i+0.5)*bbb;
        mon1=(j+0.5)*aaa;
        mon3=(k+0.5)*ccc;
        ax[precounter4]=j*aaa;
        ay[precounter4]=i*bbb;
        az[precounter4]=k*ccc;
        ax[precounter4+(2*m)]=mon1;
        ay[precounter4+(2*l)]=mon2;
        az[precounter4+(2*n)]=k*ccc;
        ax[precounter4+(3*m)]=mon1;
        ay[precounter4+(3*l)]=i*bbb;
        az[precounter4+(3*n)]=mon3;
        ax[precounter4+(4*m)]=j*aaa;
        ay[precounter4+(4*l)]=mon2;
        az[precounter4+(4*n)]=mon3;
        precounter4++;
      }
    }
  }
}
//----------------------------------------------------------------
//-----function to start the cicles over sample
//----------------------------------------------------------------
void layersInZ(double arrayZ[],double miniZ[],int limits2)
{
 $piece=array();
 int layer;
 double lay;
 for(layer=0;layer<limits2;layer++)
 {
   lay=arrayZ[layer];
   if(lay)>0.1)
   {
     miniZ[layer]=lay;
   }
   else
   {
     $miniZ[layer]=0;
   }
 }
}
