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
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <valarray>
#include <cmath>
using namespace std;
std::string const digits = "0123456789abcdefghijklmnopqrstuvwxyz";
//********************************PROTOTYPES****************/
void render8(double,double,double,double,int,int);
std::string percentageToHex(double);
std::string to_base(unsigned long, int);
std::string tag(double,double,double,double,double,double);
//********************************PROTOTYPES****************/
//--------------------RENDER FUNCTIONS-------------------------------------
void render8(double x,double y,double z,double minz,int u,int v)
{
  double temporal;
  temporal=z+minz;
  u=(int)(x*16)/temporal;
  v=(int)(y*16)/temporal;
}
//percentageToHex return a string
//percentage must be 1(100%) to 0
std::string percentageToHex(double percentage)
{
  double maximum,toHex;
  maximum=16777215;
  toHex=maximum*percentage;
  return to_base(toHex,16);
}
std::string to_base(unsigned long num, int base)
{
  if (num == 0)
    return "0";

  std::string result;
  while (num > 0)
  {
    std::ldiv_t temp = std::div(num, (long)base);
    result += digits[temp.rem];
    num = temp.quot;
  }
  std::reverse(result.begin(), result.end());
  return result;
}
std::string tag(double pa,double pb,double pc,double anga,double angb,double angg)
{
  double anga4,angb4,angg4,cosquada,cosquadb,cosquadg;
  double part1v,part2v,part3v,part4v,part5v,volume,arecip;
  std::string leyend="",tok1=" A-1";
  anga4=convertToRadians(anga);
  angb4=convertToRadians(angb);
  angg4=convertToRadians(angg);
  cosquada=pow(cos(anga4),2);
  cosquadb=pow(cos(angb4),2);
  cosquadg=pow(cos(angg4),2);
  part1v=pa*pb*pc;
  part2v=1-cosquada-cosquadb-cosquadg;
  part3v=(2*cos(anga4)*cos(angb4))+(2*cos(anga4)*cos(angg4))
  +(2*cos(angb4)*cos(angg4));
  part4v=part2v+part3v;
  part5v=sqrt(part4v);
  volume=part1v*part5v;
  arecip=(pb*pc*sin(anga4))/$volume;
  std::string varAsString = std::to_string(arecip);
  leyend=varAsString+tok1;
  return leyend;
}
/**************************************************************************/
//*************************************************************************/
// electro-Chemical potential library
/*
void getAtomScatt(string atomicNumber,string scattFile,double atoms[])
{
	$handle = fopen($scattFile,'r');
  $row = 0;
  $arr = array();
  while($line = fgetcsv($handle,1000,";"))
  {
    $arr[] = $line;
  }
  $integrate1=$arr[$atomicNumber][0];
  $integrate2=$arr[$atomicNumber][1];
  $integrate3=$arr[$atomicNumber][2];
  $integrate4=$arr[$atomicNumber][3];
  $integrate5=$arr[$atomicNumber][4];
  $integrate6=$arr[$atomicNumber][5];
  $integrate7=$arr[$atomicNumber][6];
  $integrate8=$arr[$atomicNumber][7];
  $integrate9=$arr[$atomicNumber][8];
  $integrate10=$arr[$atomicNumber][9];
  $integrate11=$arr[$atomicNumber][10];
  $integrate12=$arr[$atomicNumber][11];
  $integrate=array($integrate1,$integrate2,$integrate3,$integrate4,
  $integrate5,$integrate6,$integrate7,$integrate8,$integrate9,
  $integrate10,$integrate11,$integrate12);
  return $integrate;
}
*/
