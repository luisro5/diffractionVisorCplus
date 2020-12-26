/**
 *
 *
 * This source code is licensed under the MIT license found in the
 * MIT-LICENSE file in the root directory of this source tree.
 *
 *

COLUMNS
DATA TYPE
FIELD
DEFINITION
-------------------------------------------------------------------------------------
1 - 6 Record name "ATOM "
7 - 11 Integer serial Atom serial number.
13 - 16 Atom name Atom name.
17 Character altLoc Alternate location indicator.
18 - 20 Residue name resName Residue name.
22 Character chainID Chain identifier.
23 - 26 Integer resSeq Residue sequence number.
27 AChar iCode Code for insertion of residues.
31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms.
39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms.
47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms.
55 - 60 Real(6.2) occupancy Occupancy.
61 - 66 Real(6.2) tempFactor Temperature factor.
77 - 78 LString(2) element Element symbol, right-justified.
79 - 80 LString(2) charge Charge on the atom.
*/
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;
/** Return the current date in a string, formatted as either ISO-8601
 *  or "Weekday-name, Month-name Day, Year".
 *
 *  The date is initialized when the object is created and will return
 *  the same date for the lifetime of the object.  The date returned
 *  is the date in the local timezone.
 */
 class Date
 {
     struct tm ltime;

 public:
     /// Default constructor.
     Date()
     {
         time_t t = time(0);
         localtime_r(&t, &ltime);
     }

     /** Return the date based on a format string.  The format string is
      *  fed directly into strftime().  See the strftime() documentation
      *  for information on the proper construction of format strings.
      *
      *  @param[in] fmt is a valid strftime() format string.
      *
      *  @return a string containing the formatted date, or a blank string
      *      if the format string was invalid or resulted in a string that
      *      exceeded the internal buffer length.
      */
     std::string getDate(const char* fmt)
     {
         char out[200];
         size_t result = strftime(out, sizeof out, fmt, &ltime);
         return std::string(out, out + result);
     }

     /** Return the date in ISO-8601 date format.
      *
      *  @return a string containing the date in ISO-8601 date format.
      */
     std::string getISODate()
     {
       return getDate("%F");
   }

     /** Return the date formatted as "Weekday-name, Month-name Day, Year".
      *
      *  @return a string containing the date in the specified format.
      */
     std::string getTextDate()
     {
       return getDate("%F_%X");
     }
 };
//********************************PROTOTYPES****************/
std::string getDateFile();
void writedata2(const char*,double[],double[],double[],int);
std::string convertElementText(int);
double cutter(double number);
std::string cellFill4(std::string,int);
//********************************PROTOTYPES****************/
std::string getDateFile()
   {
     Date d;
     return d.getTextDate();
   }
void writedata2(const char* elements,double theX[],double theY[],
  double theZ[],int sizeElem)
{
  std::string current=getDateFile()+".pdb";
  ofstream file(current);
  int indexFile;
  std::string elementalFile,occu="1.000";
  elementalFile.append(elements);
  std::string currentX,currentY,currentZ;
  std::string chamber1,chamber2,chamber3,chamber4,chamber5,chamber6;
  std::string linesOfFile;
  for(indexFile=0;indexFile<sizeElem;indexFile++)
  {
    linesOfFile="";
    currentX=to_string(cutter(theX[indexFile]));
    currentY=to_string(cutter(theY[indexFile]));
    currentZ=to_string(cutter(theZ[indexFile]));
    chamber1=cellFill4("ATOM",6);
    chamber2=cellFill4(indexFile,4);
    chamber3=cellFill4(elements,3);
    chamber4=cellFill4(currentX,7);
    chamber5=cellFill4(currentY,7);
    chamber6=cellFill4(currentZ,7);
    linesOfFile=chamber1+chamber2+chamber3+chamber4+chamber5+chamber6+occu;
    file <<linesOfFile<<endl;
  }
    file.close();
}
//this function converts the atomic number to the element
std::string convertElementText(int atomicNumber)
{
  std::string elemental="";
  switch (atomicNumber)
  {
    case 1:
    elemental="H";
    break;
    case 2:
    elemental="He";
    break;
    case 3:
    elemental="Li";
    break;
    case 4:
    elemental="Be";
    break;
    case 5:
    elemental="B";
    break;
    case 6:
    elemental="C";
    break;
    case 7:
    elemental="N";
    break;
    case 8:
    elemental="O";
    break;
    case 9:
    elemental="F";
    break;
    case 10:
    elemental="Ne";
    break;
    case 11:
    elemental="Na";
    break;
    case 12:
    elemental="Mg";
    break;
    case 13:
    elemental="Al";
    break;
    case 14:
    elemental="Si";
    break;
    case 15:
    elemental="P";
    break;
    case 16:
    elemental="S";
    break;
    case 17:
    elemental="Cl";
    break;
    case 18:
    elemental="Ar";
    break;
    case 19:
    elemental="K";
    break;
    case 20:
    elemental="Ca";
    break;
    case 21:
    elemental="Sc";
    break;
    case 22:
    elemental="Ti";
    break;
    case 23:
    elemental="V";
    break;
    case 24:
    elemental="Cr";
    break;
    case 25:
    elemental="Mn";
    break;
    case 26:
    elemental="Fe";
    break;
    case 27:
    elemental="Co";
    break;
    case 28:
    elemental="Ni";
    break;
    case 29:
    elemental="Cu";
    break;
    case 30:
    elemental="Zn";
    break;
    case 31:
    elemental="Ga";
    break;
    case 32:
    elemental="Ge";
    break;
    case 33:
    elemental="As";
    break;
    case 34:
    elemental="Se";
    break;
    case 35:
    elemental="Br";
    break;
    case 36:
    elemental="Kr";
    break;
    case 37:
    elemental="Rb";
    break;
    case 38:
    elemental="Sr";
    break;
    case 39:
    elemental="Y";
    break;
    case 40:
    elemental="Zr";
    break;
    case 41:
    elemental="Nb";
    break;
    case 42:
    elemental="Mo";
    break;
    case 43:
    elemental="Tc";
    break;
    case 44:
    elemental="Ru";
    break;
    case 45:
    elemental="Rh";
    break;
    case 46:
    elemental="Pd";
    break;
    case 47:
    elemental="Ag";
    break;
    case 48:
    elemental="Cd";
    break;
    case 49:
    elemental="In";
    break;
    case 50:
    elemental="Sn";
    break;
    case 51:
    elemental="Sb";
    break;
    case 52:
    elemental="Te";
    break;
    case 53:
    elemental="I";
    break;
    case 54:
    elemental="Xe";
    break;
    case 55:
    elemental="Cs";
    break;
    case 56:
    elemental="Ba";
    break;
    case 57:
    elemental="La";
    break;
    case 58:
    elemental="Ce";
    break;
    case 59:
    $elemental="Pr";
    break;
    case 60:
    elemental="Nd";
    break;
    case 61:
    elemental="Pm";
    break;
    case 62:
    elemental="Sm";
    break;
    case 63:
    $elemental="Eu";
    break;
    case 64:
    elemental="Gd";
    break;
    case 65:
    elemental="Tb";
    break;
    case 66:
    elemental="Dy";
    break;
    case 67:
    elemental="Ho";
    break;
    case 68:
    elemental="Er";
    break;
    case 69:
    elemental="Tm";
    break;
    case 70:
    elemental="Yb";
    break;
    case 71:
    elemental="Lu";
    break;
    case 72:
    elemental="Hf";
    break;
    case 73:
    $elemental="Ta";
    break;
    case 74:
    elemental="W";
    break;
    case 75:
    elemental="Re";
    break;
    case 76:
    elemental="Os";
    break;
    case 77:
    elemental="Ir";
    break;
    case 78:
    elemental="Pt";
    break;
    case 79:
    elemental="Au";
    break;
    case 80:
    elemental="Hg";
    break;
    case 81:
    elemental="Tl";
    break;
    case 82:
    elemental="Pb";
    break;
    case 83:
    $elemental="Bi";
    break;
    case 84:
    elemental="Po";
    break;
    case 85:
    elemental="At";
    break;
    case 86:
    elemental="Rn";
    break;
    case 87:
    elemental="Fr";
    break;
    case 88:
    elemental="Ra";
    break;
    case 89:
    elemental="Ac";
    break;
    case 90:
    elemental="Th";
    break;
    case 91:
    elemental="Pa";
    break;
    case 92:
    elemental="U";
    break;
    case 93:
    elemental="Np";
    break;
    case 94:
    elemental="Pu";
    break;
    case 95:
    elemental="Am";
    break;
    case 96:
    elemental="Cm";
    break;
    case 97:
    elemental="Bk";
    break;
    case 98:
    elemental="Cf";
    break;
  }
  return elemental;
}
double cutter(double number)
{
  if(number<0.01)
  {
    return 0;
  }
  else
  {
    return number;
  }
}
std::string cellFill4(std::string statement,int field)
{

  std::string space=" ";
  std::string cellField=statement;
  //fill the data with spaces depending on the field
  int counterOfChar;
  int old=statement.length();
  counterOfChar=old;
  if(old<field)
  {
    while(counterOfChar<old)
    {
      cellField.append(space);
      counterOfChar++;
    }
  }
  return cellField;
}
