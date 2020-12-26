/**
 *
 *
 * This source code is licensed under the MIT license found in the
 * MIT-LICENSE file in the root directory of this source tree.
 *
 *
 */
#include "latticeLibrary.h"
//runnerLattice
//g++ -o test11 runnerLattice.cpp -lm
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
int main()
{
    /********************************************/
    //test of rotate matrix
    int test5[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};


		int baz2[large] = {};
		initial(baz2,16);
    fftshift2d(test3,baz2,4);
    std::cout << "fft 2d shift" << std::endl;
    for (l = 0; l < large; l++)
    {
        std::cout << baz2[l] << std::endl;
    }
    /********************************************/
    //test of  x0Filter
    /*
		int test6[large] = {};
    x0Filter(test4,4,4,1);
    std::cout << "x0Filter" << std::endl;
    for (l = 0; l < large; l++)
    {
        std::cout << test4[l] << std::endl;
    }
    */
    return 0;
}
