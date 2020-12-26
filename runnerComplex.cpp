/**
 *
 *
 * This source code is licensed under the MIT license found in the
 * MIT-LICENSE file in the root directory of this source tree.
 *
 *
 */
#include "complexLibrary.h"
//runnerComplex
//g++ -o test10 runnerComplex.cpp -lm
int main()
{
   /********************************************/
   //test of fft
    const Complex test[] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
    CArray data(test, 8);

    // forward fft
    fft(data);

    std::cout << "fft" << std::endl;
    for (int i = 0; i < 8; ++i)
    {
        std::cout << data[i] << std::endl;
    }
    /********************************************/
    // inverse fft
    ifft(data);

    std::cout << std::endl << "ifft" << std::endl;
    for (int i = 0; i < 8; ++i)
    {
        std::cout << data[i] << std::endl;
    }
		/********************************************/
    //test of  transposeShift
		int test2[] = { 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0};
		int large=16,l;
		int baz[large] = {};
		initial(baz,16);
		transposeShift(test2,baz,4);
    std::cout << "Transpose shift" << std::endl;
		for (l = 0; l < large; l++)
    {
        std::cout << baz[l] << std::endl;
    }
    /********************************************/
    //test of  fftshift2d
    int test3[] = { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
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
		int test4[large] = {};
    x0Filter(test4,4,4,1);
    std::cout << "x0Filter" << std::endl;
    for (l = 0; l < large; l++)
    {
        std::cout << test4[l] << std::endl;
    }
    return 0;
}
