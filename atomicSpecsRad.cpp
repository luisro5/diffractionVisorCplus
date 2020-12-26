/**
 *
 *
 * This source code is licensed under the MIT license found in the
 * MIT-LICENSE file in the root directory of this source tree.
 *
 *
 */
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <exception>
#include <fstream>
#include <vector>
#include <algorithm>

struct confi
{
	std::vector<std::string> otherfamily;
};

std::string read_config(std::ifstream& in, confi& out,std::string elementSelected)
{
	in.open("scatteringTables40.txt");
	std::string str;
	while(!in.eof())
	{
		while(getline(in,str))
		{
			std::string::size_type begin = str.find_first_not_of(" \f\t\v");
			//Skips blank lines
			if(begin == std::string::npos)
			continue;
			//Skips #
			if(std::string("#").find(str[begin]) != std::string::npos)
			continue;
			std::string firstWord;
			try
			{
				firstWord = str.substr(0,str.find(" "));
			}
			catch(std::exception& e)
			{
				firstWord = str.erase(str.find_first_of(" "),str.find_first_not_of(" "));
			}
			if(firstWord == elementSelected)
			{
				size_t found = str.find(" ");
				if(found != std::string::npos)
				{
					out.otherfamily.push_back(str.substr(str.find_first_of(" ")+1,found-str.find_first_of(" ")-1));
					out.otherfamily.push_back(str.substr(found+2,str.length()));
				}
			}

		}
	}
	return out.otherfamily[0];
}
int main()
{
	std::ifstream inp;
	confi outp;
	std::string gest="";
	gest=read_config(inp,outp,"As");
	std::cout << gest<<std::endl;
	return 0;
}
