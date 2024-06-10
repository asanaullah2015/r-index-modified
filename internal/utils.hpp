// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#ifndef UTILS_RI_HPP_
#define UTILS_RI_HPP_

#include <sstream>

//using namespace std;

using ulint = uint64_t;

std::string get_time(uint64_t time){

	std::stringstream ss;

	if(time>=3600){

		uint64_t h = time/3600;
		uint64_t m = (time%3600)/60;
		uint64_t s = (time%3600)%60;

		ss  << time << " seconds. ("<< h << "h " << m << "m " << s << "s" << ")";

	}else if (time>=60){

		uint64_t m = time/60;
		uint64_t s = time%60;

		ss << time << " seconds. ("<< m << "m " << s << "s" << ")";

	}else{

		ss << time << " seconds.";

	}

	return ss.str();

}

uint8_t bitsize(uint64_t x){

	if(x==0) return 1;
	return 64 - __builtin_clzll(x);

}

//parse pizza&chilli patterns header:
void header_error(){
	std::cout << "Error: malformed header in patterns file" << std::endl;
	std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
	exit(0);
}

ulint get_number_of_patterns(std::string header){

	ulint start_pos = header.find("number=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

ulint get_patterns_length(std::string header){

	ulint start_pos = header.find("length=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}



#endif /* UTILS_RI_HPP_ */
