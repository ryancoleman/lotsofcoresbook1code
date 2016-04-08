//	$Id: file_exists.cpp 19708 2010-10-29 18:04:21Z d3y133 $

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
  ifstream is(argv[1]);
  return bool(is);
}
