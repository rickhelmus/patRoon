/*
 *  main.cpp, part of GenForm by M. Meringer Copyright (C) 2015
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
*/

#include <iostream>
#include <map>

using namespace std;

extern int GenFormMatchIsotopeMsMs(const string& strProgName,map<string,string>& mapArgValue);

void ParseCommandLine(int argc, char *argv[],map<string,string>& mapParamValue)
{
	while(--argc>0)
	{
		string strArg=*++argv;
		unsigned int iPos=strArg.find('=');
		mapParamValue[strArg.substr(0,iPos)]=iPos==string::npos?"":strArg.substr(iPos+1);
	}
}

int main(int argc, char *argv[]) {
	map<string,string> mapParamValue;
	ParseCommandLine(argc,argv,mapParamValue);
	return GenFormMatchIsotopeMsMs(*argv,mapParamValue);
}
