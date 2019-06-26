#ifndef ICONS_H
#define ICONS_H

#include "icons_repo.h"

void writeIcon_binary(MString path, MString filename, const unsigned char output[], int char_size)
{

	ofstream myfile;

	myfile.open(MString(path + filename).asChar(), ios::out | ios::binary);
	myfile.write((char*)output, char_size);
	myfile.close();
	MGlobal::displayInfo(MString() + "[ShellMod] Created icon: " + path + filename);


}



void icons_data_write()
{


	MString path;
	MGlobal::executeCommand("internalVar -userBitmapsDir", path);
	MGlobal::displayInfo(MString() + "[ShellMod] Icons path: " + path);

	std::string c_path = path.asChar();

	writeIcon_binary(path, "primitiveGenerator_muscle.png", _acprimitiveGenerator_muscle, 16432UL + 1);
	writeIcon_binary(path, "primitiveGenerator_CCLogo.png", _acprimitiveGenerator_CCLogo, 19052UL + 1);
	writeIcon_binary(path, "primitiveGenerator.png", _acprimitiveGenerator, 1567UL + 1);
	writeIcon_binary(path, "out_PrimitiveGeneratorLoc.png", _acout_PrimitiveGeneratorLoc, 21932UL + 1);

}

#endif