//
//  PluginMain.cpp
//  hairDrawTool
//
//  Created by Hunyadi János on 2014. 11. 05..
//  Copyright (c) 2014. Janos Hunyadi. All rights reserved.
//


#include "PrimitiveGeneratorNode.h"
#include "PrimitiveGeneratorLoc.h"
#include "PrimitiveGeneratorCommand.h"
#include "AETemplate.h"
#include "icons.h"
//#include "LicCheck.h"


#include <maya/MFnPlugin.h>
#include <maya/MCommonSystemUtils.h>

MStatus initializePlugin( MObject obj )
{
	MStatus		status;
	MFnPlugin	fnPlugin( obj, "Creative Case", "2.0 beta for Capcom", "Any");


	//bool check_lic = false;

	//// Licence Check

	//if (check_lic)
	//{



	//	MString plugin_name = "PrimGen";

	//	// MString lic_key = MString(std::getenv("PRIMGEN_LICENSE_KEY"));
	//	MString lic_key = "D0F9F4D8-50FD4FE3-A411C173-FB6C77C7";

	//	if (!checkLic(lic_key))
	//	{
	//		MGlobal::displayInfo( "--------------------------------------------");
	//		MGlobal::displayInfo("["+plugin_name+"] No licence found!");
	//		MGlobal::displayInfo( "Make sure the Maya.env file has your licence");
	//		MGlobal::displayInfo( "--------------------------------------------");

	//		MGlobal::executeCommand("inViewMessage -smg \"["+plugin_name+"] No licence found!\" -pos midCenter -bkc 0x00750000 -fade;");

	//		return MStatus::kFailure;

	//	} 

	//	else
	//	{
	//		MGlobal::displayInfo("["+plugin_name+"] Licence found!");
	//	}


	//}

	// Icon / UI rebuild check

	MString rebuild_icons = MCommonSystemUtils::getEnv("PRIMGEN_REBUILD_ICONS", &status);

	if( !rebuild_icons.asShort() )
	{
		icons_data_write();
	}

	MString rebuild_shelf = MCommonSystemUtils::getEnv("PRIMGEN_REBUILD_SHELF", &status);

	if( !rebuild_shelf.asShort() )
	{
		MGlobal::executeCommand( mel_createShelf() );
	}



	MGlobal::executeCommand( mel_AETemplate() );

	status = fnPlugin.registerCommand( "primitiveGeneratorCommand", primitiveGeneratorCommand::creator, primitiveGeneratorCommand::newSyntax );
	CHECK_MSTATUS_AND_RETURN_IT( status );

	status = fnPlugin.registerNode( "primitiveGenerator", primitiveGenerator::id, primitiveGenerator::creator, primitiveGenerator::initialize );
	CHECK_MSTATUS_AND_RETURN_IT( status );

	// Locator
	status = fnPlugin.registerNode( "PrimitiveGeneratorLoc", PrimitiveGeneratorLoc::id, &PrimitiveGeneratorLoc::creator, &PrimitiveGeneratorLoc::initialize, MPxNode::kLocatorNode, &PrimitiveGeneratorLoc::drawDbClassification);
	CHECK_MSTATUS_AND_RETURN_IT(status);
	status = MHWRender::MDrawRegistry::registerDrawOverrideCreator( PrimitiveGeneratorLoc::drawDbClassification, PrimitiveGeneratorLoc::drawRegistrantId, PrimitiveGeneratorLocOverride::Creator);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	return status;
}

MStatus uninitializePlugin( MObject obj )
{
	MStatus		status;
	MFnPlugin	fnPlugin( obj );

	status = fnPlugin.deregisterCommand( "primitiveGeneratorCommand" );
	CHECK_MSTATUS_AND_RETURN_IT( status );

	status = fnPlugin.deregisterNode( primitiveGenerator::id );
	CHECK_MSTATUS_AND_RETURN_IT( status );

	// Locator 
	status = MHWRender::MDrawRegistry::deregisterDrawOverrideCreator( PrimitiveGeneratorLoc::drawDbClassification, PrimitiveGeneratorLoc::drawRegistrantId);
	CHECK_MSTATUS_AND_RETURN_IT(status);
	status = fnPlugin.deregisterNode( PrimitiveGeneratorLoc::id );
	CHECK_MSTATUS_AND_RETURN_IT(status);


	return status;
}

