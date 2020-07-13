//
//  PluginMain.cpp
//  hairDrawTool
//
//  Created by Hunyadi János on 2014. 11. 05..
//  Copyright (c) 2014. Janos Hunyadi. All rights reserved.
//

//#include <maya/MCppCompat.h>
#include "PrimitiveGeneratorNode.h"
#include "PrimitiveGeneratorLoc.h"
#include "PrimitiveGeneratorCommand.h"
#include "PrimitiveGeneratorManip.h"
#include "AETemplate.h"
#include "icons.h"
//#include "LicCheck.h"


#include <maya/MFnPlugin.h>
#include <maya/MCommonSystemUtils.h>

MStatus initializePlugin( MObject obj )
{
	MStatus		status;
	MFnPlugin	fnPlugin( obj, "Creative Case", "2.12", "Any");


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

	MStringArray aeTemplateA = mel_AETemplate();

	for (int i = 0; i < aeTemplateA.length(); i++)
	{
		MGlobal::executeCommand(aeTemplateA[i]);
	}


	status = fnPlugin.registerCommand( "primitiveGeneratorCommand", primitiveGeneratorCommand::creator, primitiveGeneratorCommand::newSyntax );
	CHECK_MSTATUS_AND_RETURN_IT( status );

	status = fnPlugin.registerNode( "primitiveGenerator", primitiveGenerator::id, primitiveGenerator::creator, primitiveGenerator::initialize );
	CHECK_MSTATUS_AND_RETURN_IT( status );

	// Locator
	status = fnPlugin.registerNode( "PrimitiveGeneratorLoc", PrimitiveGeneratorLoc::id, &PrimitiveGeneratorLoc::creator, &PrimitiveGeneratorLoc::initialize, MPxNode::kLocatorNode, &PrimitiveGeneratorLoc::drawDbClassification);
	CHECK_MSTATUS_AND_RETURN_IT(status);
	status = MHWRender::MDrawRegistry::registerDrawOverrideCreator( PrimitiveGeneratorLoc::drawDbClassification, PrimitiveGeneratorLoc::drawRegistrantId, PrimitiveGeneratorLocOverride::Creator);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// Manipulator

	//status = fnPlugin.registerContextCommand("rotateContext",&rotateContext::creator);
	//if (!status) {
	//	MGlobal::displayError("Error registering rotateContext command");
	//	return status;
	//}

	//status = fnPlugin.registerNode("exampleRotateManip", exampleRotateManip::id,
	//	&exampleRotateManip::creator, &exampleRotateManip::initialize,
	//	MPxNode::kManipContainer);
	//if (!status) {
	//	MGlobal::displayError("Error registering rotateManip node");
	//	return status;
	//}

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

	// Manipulator

	//status = fnPlugin.deregisterContextCommand("rotateContext");
	//if (!status) {
	//	MGlobal::displayError("Error deregistering rotateContext command");
	//	return status;
	//}

	//status = fnPlugin.deregisterNode(exampleRotateManip::id);
	//if (!status) {
	//	MGlobal::displayError("Error deregistering RotateManip node");
	//	return status;
	//}

	return status;
}



