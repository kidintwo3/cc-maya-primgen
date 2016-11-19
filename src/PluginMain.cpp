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

#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
{
	MStatus		status;
	MFnPlugin	fnPlugin( obj, "Creative Case", "1.6", "Any");

	icons_data_write();

	MGlobal::executeCommand( mel_createShelf() );
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

