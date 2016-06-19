//
//  PrimitiveGeneratorLoc.cpp
//  PrimitiveGeneratorLoc
//
//  Created by Hunyadi Janos on 31/01/15.
//  Copyright (c) 2015 Janos Hunyadi. All rights reserved.
//

#include "PrimitiveGeneratorLoc.h"

// Plugin

MTypeId PrimitiveGeneratorLoc::id( 0x00123944 );



MString	PrimitiveGeneratorLoc::drawDbClassification("drawdb/geometry/PrimitiveGeneratorLoc");
MString	PrimitiveGeneratorLoc::drawRegistrantId("PrimitiveGeneratorLocPlugin");

PrimitiveGeneratorLoc::PrimitiveGeneratorLoc() {}
PrimitiveGeneratorLoc::~PrimitiveGeneratorLoc() {}

void* PrimitiveGeneratorLoc::creator() { return new PrimitiveGeneratorLoc(); }

PrimitiveGeneratorLocOverride::PrimitiveGeneratorLocOverride(const MObject& obj) : MHWRender::MPxDrawOverride(obj, PrimitiveGeneratorLocOverride::draw){}
PrimitiveGeneratorLocOverride::~PrimitiveGeneratorLocOverride()
{






}

// ----------------------------------------------------------------------------------------------------------------------------------------------------

// VP 1.0 functions

MStatus PrimitiveGeneratorLoc::compute( const MPlug& plug, MDataBlock& data )
{

	return MS::kUnknownParameter;
}


void PrimitiveGeneratorLoc::draw( M3dView & view, const MDagPath & path, M3dView::DisplayStyle style,  M3dView::DisplayStatus status )
{
	MObject thisNode = thisMObject();

	MFnDependencyNode fnDepCloner( thisNode );

	// Draw locator
	view.beginGL();

	glPushAttrib( GL_CURRENT_BIT );

	if ( status == M3dView::kActive ) {
		view.setDrawColor( 13, M3dView::kActiveColors );
	} else {
		view.setDrawColor( 13, M3dView::kDormantColors );
	}



	float r = 1.0f;

	int lats = 10;
	int longs = 10;

	for(int i = 0; i <= lats; i++) {
		double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
		double z0  = sin(lat0);
		double zr0 =  cos(lat0);
		z0 *= r*0.5;
		zr0 *= r*0.5;

		double lat1 = M_PI * (-0.5 + (double) i / lats);
		double z1 = sin(lat1);
		double zr1 = cos(lat1);
		z1 *= r*0.5;
		zr1 *= r*0.5;

		glBegin(GL_QUAD_STRIP);
		for(int j = 0; j <= longs; j++) 
		{
			double lng = 2 * M_PI * (double) (j - 1) / longs;
			double x = cos(lng);
			double y = sin(lng);

			MPoint a(float(x) * float(zr0), float(y) * float(zr0), float(z0));
			MPoint b(float(x) * float(zr1), float(y) * float(zr1), float(z1));


			glVertex3f(float(a.x),float(a.y),float(a.z));
			glVertex3f(float(b.x),float(b.y),float(b.z));
		}
		glEnd();

	}





	// view.drawText( "Switch to VP 2.0", MPoint::origin, M3dView::kCenter );



	glPopAttrib();

	view.endGL();

}


bool PrimitiveGeneratorLoc::isBounded() const
{
	return true;
}

MBoundingBox PrimitiveGeneratorLoc::boundingBox() const
{

	MBoundingBox bbox;

	MFnDependencyNode fnDepLocNode( thisMObject() );

	MPoint corner1( -1.0, 1.0, -1.0 );
	MPoint corner2( 1.0, 0.0, 1.0 );
	return MBoundingBox( corner1, corner2 );


}




// ----------------------------------------------------------------------------------------------------------------------------------------------------

// VP 2.0 Override functions

MHWRender::DrawAPI PrimitiveGeneratorLocOverride::supportedDrawAPIs() const
{

#if MAYA_API_VERSION > 201600

	return (MHWRender::kOpenGL | MHWRender::kDirectX11 | MHWRender::kOpenGLCoreProfile );

#else
	return (MHWRender::kOpenGL | MHWRender::kDirectX11 );
#endif

}

bool PrimitiveGeneratorLocOverride::isBounded(const MDagPath& /*objPath*/, const MDagPath& /*cameraPath*/) const
{
	return true;
}

MBoundingBox PrimitiveGeneratorLocOverride::boundingBox( const MDagPath& objPath, const MDagPath& cameraPath) const
{


	MStatus status;
	MObject CLonerMultiNode = objPath.node(&status);


	MPoint corner1( -1.0, 1.0, -1.0 );
	MPoint corner2( 1.0, 0.0, 1.0 );

	return MBoundingBox( corner1, corner2 );

	//return m_bbP;
}


// Called by Maya each time the object needs to be drawn.
MUserData* PrimitiveGeneratorLocOverride::prepareForDraw( const MDagPath& objPath, const MDagPath& cameraPath, const MHWRender::MFrameContext& frameContext, MUserData* oldData)
{

	// Get outside data from plugs
	MStatus status;
	MObject CLonerMultiNode = objPath.node(&status);


	//

	// Add data
	PrimitiveGeneratorLocData* data = dynamic_cast<PrimitiveGeneratorLocData*>(oldData);
	if (!data)
	{
		data = new PrimitiveGeneratorLocData();
	}




	data->m_inLoc_mat = objPath.exclusiveMatrix();

	// get correct color based on the state of object, e.g. active or dormant
	data->m_locColor = MHWRender::MGeometryUtilities::wireframeColor(objPath);

	return data;
}


void PrimitiveGeneratorLocOverride::addUIDrawables( const MDagPath& objPath, MHWRender::MUIDrawManager& drawManager, const MHWRender::MFrameContext& frameContext, const MUserData* data)
{
	PrimitiveGeneratorLocData* pLocatorData = (PrimitiveGeneratorLocData*)data;
	if (!pLocatorData)
	{
		return;
	}

	drawManager.beginDrawable();



	M3dView view = M3dView::active3dView();
	short ox, oy;


	MColor fillCol, lineCol;

	fillCol = MColor( 1.0f, 0.0f, 0.0f, 0.5f );
	lineCol = MColor( 1.0f, 1.0f, 1.0f, 1.0f );
	


	if ( frameContext.getDisplayStyle() & MHWRender::MDrawContext::kWireFrame ) {


		fillCol = MColor( pLocatorData->m_locColor.r, pLocatorData->m_locColor.g, pLocatorData->m_locColor.b, 0.0 );

	}

	if ( MHWRender::MGeometryUtilities::displayStatus(objPath) == M3dView::kLead ) {


		fillCol = MColor( pLocatorData->m_locColor.r, pLocatorData->m_locColor.g, pLocatorData->m_locColor.b, fillCol.a );
		
	}

	if ( MHWRender::MGeometryUtilities::displayStatus(objPath) == M3dView::kActive ) {

		fillCol = MColor( pLocatorData->m_locColor.r, pLocatorData->m_locColor.g, pLocatorData->m_locColor.b, fillCol.a );

	}

	if ( MHWRender::MGeometryUtilities::displayStatus(objPath) == M3dView::kTemplate ) {

		fillCol = MColor( pLocatorData->m_locColor.r, pLocatorData->m_locColor.g, pLocatorData->m_locColor.b, fillCol.a );
		
	}


	


	MPoint p = MPoint::origin;
	p *= pLocatorData->m_inLoc_mat;

	view.worldToView(p, ox, oy);

	drawManager.setColor( MColor(fillCol.r,fillCol.g,fillCol.b, 0.8f) );
	drawManager.circle2d(MPoint(ox,oy),4.0, true );

	drawManager.setColor( MColor(lineCol.r,lineCol.g,lineCol.b, 0.8f) );
	drawManager.circle2d(MPoint(ox,oy),4.0, false );





	drawManager.endDrawable();
}

MStatus PrimitiveGeneratorLoc::initialize()
{
	MFnMatrixAttribute		mAttr;
	MFnNumericAttribute		nAttr;
	MStatus					status;



	return MS::kSuccess;
}
