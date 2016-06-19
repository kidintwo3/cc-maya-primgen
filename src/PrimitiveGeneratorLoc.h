//
//  PrimitiveGeneratorLoc.h
//  PrimitiveGeneratorLoc
//
//  Created by Hunyadi Janos on 31/01/15.
//  Copyright (c) 2015 Janos Hunyadi. All rights reserved.
//

#ifndef PrimitiveGeneratorLoc_H
#define PrimitiveGeneratorLoc_H

#include <maya/MPxLocatorNode.h>
#include <maya/M3dView.h>
#include <maya/MGlobal.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MSelectionList.h>
#include <maya/MFnCamera.h>

#include <maya/MColorArray.h>
#include <maya/MTime.h>
#include <maya/MAnimControl.h>

#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnGenericAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnStringData.h>
#include <maya/MStringArray.h>
#include <maya/MEulerRotation.h>

#include <maya/MFnDependencyNode.h>

#include <maya/MFileIO.h>

#include <maya/MPlugArray.h>
#include <maya/MPointArray.h>

// Viewport 2.0 includes
#include <maya/MDrawRegistry.h>
#include <maya/MPxDrawOverride.h>
#include <maya/MUserData.h>
#include <maya/MDrawContext.h>
#include <maya/MHWGeometryUtilities.h>

//

#include <vector>
#include <sstream>
#include <ctime>

using namespace std;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Classic 1.0 viewport implementation
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class PrimitiveGeneratorLoc : public MPxLocatorNode
{
public:
	PrimitiveGeneratorLoc();
	virtual ~PrimitiveGeneratorLoc();

    virtual MStatus   		compute( const MPlug& plug, MDataBlock& data );

	virtual void            draw( M3dView & view, const MDagPath & path, M3dView::DisplayStyle style, M3dView::DisplayStatus status );

	virtual bool            isBounded() const;
	virtual MBoundingBox    boundingBox() const;

	static  void *          creator();
	static  MStatus         initialize();

public:
	static	MTypeId			id;

	// Vp 2.0 stuff
	static	MString			drawDbClassification;
	static	MString			drawRegistrantId;

private:
	
	MBoundingBox			m_bbP;
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Viewport 2.0 override implementation
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


// data for override

class PrimitiveGeneratorLocData : public MUserData
{
public:
	PrimitiveGeneratorLocData() : MUserData(false) {} // don't delete after draw
	virtual ~PrimitiveGeneratorLocData() {}


	MMatrix					m_inLoc_mat;

	MColor					m_locColor;
	MPoint					m_controllerPos;

};


// override

class PrimitiveGeneratorLocOverride : public MHWRender::MPxDrawOverride
{
public:
	static MHWRender::MPxDrawOverride* Creator(const MObject& obj)
	{
		return new PrimitiveGeneratorLocOverride(obj);
	}

	virtual ~PrimitiveGeneratorLocOverride();
	virtual MHWRender::DrawAPI supportedDrawAPIs() const;

	virtual bool isBounded( const MDagPath& objPath, const MDagPath& cameraPath) const;
	virtual MBoundingBox boundingBox( const MDagPath& objPath, const MDagPath& cameraPath) const;

	virtual MUserData* prepareForDraw( const MDagPath& objPath, const MDagPath& cameraPath, const MHWRender::MFrameContext& frameContext, MUserData* oldData);

	virtual bool hasUIDrawables() const { return true; }
	virtual void addUIDrawables( const MDagPath& objPath, MHWRender::MUIDrawManager& drawManager, const MHWRender::MFrameContext& frameContext, const MUserData* data);

	static void draw(const MHWRender::MDrawContext& context, const MUserData* data) {};

private:
	PrimitiveGeneratorLocOverride(const MObject& obj);

};


#endif
