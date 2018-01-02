//-
// ==========================================================================
// Copyright 2015 Autodesk, Inc.  All rights reserved.
//
// Use of this software is subject to the terms of the Autodesk
// license agreement provided at the time of installation or download,
// or which otherwise accompanies this software in either electronic
// or hard copy form.
// ==========================================================================
//+

/*

This is an example to demonstrate the use of a rotation manipulator through
a rotation tool and context.  This example uses three classes to accomplish
this task: First, a context command (rotateContext) is provided to create
instances of the context.  Next, a custom selection context
(RotateManipContext) is created to manage the rotation manipulator.
Finally, the rotation manipulator is provided as a custom node class.

Loading and unloading:
----------------------

The rotate manipulator context can be created with the
following mel commands:

rotateContext;
setToolTo rotateContext1;

If the preceding commands were used to create the manipulator context,
the following commands can destroy it:

deleteUI rotateContext1;
deleteUI rotateManip;

If the plugin is loaded and unloaded frequently (eg. during testing),
it is useful to make these command sequences into shelf buttons.

How to use:
-----------

Once the tool button has been created using the script above, select the
tool button then click on an object.  The rotate manipulator should appear
at the center of the selected object.  The rotate manipulator can be used
much like the built-in rotate manipulator.  In addition, the plugin
produces a state manipulator that can be used to control the modes for the
rotation manipulator.  The state manipulator should be displayed 2 units
along the X-axis from the object.

*/

#include <maya/MIOStream.h>
#include <stdio.h>
#include <stdlib.h>

#include <maya/MFn.h>
#include <maya/MPxNode.h>
#include <maya/MPxManipContainer.h>
#include <maya/MPxSelectionContext.h>
#include <maya/MPxContextCommand.h>
#include <maya/MModelMessage.h>
#include <maya/MGlobal.h>
#include <maya/MItSelectionList.h>
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <maya/MDagPath.h>
#include <maya/MManipData.h>
#include <maya/MEulerRotation.h>

// Manipulators
#include <maya/MFnRotateManip.h>
#include <maya/MFnStateManip.h>


/////////////////////////////////////////////////////////////
//
// rotateContext
//
// This is the command that will be used to create instances
// of our context.
//
/////////////////////////////////////////////////////////////

class rotateContext : public MPxContextCommand
{
public:
	rotateContext() {};
	virtual MPxContext * makeObj();

public:
	static void* creator();
};


class RotateManipContext : public MPxSelectionContext
{
public:
	RotateManipContext();
	virtual void	toolOnSetup(MEvent &event);
	virtual void	toolOffCleanup();

	// Callback issued when selection list changes
	static void updateManipulators(void * data);

private:
	MCallbackId id1;
};




/////////////////////////////////////////////////////////////
//
// exampleRotateManip
//
// This class implements the example rotate manipulator.
//
/////////////////////////////////////////////////////////////

class exampleRotateManip : public MPxManipContainer
{
public:
	exampleRotateManip();
	virtual ~exampleRotateManip();

	static void * creator();
	static MStatus initialize();
	virtual MStatus createChildren();
	virtual MStatus connectToDependNode(const MObject &node);

	virtual void draw(M3dView &view,
		const MDagPath &path,
		M3dView::DisplayStyle style,
		M3dView::DisplayStatus status);

	// Callback function
	MManipData rotationChangedCallback(unsigned index);

public:
	static MTypeId id;

private:
	MDagPath fRotateManip;
	MDagPath fStateManip;

	unsigned rotatePlugIndex;
};
