#ifndef PRIMITIVEGENDCOMMAND_H
#define PRIMITIVEGENDCOMMAND_H

#ifdef __linux__
#include <maya/MArgDatabase.h>
#else
#include <maya/MArgDataBase.h>
#endif

#include <maya/MDagPath.h>
#include <maya/MDGModifier.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMesh.h>
#include <maya/MGlobal.h>
#include <maya/MIntArray.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MMeshIntersector.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MPointArray.h>
#include <maya/MPxCommand.h>
#include <maya/MSelectionList.h>
#include <maya/MSyntax.h>
#include <maya/MDagModifier.h>
#include <maya/MFnSet.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnTransform.h>
#include <maya/MRampAttribute.h>

class primitiveGeneratorCommand : public MPxCommand
{
public:
    primitiveGeneratorCommand();
    virtual MStatus doIt( const MArgList& argList );
    virtual MStatus redoIt();
    virtual MStatus undoIt();
    virtual bool isUndoable() const;
    static void* creator();
    static MSyntax newSyntax();

	

private:

	MStatus createPrimGenFromCurves(MDagPathArray p_currSelTrA, MDagPathArray p_currSelShapeA);
	MStatus createPrimGenFromLocators();

	MDagModifier m_DAGMod;
	MDGModifier m_DGMod;
	MFnDependencyNode m_DEPNode;

	MObjectArray o_outputMeshA;
	MObjectArray o_primGenNodeA;
	MObject		 o_locA;
	MObject		 o_locB;

	bool		m_muscle;

	

};


#endif