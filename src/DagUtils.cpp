//
//  connectToNode.cpp
//  hairDrawTool
//
//  Created by Hunyadi János on 2015. 01. 03..
//  Copyright (c) 2015. Janos Hunyadi. All rights reserved.
//

#include "DagUtils.h"

MSelectionList listNodeType(MFn::Type mfn_nodeType)
{
	MSelectionList list;
	MItDag iter( MItDag::kDepthFirst, mfn_nodeType);

	for ( ; !iter.isDone(); iter.next() )
	{
		MDagPath dagPath;
		iter.getPath(dagPath);

		dagPath.pop(0);
		list.add(dagPath);



		iter.next();
	}

	return list;
}

bool checkNodeExist(MSelectionList selList, MString searchString)
{
	MItSelectionList iter( selList );
	MStringArray stringArray;

	for ( ; !iter.isDone(); iter.next() )
	{

		iter.getStrings(stringArray);

		if (stringArray[0]  == searchString)
		{
			return true;
			break;
		}
		iter.next();

	}

	return false;
}

MDagPath getDagFromString(MString nodeName)
{
	MSelectionList selList;
	selList.add(nodeName);
	MDagPath mDagPath;
	selList.getDagPath(0, mDagPath);

	return mDagPath;
}

MDagPath getCurrSelectionDAG()
{
	MSelectionList selectedObjects;
	MDagPath mDagPath;
	MGlobal::getActiveSelectionList(selectedObjects);

	for (unsigned int i = 0; i < selectedObjects.length(); i++)
	{
		selectedObjects.getDagPath(i, mDagPath);

	/*	if (mDagPath.apiType() == currType)
		{*/
		return mDagPath;
		/*}*/
	}

	return mDagPath;
}

MDagPathArray getCurrSelectionDAGArray()
{
	MSelectionList selectedObjects;
	MDagPath mDagPath;
	MDagPathArray retDagArray;
	MGlobal::getActiveSelectionList(selectedObjects);

	for (unsigned int i = 0; i < selectedObjects.length(); i++)
	{
		selectedObjects.getDagPath(i, mDagPath);
		retDagArray.append(mDagPath);


	}

	return retDagArray;
}

MPlug findPlug(MDagPath pathName, MString plugName)
{
	MDagPath m_currShape = pathName;
	m_currShape.extendToShape();

	MFnDagNode fn_DagMod( m_currShape );
	MPlug plug_currPlugArray(fn_DagMod.findPlug( plugName ));

	MPlug plug_currPlug;

	if (plug_currPlugArray.isArray() == true){
		int count = plug_currPlugArray.numConnectedElements();
		plug_currPlug = plug_currPlugArray.elementByLogicalIndex(count);
	}

	if (plug_currPlugArray.isArray() == false){
		plug_currPlug = plug_currPlugArray;
	}

	return plug_currPlug;
}

MObject createNodeMaya(MFnDependencyNode& m_DEPNode, MString id)
{

	MObject obj = m_DEPNode.create(id);

	return  obj;
}


MObject createNodeCustom(MDGModifier& m_DAGMod, MString id)
{

	//	MDGModifier m_DAGMod;
	MObject obj = m_DAGMod.createNode(id);
	m_DAGMod.doIt();

	return  obj;
}

MStatus setPlugs(MObject obj, MString plug, MString value){
	MStatus status;

	MPlug plugTarget;
	MDagPath pathTarget;

	if (obj.apiType() == MFn::kPluginDependNode) {
		MFnDependencyNode fnDepSource( obj );
		plugTarget = fnDepSource.findPlug( plug );
	}

	if (obj.apiType() == MFn::kTransform) {

		MFnDagNode fnDagT(obj);
		fnDagT.getPath(pathTarget);
		pathTarget.extendToShape();

		plugTarget = findPlug(pathTarget, plug);

	}

	if (obj.apiType() == MFn::kMesh || obj.apiType() == MFn::kNurbsCurve) {
		MFnDagNode fnDagT(obj);
		fnDagT.getPath(pathTarget);

		plugTarget = findPlug(pathTarget, plug);


	}

	if (value.isInt() && value != "true" && value != "false" ) {
		plugTarget.setInt(value.asInt());
	}

	if (value.isDouble() && value != "true" && value != "false" ) {
		plugTarget.setDouble(value.asDouble());
	}

	if (value == "true" && value != "1") {
		plugTarget.setBool(true);
	}

	if (value == "false" && value != "0") {
		plugTarget.setBool(false);
	}


	return status;
}

MStatus connectPlug(MDagModifier& m_DAGMod, MObject objSource, MObject objTarget, MString plugSourceStr, MString plugTargetStr){

	MStatus status;

	//MDagModifier m_dagMod;

	MPlug plugSource, plugTarget;
	MDagPath pathSource, pathTarget;

	// Source

	if (objSource.apiType() == MFn::kPluginDependNode) {
		MFnDependencyNode fnDepSource( objSource );
		plugSource = fnDepSource.findPlug( plugSourceStr );
	}

	if (objSource.apiType() == MFn::kTransform) {

		MFnDagNode fnDagS(objSource);
		fnDagS.getPath(pathSource);
		pathSource.extendToShape();

		plugSource = findPlug(pathSource, plugSourceStr);
	}

	if (objSource.apiType() == MFn::kMesh || objSource.apiType() == MFn::kNurbsCurve) {

		MFnDagNode fnDagS(objSource);
		fnDagS.getPath(pathSource);

		plugSource = findPlug(pathSource, plugSourceStr);
	}



	// Target

	if (objTarget.apiType() == MFn::kPluginDependNode) {
		MFnDependencyNode fnDepTarget( objTarget );
		plugTarget = fnDepTarget.findPlug( plugTargetStr );
	}

	if (objTarget.apiType() == MFn::kTransform) {

		MFnDagNode fnDagT(objTarget);
		fnDagT.getPath(pathTarget);
		pathTarget.extendToShape();

		plugTarget = findPlug(pathTarget, plugTargetStr);
	}

	if (objTarget.apiType() == MFn::kMesh || objTarget.apiType() == MFn::kNurbsCurve) {

		MFnDagNode fnDagT(objTarget);
		fnDagT.getPath(pathTarget);

		plugTarget = findPlug(pathTarget, plugTargetStr);
	}

	// Connect


	m_DAGMod.connect(plugSource, plugTarget);
	m_DAGMod.doIt();

	return status;
}


MStatus deleteNode(MObject obj)
{

	MStatus status;

	MDagPath pathTarget;
	MDagPath pathTargetShape;

	pathTarget = MDagPath::getAPathTo(obj);
	pathTargetShape = pathTarget;
	pathTargetShape.extendToShape();


	if ( pathTarget.isValid() )
	{
		char buffer[512];
		sprintf( buffer, "delete %s", pathTarget.partialPathName().asChar() );
		MGlobal::executeCommand( buffer );

	}

	if ( pathTargetShape.isValid() )
	{
		char buffer[512];
		sprintf( buffer, "delete %s", pathTargetShape.partialPathName().asChar() );
		MGlobal::executeCommand( buffer );

	}



	return status;
}

MStatus getShapeNodeFromTransformDAG(MDagPath& path)
{
	MStatus status;

	if (path.apiType() == MFn::kMesh)
	{
		return MS::kSuccess;
	}

	unsigned int numShapes;
	status = path.numberOfShapesDirectlyBelow(numShapes);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	for (unsigned int i = 0; i < numShapes; ++i)
	{
		status = path.extendToShapeDirectlyBelow(i);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		if (!path.hasFn(MFn::kMesh))
		{
			path.pop();
			continue;
		}

		MFnDagNode fnNode(path, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);
		if (!fnNode.isIntermediateObject())
		{
			return MS::kSuccess;
		}
		path.pop();
	}

	MGlobal::displayWarning(MString() + "Selection is not a mesh");

	return MS::kFailure;

}

bool checkMatExist(MString matName){

	// search for material

	MItDependencyNodes mat_ItDag(MFn::kShadingEngine);

	while ( !mat_ItDag.isDone() )
	{

		MObject o_currMat = mat_ItDag.thisNode();
		if (MFnDependencyNode(o_currMat).name() == matName) {

			return true;

			break;
		}


		mat_ItDag.next();
	}

	return false;
}

MStatus assignSameMaterial(MDagPath& inputShapeDagPath, MObject& outputShapeDagPath)
{

	MStatus status;

	MString sMaterial;

	if (inputShapeDagPath.hasFn(MFn::kMesh))
	{
		// Find the Shading Engines Connected to the SourceNode 
		MFnMesh fnMesh(inputShapeDagPath.node());

		// A ShadingGroup will have a MFnSet 
		MObjectArray sets, comps;
		fnMesh.getConnectedSetsAndMembers(inputShapeDagPath.instanceNumber(), sets, comps, true);

		// Each set is a Shading Group. Loop through them
		for (unsigned int i = 0; i < sets.length(); ++i)
		{
			MFnDependencyNode fnDepSGNode(sets[i]);

			sMaterial = fnDepSGNode.name();
		}
	}


	MGlobal::displayInfo(MString() + "Initial SG: " + sMaterial);

	MSelectionList sList;
	MGlobal::getSelectionListByName(sMaterial, sList);
	MObject oInitialShadingGroup;
	status = sList.getDependNode(0, oInitialShadingGroup);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	MFnSet fnShadingGroup(oInitialShadingGroup, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	status = fnShadingGroup.addMember(outputShapeDagPath);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	return MS::kSuccess;
}

MStatus assignInitialShadingGroup(MObject& oMesh)
{
	MStatus status;

	MSelectionList sList;
	MGlobal::getSelectionListByName("initialShadingGroup", sList);
	MObject oInitialShadingGroup;
	status = sList.getDependNode(0, oInitialShadingGroup);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	MFnSet fnShadingGroup(oInitialShadingGroup, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	status = fnShadingGroup.addMember(oMesh);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	return MS::kSuccess;
}