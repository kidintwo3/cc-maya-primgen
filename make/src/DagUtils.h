//
//  connectToNode.h
//  hairDrawTool
//
//  Created by Hunyadi János on 2015. 01. 03..
//  Copyright (c) 2015. Janos Hunyadi. All rights reserved.
//

#ifndef DAGUTILS_H
#define DAGUTILS_H

//#include <stdio.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MItDag.h>
#include <maya/MString.h>
#include <maya/MGlobal.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MPlug.h>
#include <maya/MFnDagNode.h>
#include <maya/MDGModifier.h>
#include <maya/MDagModifier.h>
#include <maya/MItDependencyNodes.h>
#include <maya/MFnMesh.h>
#include <maya/MFnSet.h>

//


// Create Dag
MObject             createNodeCustom(MDGModifier& m_DAGMod, MString id);
MObject             createNodeMaya(MFnDependencyNode& m_DEPNode, MString id);

// Check Dag
MSelectionList      listNodeType(MFn::Type mfn_nodeType);
bool                checkNodeExist(MSelectionList selList, MString searchString);

//Paths
MStatus				getShapeNodeFromTransformDAG(MDagPath& path);
MDagPath            getDagFromString(MString nodeName);
MDagPath            getCurrSelectionDAG();
MDagPathArray		getCurrSelectionDAGArray();

// Plugs
MPlug               findPlug(MDagPath pathName, MString plugName);
MStatus             setPlugs(MObject obj, MString plug, MString value);
MStatus             connectPlug(MDagModifier& m_DAGMod, MObject objSource, MObject objTarget, MString plugSourceStr, MString plugTargetStr);
MStatus             deleteNode(MObject obj);


// Material
bool                checkMatExist(MString matName);
MStatus				assignSameMaterial(MDagPath& inputShapeDagPath, MObject& outputShapeDagPath);
MStatus				assignInitialShadingGroup(MObject& oMesh);

#endif
