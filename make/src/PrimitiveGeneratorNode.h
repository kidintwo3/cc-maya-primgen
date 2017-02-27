//
//  PrimitiveGenerator.h
//  hairDrawTool
//
//  Created by Hunyadi János on 2014. 12. 21..
//  Copyright (c) 2014. Janos Hunyadi. All rights reserved.
//

#ifndef hairDrawTool_PrimitiveGenerator_h
#define hairDrawTool_PrimitiveGenerator_h

#include <maya/MPxNode.h>

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnComponentListData.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MRampAttribute.h>

#include <maya/MFnNurbsCurve.h>

#include <maya/MDagPath.h>
#include <maya/MFnData.h>
#include <maya/MDataHandle.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MGlobal.h>
#include <maya/MRampAttribute.h>

#include <maya/MPointArray.h>
#include <maya/MFloatPointArray.h>
#include <maya/MMatrixArray.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MVectorArray.h>
#include <maya/MEulerRotation.h>
#include <maya/MMatrix.h>
#include <maya/MFloatMatrix.h>
#include <maya/MPlugArray.h>

#include <maya/MFnDagNode.h>
#include <maya/MFnTransform.h>
#include <maya/MFnMesh.h>

#include <maya/MNodeMessage.h>
#include <maya/MCallbackIdArray.h>

#include <maya/MTime.h>

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <vector>




//#define M_PI = 3.14159265358979323846;

class primitiveGenerator : public MPxNode
{
public:

	primitiveGenerator();
	virtual ~primitiveGenerator();

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static  void*		creator();
	static  MStatus		initialize();
	void				postConstructor();

	MObject		        generateTubes();
	MObject             generateStrips();
	MPoint				rotate_point(float cx,float cy,float angle, MPoint p);
	std::vector<MMatrixArray>		calculateMatrix();
	MFloatArray			storeProfileCurveData(MRampAttribute a_segmentsAttribute, int segments, int segmentsLoop);
	
	MStatus				jiggle_calculate(MFloatVector goal);

	MStatus				displayOverride();

	static void			aboutToDeleteCB( MObject& node, MDGModifier& modifier, void* clientData );


	static MTypeId      id;

	static MObject		aOutMesh;
	static MObject      aInCurve;
	static MObject      aRefCurve;

	static MObject      aRadius;
	static MObject		aWidth;
	static MObject		aHeight;
	static MObject		aStrandWidth;
	static MObject		aStrandHeight;
	static MObject      aRotate;
	static MObject      aTwist;

	static MObject      aSides;
	static MObject      aSegments;
	static MObject		aSegmentsLoop;
	

	static MObject      aNumstrands;

	static MObject      aUseInputCurve;
	static MObject      aSmoothNormals;
	static MObject      aCapTop;
	static MObject		aAlingToUpVector;


	static MObject		aProfilePresets;
	static MObject		aStrandPresets;

	static MObject      aAutoSegments;
	static MObject      aAutoSegmentsRes;
	static MObject      aOnlyKnotSegmentsRes;

	static MObject      aCurveZOffset;

	static MObject      aStrandOffset;
	static MObject      aStrandOffsetRandom;
	static MObject		aRotationRandom;
	static MObject      aStrandThinning;
	static MObject      aStrandThinningRandomness;
	static MObject      aStrandThinningSeed;
	static MObject		aStrandCurl;
	static MObject		aStrandCurlWave;
	

	static MObject      aInLocAPos;
	static MObject      aInLocBPos;

	static MObject      aFirstUpVecX;
	static MObject      aFirstUpVecY;
	static MObject      aFirstUpVecZ;
	static MObject      aFirstUpVec;


	static MObject		aSegmentRamp;
	static MObject		aStrandOffsetRamp;
	static MObject		aTwistRamp;
	static MObject		aStrandCurlRamp;


	// UV
	static MObject      aCapUVSize;
	static MObject      aUOffset;
	static MObject      aVOffset;
	static MObject      aUOffsetCap;
	static MObject      aVOffsetCap;
	static MObject      aUWidth;
	static MObject      aVWidth;
	static MObject      aUVRotate;
	static MObject      aUVAutoV;
	static MObject      aUVAutoVMultiplier;

	// Overrides
	static MObject		aDisableBaseMeshOverride;

	// Jiggle

	static MObject		aJiggleEnabled;

	static MObject		aTime;
	static MObject		aStartFrame;
	static MObject		aJiggleAmount;
	static MObject		aDamping;
	static MObject		aStiffness;

private:

	double				m_r, m_width, m_height, m_strandWidth, m_strandHeight, m_rotate, m_twist, m_zOffset, m_strandOffset, m_strandOffsetRandom, m_rotationRandom, m_strandThinning, m_strandCurl, m_strandCurlWave, m_autoUV_V_mult, m_trandThinningRandomness;
	double				m_capUVsize, m_uWidth, m_vWidth, m_uOffset, m_vOffset, m_uOffsetCap, m_vOffsetCap, m_uvRotate;
	int					m_sides, m_segmentsLoop, m_segments, m_autoSegRes, m_type, m_numstrands, m_strandsLoop, m_strandThinningSeed;
	bool				m_autoSeg, m_smoothNorm, m_capTop, m_useProfile, m_segOnlyKnots, m_alingToUpVector, m_autoUV_V;
	short				m_profilePreset, m_strandPreset;

	MVector				m_firstUpVec;
	MObject				m_o_curve;
	MObject				m_o_curve_ref;

	MMatrix				m_curveMatrix;

	MDoubleArray		m_rndAr;
	MDoubleArray		m_lengthAr;

	MPointArray			m_profilePointsA;
	MFloatArray			m_segmentsProfileA;
	MFloatArray			m_strandOffsetProfileA;
	MFloatArray			m_twistProfileA;
	MFloatArray			m_stranCurlProfileA;

	MPoint				m_inLocA_pos;
	MPoint				m_inLocB_pos;

	// Jiggle

	bool				m_jiggleEnabled;

	MTime				m_currentTime;
	int					m_startFrame;
	float				m_damping;
	float				m_stiffness;
	MMatrix				m_parentInverse;
	float				m_jiggleAmount;

	MTime				m_previousTime;
	MPoint				m_currentPosition;
	MPoint				m_previousPosition;
	bool				m_init;
	bool				m_disableBaseMeshOverride;

	MFloatVector		m_jiggleVector;

	MCallbackIdArray	m_callbackIDs;

	std::vector<MMatrixArray> trMatrixA_vec;


};

#endif
