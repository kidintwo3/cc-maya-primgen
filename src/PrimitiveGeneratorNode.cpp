//
//  PrimitiveGenerator.cpp
//  hairDrawTool
//
//  Created by Hunyadi János on 2014. 12. 21..
//  Copyright (c) 2014. Janos Hunyadi. All rights reserved.
//

#include "PrimitiveGeneratorNode.h"
#include "PrimitiveGeneratorProfiles.h"

MTypeId     primitiveGenerator::id( 0x00123941 );
MObject     primitiveGenerator::aOutMesh;
MObject     primitiveGenerator::aInCurve;
MObject     primitiveGenerator::aRefCurve;
MObject     primitiveGenerator::aRadius;
MObject     primitiveGenerator::aWidth;
MObject     primitiveGenerator::aHeight;
MObject     primitiveGenerator::aRotate;
MObject     primitiveGenerator::aTwist;
MObject     primitiveGenerator::aSides;
MObject     primitiveGenerator::aSegments;
MObject		primitiveGenerator::aSegmentsLoop;
MObject     primitiveGenerator::aUseInputCurve;
MObject     primitiveGenerator::aSmoothNormals;
MObject     primitiveGenerator::aCapTop;
MObject     primitiveGenerator::aAlingToUpVector;

MObject     primitiveGenerator::aAutoSegments;
MObject     primitiveGenerator::aAutoSegmentsRes;
MObject     primitiveGenerator::aOnlyKnotSegmentsRes;

MObject		primitiveGenerator::aProfilePresets;

MObject     primitiveGenerator::aInLocAPos;
MObject     primitiveGenerator::aInLocBPos;

MObject     primitiveGenerator::aFirstUpVec;
MObject     primitiveGenerator::aFirstUpVecX;
MObject     primitiveGenerator::aFirstUpVecY;
MObject     primitiveGenerator::aFirstUpVecZ;

MObject		primitiveGenerator::aCurveZOffset;

MObject     primitiveGenerator::aSegmentRamp;

// UV
MObject     primitiveGenerator::aCapUVSize;
MObject     primitiveGenerator::aUOffset;
MObject     primitiveGenerator::aVOffset;
MObject     primitiveGenerator::aUOffsetCap;
MObject     primitiveGenerator::aVOffsetCap;
MObject     primitiveGenerator::aUWidth;
MObject     primitiveGenerator::aVWidth;
MObject     primitiveGenerator::aUVRotate;

// Overrides

MObject		primitiveGenerator::aDisableBaseMeshOverride;

// Jiggle


MObject     primitiveGenerator::aTime;
MObject		primitiveGenerator::aStartFrame;
MObject     primitiveGenerator::aJiggleEnabled;
MObject     primitiveGenerator::aJiggleAmount;
MObject     primitiveGenerator::aDamping;
MObject     primitiveGenerator::aStiffness;


// local attributes
//



primitiveGenerator::primitiveGenerator() 
{

}

primitiveGenerator::~primitiveGenerator() 
{
	MMessage::removeCallbacks(m_callbackIDs);
	m_callbackIDs.clear();
}

void* primitiveGenerator::creator()
{

	return new primitiveGenerator();
}

void primitiveGenerator::postConstructor(){
	m_init = false;

	MObject oSelf = thisMObject();

		// delete callback
	MCallbackId callbackID;
	callbackID = MNodeMessage::addNodeAboutToDeleteCallback(oSelf, aboutToDeleteCB, this);
	m_callbackIDs.append(callbackID);

}

// -------------------------------------------------

void primitiveGenerator::aboutToDeleteCB(MObject& node, MDGModifier& modifier, void* pUserPtr) 
{
	MFnDependencyNode nodeFn(node);
	MGlobal::displayInfo(MString("[PrimGen] About to delete callback for node: ") + nodeFn.name());


	// Find the output mesh connected to the node
	MPlug worldP = nodeFn.findPlug( "outMesh" );

	MFnDagNode mfDgN(worldP.node());

	MPlugArray destPlugs;
	worldP.connectedTo(destPlugs, false, true);

	if (destPlugs.length() != 0)
	{
		MPlug destPlug = destPlugs[0];

		mfDgN.setObject(destPlug.node());

		MFnDagNode mfDgN_transform(mfDgN.parent(0));


		MPlug p_out_overrideEnabled = mfDgN_transform.findPlug("overrideEnabled", false);
		p_out_overrideEnabled.setBool(false);

		MGlobal::displayInfo(MString("[PrimGen] Deleting / Setting output mesh overrides: ") + p_out_overrideEnabled.name() );
	}

	else
	{
		MGlobal::displayInfo(MString()+ "[PrimGen] Deleting / No connection, or wrong connection to output mesh: " );
	}

}

// -----------------------------------------------

MFloatArray primitiveGenerator::storeProfileCurveData(MRampAttribute a_segmentsAttribute, int segments, int segmentsLoop)
{
	MStatus status;

	MFloatArray curve_segments_values, curve_segments_values_loop;

	for (int i = 0; i < segments+1; i++)
	{
		float rampPosition = (1.0f / float(segments)) * float(i);
		float curveRampValue;
		a_segmentsAttribute.getValueAtPosition(rampPosition, curveRampValue, &status);
		CHECK_MSTATUS(status);
		curve_segments_values.append(curveRampValue);

	}

	return curve_segments_values;

}



MMatrixArray primitiveGenerator::calculateMatrix()
{


	MStatus status;
	MMatrixArray trMatrixA;

	// Spline
	if (m_type == 0)
	{



		MFnNurbsCurve curveFn(m_o_curve);


		double length = (curveFn.length() / double(m_segments) );

		// extra attributes

		MPointArray pA;

		MVector currentNormal = m_firstUpVec;
		currentNormal.normalize();

		MVector prevNormal;
		MMatrix rotMatrix;



		MDoubleArray knotsA;
		curveFn.getKnots(knotsA);

		MPointArray cvA;
		curveFn.getCVs(cvA, MSpace::kWorld);


		// If full spline


		for(int i=0; i < m_segments+1; i++){

			MPoint p;

			double param = curveFn.findParamFromLength( double(i) * length );
			curveFn.getPointAtParam(param, p, MSpace::kWorld );

			if (m_segOnlyKnots)
			{
				curveFn.getParamAtPoint(cvA[i], param, 1.0, MSpace::kWorld);
				p = cvA[i];
			}



			MVector tan = curveFn.tangent(param , MSpace::kWorld);
			tan.normalize();

			MVector cross1 = currentNormal^tan;
			cross1.normalize() ;



			MVector cross2 =  tan^cross1;

			if(m_alingToUpVector)
			{
				cross2 = m_firstUpVec;
			}


			cross2.normalize();
			currentNormal = cross2;

			p += cross1 * m_zOffset;

			double m[4][4] = {{tan.x, tan.y , tan.z, 0.0},
			{ cross1.x, cross1.y , cross1.z, 0.0},
			{cross2.x, cross2.y , cross2.z, 0.0},
			{p.x, p.y, p.z, 1.0}};



			rotMatrix = m;

			//rotMatrix += m_curveMatrix;
			//p *= m_curveMatrix;

			pA.append( p );

			// put everything back

			trMatrixA.append(rotMatrix);

		}





	}

	// A to B
	if (m_type == 1)
	{


		MMatrix rotMatrix;

		MVector upVec(0.0, 1.0, 0.0);


		MVector ab = m_inLocB_pos - m_inLocA_pos;

		double abLength = ab.length();
		double step = abLength / double(m_segments);

		MVector ab_dir = ab;
		ab_dir.normalize();

		MVector xDir = ab_dir;
		MVector yDir = upVec ^ xDir;
		yDir.normalize();
		MVector zDir = xDir ^ yDir;
		zDir.normalize();


		// Calculate jiggle
		//
		if (m_jiggleEnabled)
		{

			MFloatVector cp = ab;
			cp.normalize();
			MFloatVector cv = m_inLocA_pos + (0.5f * cp);

			jiggle_calculate(cv);


			for(int i=0; i <= m_segments; i++)
			{
				MPoint p(m_inLocA_pos + xDir * step * double(i));

				MPoint p_jiggle( m_jiggleVector + xDir * step * double(i) );
				MVector pOff = p_jiggle - p;

				double mult = 1.0;

				if (i < int(m_segmentsProfileA.length()))
				{
					mult = m_segmentsProfileA[i];
				}

				MPoint pFin = p - (pOff * mult);

				if (!m_init)
				{
					pFin = p;
				}


				double m[4][4] = {{xDir.x, xDir.y, xDir.z, 0.0}, {yDir.x, yDir.y, yDir.z, 0.0}, {zDir.x, zDir.y, zDir.z, 0.0}, {pFin.x, pFin.y, pFin.z, 1.0}};

				rotMatrix = m;

				// put everything back

				trMatrixA.append(rotMatrix);

			}

		}

		// No jiggle
		//
		if (!m_jiggleEnabled)
		{

			for(int i=0; i <= m_segments; i++)
			{

				MPoint p(  m_inLocA_pos + xDir * step * double(i)   );	

				double m[4][4] = {{xDir.x, xDir.y, xDir.z, 0.0}, {yDir.x, yDir.y, yDir.z, 0.0}, {zDir.x, zDir.y, zDir.z, 0.0}, {p.x, p.y, p.z, 1.0}};

				rotMatrix = m;

				// put everything back

				trMatrixA.append(rotMatrix);

			}

		}
	}

	return trMatrixA;

}


MStatus primitiveGenerator::jiggle_calculate(MFloatVector goal)
{

	if (!m_init)
	{
		m_previousTime = m_currentTime;
		m_currentPosition = goal;
		m_previousPosition = goal;
		m_init = true;
	}


	float timeDifference = float(m_currentTime.value() - m_previousTime.value());

	if (timeDifference > 1.0 || timeDifference < 0.0 || m_currentTime.value() < m_startFrame) 
	{

		// MGlobal::displayInfo(MString() + "m_currentTime.value(): " + m_currentTime.value());

		m_init = false;
		m_previousTime = m_currentTime;

		m_jiggleVector = MFloatVector(0.0,0.0,0.0);

		return MStatus::kSuccess;
	}



	MVector velocity = (m_currentPosition - m_previousPosition) * (1.0 - m_damping);
	MPoint newPosition = m_currentPosition + velocity;
	MVector goalForce = (goal - newPosition) * m_stiffness;
	newPosition += goalForce;

	//Store the states for the next computation
	m_previousPosition = MPoint(m_currentPosition);
	m_currentPosition = MPoint(newPosition);
	m_previousTime = MTime(m_currentTime);

	newPosition = goal + ((MVector(newPosition) - goal) * m_jiggleAmount);


	// Put in the output local space
	//newPosition *= parentInverse;


	m_jiggleVector = MFloatVector(newPosition.x, newPosition.y, newPosition.z);


	return MStatus::kSuccess;
}

MObject primitiveGenerator::generateStrips(){

	MStatus status;
	MMatrixArray trMatrixA = calculateMatrix();

	MDoubleArray profile;

	MPointArray pA;

	int num_verts = 0;
	int num_faces = 0;

	num_verts = ( 2 * m_segments) + 2 ;
	num_faces = m_segments;


	// double angleRot = m_rotate / 180.0 * M_PI;




	for (int i=0; i < m_segments+1; i++) 
	{


		double angleRot = m_rotate / 180.0 * M_PI;
		angleRot += m_twist*i / double(m_segments);

		MTransformationMatrix trM(trMatrixA[i]);
		double scale[3] = {1.0,m_width,0.0};
		trM.rotateBy(MEulerRotation(angleRot,0.0,0.0),MSpace::kObject);
		trM.setScale(scale,MSpace::kObject);

		double rad = m_r;
		rad *= m_segmentsProfileA[i];

		MPoint p1 = MPoint( 0.0, -rad, 0.0 );
		MPoint p2 = MPoint( 0.0, rad, 0.0 );

		pA.append(  MFloatPoint( p1 * trM.asMatrix()) );
		pA.append(  MFloatPoint( p2 * trM.asMatrix()) );

	}

//	int len = pA.length();

	// Facecounts CAPS
	MIntArray faceCounts;
	MIntArray faceConnects;


	for (int i=0; i < m_segments; i++) {

		faceCounts.append(4);
	}


	int v00 = 0;
	int v01 = 2;
	int v02 = 3;
	int v03 = 1;

	for (int i=0; i < m_segments; i++) {

		faceConnects.append(v00);
		faceConnects.append(v01);
		faceConnects.append(v02);
		faceConnects.append(v03);
		v00 = v01;
		v01 = v01 + 2;
		v03 = v02;
		v02 = v01 + 1;
	}


	// ------------------------------------------------------------------------------------------
	// UV

	MIntArray           uvCounts;
	MIntArray           uvIds;
	MFloatArray         uArray;
	MFloatArray         vArray;

	double uv_width = m_uWidth;
	double uv_height = m_vWidth;
	double u_offset = m_uOffset;
	double v_offset = m_vOffset;

	if (uv_height == 0.0) {
		uv_height = 0.001;
	}

	if (uv_width == 0.0) {
		uv_width = 0.001;
	}

	double u = uv_height / double(m_segments);
	double vh = 0.0;
	double rotAxis = (m_uvRotate + 180.000)  * ( M_PI / 180.0 );

	MPoint rotUVP;

	for (int i=0; i<m_segments+1; i++) {

		float uO = uv_width + u_offset;
		float vO = 0.0 + v_offset + vh;

		rotUVP = rotate_point(uO,vO,rotAxis, MPoint(m_uOffset,m_vOffset,0.0));
		uArray.append(rotUVP.x + m_uOffset);
		vArray.append(rotUVP.y + m_vOffset);

		float u1 = 0.0 + u_offset;
		float v1 = 0.0 + v_offset + vh;

		rotUVP = rotate_point(u1,v1,rotAxis, MPoint(m_uOffset,m_vOffset,0.0));
		uArray.append(rotUVP.x + m_uOffset);
		vArray.append(rotUVP.y + m_vOffset);

		vh+=u;
	}

	uvCounts = faceCounts;
	uvIds = faceConnects;


	// ------------------------------------------------------------------------------------------
	// Create mesh

	MFnMeshData meshDataFn;
	MObject newMeshData = meshDataFn.create();
	MFnMesh meshFn;
	meshFn.create( num_verts, num_faces, pA, faceCounts, faceConnects, uArray,vArray, newMeshData, &status );
	CHECK_MSTATUS( status );

	status = meshFn.assignUVs(uvCounts,uvIds);
	CHECK_MSTATUS( status );

	for (int i = 0; i < meshFn.numEdges(); i++)
	{
		if (m_smoothNorm) { meshFn.setEdgeSmoothing(i, true);	}
		if (!m_smoothNorm) { meshFn.setEdgeSmoothing(i, false);	}
	}

	return newMeshData;

}







MObject primitiveGenerator::generateTubes()
{

	MStatus status;

	MMatrixArray trMatrixA = calculateMatrix();

	// ------------------------------------------------------------------------------------------
	//create verts
	//

	if (m_sides == 0) 
	{
		m_sides = 2;
	}

	MPointArray pA;

	int num_verts = 0;
	double x,z;
	double deg = 360.0 / double(m_sides);

	for (int i=0; i < m_segments+1; i++) {


		for (int j=0; j< m_sides; j++) {

			double angle = deg * j / 180.0  * M_PI;
			double angleRot = m_rotate / 180.0 * M_PI;
			angleRot += m_twist*i / double(m_segments);

			if (m_useProfile)
			{
				x = m_profilePointsA[j].x * m_r;
				z = m_profilePointsA[j].z * m_r;

				x *= m_segmentsProfileA[i];
				z *= m_segmentsProfileA[i];

				MPoint pnt( 0.0, x, z );
				MTransformationMatrix trM;
				trM.rotateBy(MEulerRotation(angleRot,0.0,0.0),MSpace::kObject);

				double scale[3] = {1.0,m_width,m_height};
				trM.setScale(scale,MSpace::kObject);



				pnt *= trM.asRotateMatrix();
				pnt *= trMatrixA[i];
				

				pA.append( MFloatPoint( pnt ) );
			}

			if (!m_useProfile)
			{

				x  = cos( angle ) * m_r;
				z  = sin( angle ) * m_r;


				x *= m_segmentsProfileA[i];
				z *= m_segmentsProfileA[i];

				MPoint pnt( 0.0, x, z );

				MTransformationMatrix trM(trMatrixA[i]);
				double scale[3] = {1.0,m_width,m_height};
				trM.rotateBy(MEulerRotation(angleRot,0.0,0.0),MSpace::kObject);
				trM.setScale(scale,MSpace::kObject);


				MFloatPoint outP = MFloatPoint( (pnt * trM.asMatrix()));

				pA.append( outP );
			}



			num_verts += 1;
		}

	}

	// numFaces
	int num_faces = m_sides * (m_segments) + 2;
	if (!m_capTop){ num_faces = m_sides * (m_segments); }


	// Facecounts CAPS
	MIntArray faceCounts;

	for (int i=0; i<m_segments; i++) 
	{

		for (int j=0; j<m_sides; j++) 
		{
			faceCounts.append(4);
		}

	}

	if (m_capTop)
	{ 
		faceCounts.append(m_sides);
		faceCounts.append(m_sides);
	}



	// Faceconnects SIDES
	MIntArray faceConnects;
	for (int j=0; j<m_segments; j++) 
	{

		for (int i=0; i<m_sides; i++) 
		{

			int v0 = m_sides * j + i;
			int v1 = m_sides * (j+1) + i;
			int v2 = m_sides * (j+1) + (i + 1);
			int v3 = m_sides * (j) + (i + 1);

			if (i == m_sides-1)
			{
				v2 = m_sides * j + (i + 1);
				v3 = m_sides * j;
			}

			faceConnects.append(v0);
			faceConnects.append(v3);
			faceConnects.append(v2);
			faceConnects.append(v1);

		}
	}

	if (m_capTop)
	{ 

		for (int i=0; i<m_sides; i++) {

			int f = ((m_segments) * m_sides) + i;
			faceConnects.append(f);
		}

		for (int i=m_sides-1; i > -1 ; i--) {

			faceConnects.append(i);
		}
	}

	// ------------------------------------------------------------------------------------------
	// UV

	MIntArray           uvCounts;
	MIntArray           uvIds;
	MFloatArray         uArray;
	MFloatArray         vArray;

	int v1,v2,v3,v4,lastUV;
	int counter = m_sides;

	for (int j=0; j<m_segments; j++) {

		for (int i=0; i<m_sides; i++) {

			v1 = i + 0 + counter;
			v2 = i + 1 + counter;
			v3 = i + 2 + m_sides + counter;
			v4 = i + 1 + m_sides + counter;

			uvIds.append(v1);
			uvIds.append(v2);
			uvIds.append(v3);
			uvIds.append(v4);

			lastUV = v4 + m_sides;
		}

		counter +=m_sides+1;
	}

	counter = 0;

	if (m_capTop)
	{ 

		for (int i=0; i<m_sides; i++) {
			uvIds.append(i);
		}

		uvIds.append(lastUV+1);
		for ( int i = lastUV+1-m_sides+1; i<lastUV+1; i++) {
			uvIds.append(i);
		}
	}

	// - UV uvCounts

	for (int i=0; i< m_sides * m_segments; i++) {

		uvCounts.append(4);
	}

	if (m_capTop)
	{ 
		uvCounts.append(m_sides);
		uvCounts.append(m_sides);
	}

	// - UV array

	double u,v;

	double angle;

	for (int i=0; i<m_sides; i++) {

		double deg = 360.0 / double(m_sides);

		if ( i != 0) {
			angle = deg * double(i) / 180.0  * M_PI;
		}

		else {
			angle = 0.0;
		}
		double angleRot = m_rotate / 180.0 * M_PI;

		//angleRot += twist*i / segments;

		u = cos( angle ) * m_capUVsize;
		v = sin( angle ) * m_capUVsize;

		if (!m_useProfile){

			u  = cos( angle + angleRot ) * m_capUVsize;
			v  = sin( angle + angleRot ) * m_capUVsize;
		}


		if (m_useProfile)
		{ 

			u = m_profilePointsA[i].x * m_capUVsize;
			v = m_profilePointsA[i].z * m_capUVsize;


		}


		uArray.append(u + m_uOffsetCap);
		vArray.append(v + m_vOffsetCap);

	}


	for (int i=0; i < m_segments + 1; i++) {

		for (int j=0; j < m_sides +1; j++) {
			u = double(j) / (m_sides * (1.0 / m_uWidth));
			v = double(i) / (m_segments * (1.0 / m_vWidth));

			double uO = u + m_uOffset;
			double vO = v + m_vOffset;

			double rotAxis = (m_uvRotate + 180.000)  * ( M_PI / 180.0 );

			MPoint rotUVP = rotate_point(uO,vO,rotAxis, MPoint(m_uOffset,m_vOffset,0.0));

			/*uArray.append(u + uOffset);
			vArray.append(v + vOffset);*/

			uArray.append(rotUVP.x + m_uOffset);
			vArray.append(rotUVP.y + m_vOffset);

		}
	}

	if (m_capTop)
	{ 

		for (int i=0; i<m_sides; i++) {

			double deg = 360.0 / double(m_sides);

			double angle = deg * double(i) / 180.0  * M_PI;
			double angleRot = m_rotate / 180.0 * M_PI;

			if (!m_useProfile){

				u  = cos( angle + angleRot ) * m_capUVsize;
				v  = sin( angle + angleRot ) * m_capUVsize;

			}


			if (m_useProfile)
			{ 

				u = m_profilePointsA[i].x * m_capUVsize;
				v = m_profilePointsA[i].z * m_capUVsize;
			}

			u += m_capUVsize*2.0;

			uArray.append(u + m_uOffsetCap);
			vArray.append(v + m_vOffsetCap);
		}
	}

	// ------------------------------------------------------------------------------------------
	// Create mesh

	MFnMeshData meshDataFn;
	MObject newMeshData = meshDataFn.create();
	MFnMesh meshFn;
	meshFn.create( num_verts, num_faces, pA, faceCounts, faceConnects, uArray,vArray, newMeshData, &status );
	CHECK_MSTATUS( status );
	status = meshFn.assignUVs(uvCounts,uvIds);
	CHECK_MSTATUS( status );

	for (int i = 0; i < meshFn.numEdges(); i++)
	{
		if (m_smoothNorm) { meshFn.setEdgeSmoothing(i, true);	}
		if (!m_smoothNorm) { meshFn.setEdgeSmoothing(i, false);	}
	}

	return newMeshData;

}

MPoint primitiveGenerator::rotate_point(float cx,float cy,float angle, MPoint p)
{
	float s = sin(angle);
	float c = cos(angle);

	// translate point back to origin:
	p.x -= cx;
	p.y -= cy;

	// rotate point
	float xnew =  p.x * c + p.y * s;
	float ynew = -p.x * s + p.y * c;

	// translate point back:
	p.x = xnew ;
	p.y = ynew ;

	return p;
}

MStatus primitiveGenerator::displayOverride()
{

	MStatus status;

	MPlug p_outMesh( this->thisMObject(), aOutMesh );

	if (p_outMesh.isConnected())
	{
		// -----------------------------------------------
		// Collect output plug mesh's name

		MPlugArray outputs_plugArr;
		p_outMesh.connectedTo(outputs_plugArr, false, true, &status);
		CHECK_MSTATUS_AND_RETURN_IT(status);

		if (outputs_plugArr.length() > 0)
		{
			MPlug outMeshPlug_shape = outputs_plugArr[0];
			MFnDependencyNode outMesh_dn(outMeshPlug_shape.node());


			MPlug p_out_overrideEnabled = outMesh_dn.findPlug("overrideEnabled", false, &status);
			CHECK_MSTATUS_AND_RETURN_IT(status);
			p_out_overrideEnabled.setBool(m_disableBaseMeshOverride);

			MPlug p_out_overrideDisplayType = outMesh_dn.findPlug("overrideDisplayType", false, &status);
			CHECK_MSTATUS_AND_RETURN_IT(status);
			p_out_overrideDisplayType.setInt(2);
		}

	}


	return MS::kSuccess;
}



MStatus primitiveGenerator::compute( const MPlug& plug, MDataBlock& data )
{
	MStatus status;

	MPlug p_incurve( this->thisMObject(), aInCurve );
	MPlug p_refcurve( this->thisMObject(), aRefCurve );
	MPlug p_inLocA( this->thisMObject(), aInLocAPos );
	MPlug p_inLocB( this->thisMObject(), aInLocBPos );
	

	if ( plug != aOutMesh)
	{
		return MS::kUnknownParameter;
	}

	MDataHandle hInCurve = data.inputValue(aInCurve, &status);
	CHECK_MSTATUS_AND_RETURN_IT( status );

	MDataHandle hRefCurve = data.inputValue(aRefCurve, &status);
	CHECK_MSTATUS_AND_RETURN_IT( status );

//	MDataHandle hInLocA = data.inputValue(aInLocAPos, &status);
//	CHECK_MSTATUS_AND_RETURN_IT( status );
//
//	MDataHandle hInLocB = data.inputValue(aInLocBPos, &status);
//	CHECK_MSTATUS_AND_RETURN_IT( status );

	MDataHandle hOutput = data.outputValue( aOutMesh, &status );
	CHECK_MSTATUS_AND_RETURN_IT( status );

	m_o_curve = hInCurve.asNurbsCurve();
	m_o_curve_ref = hRefCurve.asNurbsCurve();


	// ------------------------------------------------------------------------------------------
	m_r						= data.inputValue( aRadius, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_width					= data.inputValue( aWidth, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_height				= data.inputValue( aHeight, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_rotate				= data.inputValue( aRotate, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_twist					= data.inputValue( aTwist, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_zOffset				= data.inputValue( aCurveZOffset, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_sides					= data.inputValue( aSides, &status ).asInt();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_segmentsLoop			= data.inputValue( aSegmentsLoop, &status ).asInt();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_segments				= data.inputValue( aSegments, &status ).asInt();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_autoSegRes			= data.inputValue( aAutoSegmentsRes, &status).asInt();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_autoSeg				= data.inputValue( aAutoSegments, &status).asBool();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_segOnlyKnots				= data.inputValue( aOnlyKnotSegmentsRes, &status).asBool();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_smoothNorm			= data.inputValue( aSmoothNormals, &status ).asBool();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_capTop				= data.inputValue( aCapTop, &status ).asBool();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_alingToUpVector		= data.inputValue( aAlingToUpVector, &status ).asBool();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_profilePreset			= data.inputValue(aProfilePresets, &status).asShort();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_firstUpVec			= data.inputValue(aFirstUpVec, &status).asVector();
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// Uv

	m_capUVsize				= data.inputValue( aCapUVSize, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_uWidth                = data.inputValue( aUWidth, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_vWidth                = data.inputValue( aVWidth, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_uOffset               = data.inputValue( aUOffset, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_vOffset               = data.inputValue( aVOffset, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_uOffsetCap            = data.inputValue( aUOffsetCap, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_vOffsetCap            = data.inputValue( aVOffsetCap, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_uvRotate				= data.inputValue( aUVRotate, &status ).asDouble();
	CHECK_MSTATUS_AND_RETURN_IT(status);


	// Override
	m_disableBaseMeshOverride = data.inputValue(aDisableBaseMeshOverride, &status).asBool();
	CHECK_MSTATUS_AND_RETURN_IT(status);

	// Jiggle


	m_jiggleEnabled			= data.inputValue(aJiggleEnabled, &status).asBool();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_currentTime			= data.inputValue(aTime, &status).asTime();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_startFrame			= data.inputValue(aStartFrame, &status).asInt();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_damping				= data.inputValue(aDamping, &status).asFloat();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_stiffness				= data.inputValue(aStiffness, &status).asFloat();
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_jiggleAmount			= data.inputValue(aJiggleAmount, &status).asFloat();
	CHECK_MSTATUS_AND_RETURN_IT(status);

	m_jiggleVector			= MFloatVector::zero;

	m_type					= 0;



	// Failsafe check
	if (m_segments <= 1) { m_segments = 1; }


	// AB locator
	MMatrix m_inLocA_posM   = data.inputValue( aInLocAPos ).asMatrix();
	MMatrix m_inLocB_posM   = data.inputValue( aInLocBPos ).asMatrix();

	

	MTransformationMatrix m_inLocA_posMat(m_inLocA_posM);
	m_inLocA_pos = m_inLocA_posMat.getTranslation(MSpace::kWorld);
	MTransformationMatrix m_inLocB_posMat(m_inLocB_posM);
	m_inLocB_pos = m_inLocB_posMat.getTranslation(MSpace::kWorld);
	if ( p_inLocA.isConnected() && p_inLocB.isConnected()) { m_type = 1; }
	if ( m_type == 1) { if (!p_inLocA.isConnected() && !p_inLocB.isConnected()) { m_inLocA_pos = MPoint(-5.0,0.0,0.0); m_inLocB_pos = MPoint( 5.0,0.0,0.0);} }


	// Auto segments
	if (m_autoSeg)
	{

		// If curve
		if (m_type == 0)
		{
			MFnNurbsCurve mfC(m_o_curve);
			double curveLen = mfC.length();
			m_segments = curveLen * m_autoSegRes;
		}

		// If A-B
		if (m_type == 1)
		{
			double curveLen = m_inLocA_pos.distanceTo(m_inLocB_pos);
			m_segments = curveLen * m_autoSegRes;
		}

		// Failsafe check again... just to be sure
		if (m_segments <= 1) { m_segments = 1; }

	}

	if (m_segOnlyKnots && m_type == 0)
	{
		MFnNurbsCurve mfC(m_o_curve);
		m_segments = mfC.numCVs()-1;
	}


	//MObject m_inCurve = data.inputValue(aInCurve, &status).asNurbsCurve();
	//CHECK_MSTATUS_AND_RETURN_IT(status);

	//// Curve matrix
	MFnNurbsCurve curve_fn(m_o_curve);
	m_curveMatrix = curve_fn.dagPath().inclusiveMatrix();


	// Ramp attribute
	MRampAttribute a_segmentsAttribute(this->thisMObject(), aSegmentRamp, &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);
	m_segmentsProfileA = storeProfileCurveData(a_segmentsAttribute, m_segments, m_segmentsLoop);


	//
	m_profilePointsA.clear();

	if (m_profilePreset == 1) { m_sides = 36; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_cross[i][0]  , m_profile_cross[i][1]  , m_profile_cross[i][2]  ), i); }}
	else if (m_profilePreset == 2) { m_sides = 12; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_square[i][0]  , m_profile_square[i][1]  , m_profile_square[i][2]  ), i); }}
	else if (m_profilePreset == 3) { m_sides = 24; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_uShape[i][0]  , m_profile_uShape[i][1]  , m_profile_uShape[i][2]  ), i); }}
	else if (m_profilePreset == 4) { m_sides = 17; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_hexagon[i][0]  , m_profile_hexagon[i][1]  , m_profile_hexagon[i][2]  ), i); }}
	else if (m_profilePreset == 5) { m_sides = 24; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_octagon[i][0]  , m_profile_octagon[i][1]  , m_profile_octagon[i][2]  ), i); }}
	else if (m_profilePreset == 6) { m_sides = 95; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_cogWheel[i][0]  , m_profile_cogWheel[i][1]  , m_profile_cogWheel[i][2]  ), i); }}
	else if (m_profilePreset == 7) { m_sides = 16; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_leatherBeltThick[i][0]  , m_profile_leatherBeltThick[i][1]  , m_profile_leatherBeltThick[i][2]  ), i); }}
	else if (m_profilePreset == 8) { m_sides = 12; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_leatherBeltThin[i][0]  , m_profile_leatherBeltThin[i][1]  , m_profile_leatherBeltThin[i][2]  ), i); }}
	else if (m_profilePreset == 9) { m_sides = 88; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_drill[i][0]  , m_profile_drill[i][1]  , m_profile_drill[i][2]  ), i); }}
	else if (m_profilePreset == 10) { m_sides = 18; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_leatherLShape[i][0]  , m_profile_leatherLShape[i][1]  , m_profile_leatherLShape[i][2]  ), i); }}
	else if (m_profilePreset == 11) { m_sides = 26; m_profilePointsA.setLength(m_sides); for (int i = 0; i < m_sides; i++) { m_profilePointsA.set(MPoint(m_profile_ropeShape[i][0]  , m_profile_ropeShape[i][1]  , m_profile_ropeShape[i][2]  ), i); }}




	MObject newMeshData;
	if (m_profilePreset == 0) { m_useProfile = false; }
	else					  { m_useProfile = true;  }


	// Reference curve
	if (m_profilePreset == 0) 
	{
		if (p_refcurve.isConnected())
		{
			MFnNurbsCurve mfC_ref(m_o_curve_ref);
			MPointArray ref_cvA;

			mfC_ref.getCVs(ref_cvA);

			if (ref_cvA.length() < 3)
			{
				return::MStatus::kSuccess;
			}

			// reverse array
			MPointArray rev_cvA;
			int co = ref_cvA.length() -1;
			for (unsigned i = 0; i <  ref_cvA.length() ; i++)
			{
				rev_cvA.append(ref_cvA[co]);
				co -= 1;
			}

			m_useProfile = true;
			m_sides = ref_cvA.length();
			m_profilePointsA.setLength(m_sides);

			for (int i = 0; i < m_sides; i++)
			{
				m_profilePointsA.set(MPoint(rev_cvA[i][0]  , rev_cvA[i][1]  , rev_cvA[i][2]  ), i); 
			}


		}
	}

	// ------------------------------------------------------------------------------------------

	if (m_sides <= 2)
	{
		newMeshData = generateStrips();
	}

	else
	{
		newMeshData = generateTubes();
	}




	hOutput.set( newMeshData );

	status = displayOverride();
	CHECK_MSTATUS_AND_RETURN_IT(status);

	return MS::kSuccess;
}


MStatus primitiveGenerator::initialize()
{
	// local attribute initialization

	MStatus status;


	MFnTypedAttribute		tAttr;
	MFnNumericAttribute		nAttr;
	MFnMatrixAttribute		mAttr;
	MFnCompoundAttribute	cAttr;
	MRampAttribute			rAttr;
	MFnEnumAttribute		eAttr;
	MFnUnitAttribute        uAttr;

	primitiveGenerator::aOutMesh = tAttr.create( "outMesh", "outMesh", MFnData::kMesh );
	tAttr.setChannelBox(false);
	tAttr.setStorable(false);
	tAttr.setKeyable(false);
	tAttr.setChannelBox( false );
	addAttribute( primitiveGenerator::aOutMesh );

	primitiveGenerator::aInCurve = tAttr.create( "inCurve", "inCurve", MFnData::kNurbsCurve );
	tAttr.setChannelBox(false);
	tAttr.setWritable(true);
	tAttr.setReadable(false);
	tAttr.setStorable(false);
	tAttr.setKeyable(true);
	addAttribute( primitiveGenerator::aInCurve );

	primitiveGenerator::aRefCurve = tAttr.create( "refCurve", "refCurve", MFnData::kNurbsCurve );
	tAttr.setChannelBox(false);
	tAttr.setWritable(true);
	tAttr.setReadable(false);
	tAttr.setStorable(false);
	tAttr.setKeyable(true);
	addAttribute( primitiveGenerator::aRefCurve );

	primitiveGenerator::aProfilePresets = eAttr.create( "profilePresets", "profilePresets", 0);
	eAttr.setStorable(true);
	eAttr.addField("Custom", 0);
	eAttr.addField("Cross", 1);
	eAttr.addField("Square", 2);
	eAttr.addField("U Shape", 3);
	eAttr.addField("Hexagon", 4);
	eAttr.addField("Octagon", 5);
	eAttr.addField("Cog Wheel", 6);
	eAttr.addField("Leather Belt Thick", 7);
	eAttr.addField("Leather Belt Thin", 8);
	eAttr.addField("Drilled Metal", 9);
	eAttr.addField("L Shape", 10);
	eAttr.addField("Rope", 11);

	eAttr.setDefault(1);
	addAttribute( primitiveGenerator::aProfilePresets );

	primitiveGenerator::aRadius = nAttr.create( "radius", "radius", MFnNumericData::kDouble );
	nAttr.setDefault( 1.0 );
	nAttr.setMin(0.0);
	nAttr.setSoftMax(10.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aRadius );

	primitiveGenerator::aWidth = nAttr.create( "width", "width", MFnNumericData::kDouble );
	nAttr.setDefault( 1.0 );
	nAttr.setMin(0.0);
	nAttr.setSoftMax(10.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aWidth );

	primitiveGenerator::aHeight = nAttr.create( "height", "height", MFnNumericData::kDouble );
	nAttr.setDefault( 1.0 );
	nAttr.setMin(0.0);
	nAttr.setSoftMax(10.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aHeight );


	primitiveGenerator::aRotate = nAttr.create( "rotate", "rotate", MFnNumericData::kDouble );
	nAttr.setDefault( 0.0 );
	nAttr.setMin(0.0);
	nAttr.setMax(360.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aRotate );

	primitiveGenerator::aTwist = nAttr.create( "twist", "twist", MFnNumericData::kDouble );
	nAttr.setDefault( 0.0 );
	nAttr.setSoftMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aTwist );

	primitiveGenerator::aSides = nAttr.create( "sides", "sides", MFnNumericData::kInt );
	nAttr.setDefault( 3 );
	nAttr.setMin(2);
	nAttr.setSoftMax(50);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aSides );

	primitiveGenerator::aSegments = nAttr.create( "segments", "segments", MFnNumericData::kInt );
	nAttr.setDefault( 5 );
	nAttr.setMin(1);
	nAttr.setSoftMax(50);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aSegments );

	primitiveGenerator::aSegmentsLoop = nAttr.create( "segmentsLoop", "segmentsLoop", MFnNumericData::kInt );
	nAttr.setDefault( 1 );
	nAttr.setMin( 1 );
	nAttr.setSoftMax( 5 );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aSegmentsLoop );

	primitiveGenerator::aSegmentRamp = rAttr.createCurveRamp("segmentsRamp", "segmentsRamp");
	addAttribute(aSegmentRamp);

	primitiveGenerator::aUseInputCurve = nAttr.create( "useInputCurve", "useInputCurve", MFnNumericData::kBoolean );
	nAttr.setDefault( false );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aUseInputCurve );

	primitiveGenerator::aSmoothNormals = nAttr.create( "smoothNormals", "smoothNormals", MFnNumericData::kBoolean );
	nAttr.setDefault( true );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aSmoothNormals );

	primitiveGenerator::aCapTop = nAttr.create( "capTop", "capTop", MFnNumericData::kBoolean );
	nAttr.setDefault( true );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aCapTop );

	primitiveGenerator::aAlingToUpVector = nAttr.create( "alingToUpVector", "alingToUpVector", MFnNumericData::kBoolean );
	nAttr.setDefault( false );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aAlingToUpVector );

	primitiveGenerator::aOnlyKnotSegmentsRes = nAttr.create( "autoSegmentsKnotsOnly", "autoSegmentsKnotsOnly", MFnNumericData::kBoolean );
	nAttr.setDefault( false );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aOnlyKnotSegmentsRes );

	primitiveGenerator::aAutoSegments = nAttr.create( "autoSegments", "autoSegments", MFnNumericData::kBoolean );
	nAttr.setDefault( false );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aAutoSegments );

	primitiveGenerator::aAutoSegmentsRes = nAttr.create( "autoSegmentsMultiplier", "autoSegmentsMultiplier", MFnNumericData::kInt );
	nAttr.setDefault( 2 );
	nAttr.setMin( 1 );
	nAttr.setSoftMax( 5 );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aAutoSegmentsRes );

	primitiveGenerator::aCurveZOffset = nAttr.create( "curveZOffset", "curveZOffset", MFnNumericData::kDouble );
	nAttr.setDefault( 0.0 );
	nAttr.setSoftMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aCurveZOffset );


	aInLocAPos = mAttr.create( "locatorAPos", "locatorAPos", MFnMatrixAttribute::kDouble );
	mAttr.setChannelBox(false);
	mAttr.setWritable(true);
	mAttr.setReadable(false);
	mAttr.setStorable(false);
	mAttr.setKeyable(false);
	addAttribute( aInLocAPos );

	aInLocBPos = mAttr.create( "locatorBPos", "locatorBPos", MFnMatrixAttribute::kDouble );
	mAttr.setChannelBox(false);
	mAttr.setWritable(true);
	mAttr.setReadable(false);
	mAttr.setStorable(false);
	mAttr.setKeyable(false);
	addAttribute( aInLocBPos );


	aFirstUpVecX = nAttr.create("firstUpVecX","fux",MFnNumericData::kDouble,0);
	nAttr.setChannelBox(false);
	nAttr.setStorable(true);
	nAttr.setKeyable(true);
	addAttribute(aFirstUpVecX);

	aFirstUpVecY = nAttr.create("firstUpVecY","fuy",MFnNumericData::kDouble,1);
	nAttr.setChannelBox(false);
	nAttr.setStorable(true);
	nAttr.setKeyable(true);
	addAttribute(aFirstUpVecY);


	aFirstUpVecZ = nAttr.create("firstUpVecZ","fuz",MFnNumericData::kDouble,0);
	nAttr.setChannelBox(false);
	nAttr.setStorable(true);
	nAttr.setKeyable(true);
	addAttribute(aFirstUpVecZ);


	aFirstUpVec = cAttr.create("firstUpVec","fu");
	cAttr.setChannelBox(false);
	cAttr.addChild(aFirstUpVecX);
	cAttr.addChild(aFirstUpVecY);
	cAttr.addChild(aFirstUpVecZ);
	addAttribute(aFirstUpVec);

	// UV

	primitiveGenerator::aCapUVSize = nAttr.create( "capUvSize", "capUvSize", MFnNumericData::kDouble );
	nAttr.setDefault( 0.5 );
	nAttr.setMin(0.0);
	nAttr.setSoftMax(0.5);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aCapUVSize );

	primitiveGenerator::aUOffset = nAttr.create( "uOffset", "uOffset", MFnNumericData::kDouble );
	nAttr.setDefault( 0.0 );
	nAttr.setSoftMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aUOffset );

	primitiveGenerator::aVOffset = nAttr.create( "vOffset", "vOffset", MFnNumericData::kDouble );
	nAttr.setDefault( 0.0 );
	nAttr.setSoftMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aVOffset );

	primitiveGenerator::aUOffsetCap = nAttr.create( "uOffsetCap", "uOffsetCap", MFnNumericData::kDouble );
	nAttr.setDefault( 1.5 );
	nAttr.setSoftMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aUOffsetCap );

	primitiveGenerator::aVOffsetCap = nAttr.create( "vOffsetCap", "vOffsetCap", MFnNumericData::kDouble );
	nAttr.setDefault( 0.5 );
	nAttr.setSoftMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aVOffsetCap );

	primitiveGenerator::aUWidth = nAttr.create( "uWidth", "uWidth", MFnNumericData::kDouble );
	nAttr.setDefault( 1.0 );
	nAttr.setMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aUWidth );

	primitiveGenerator::aVWidth = nAttr.create( "vWidth", "vWidth", MFnNumericData::kDouble );
	nAttr.setDefault( 1.0 );
	nAttr.setMin(0.0);
	nAttr.setSoftMax(1.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aVWidth );

	primitiveGenerator::aUVRotate = nAttr.create( "uvRotate", "uvRotate", MFnNumericData::kDouble );
	nAttr.setDefault( 0.0 );
	nAttr.setMin(0.0);
	nAttr.setMax(360.0);
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aUVRotate );

	// Overrides

	primitiveGenerator::aDisableBaseMeshOverride = nAttr.create("baseMeshDisplayOverride", "baseMeshDisplayOverride", MFnNumericData::kBoolean);
	nAttr.setStorable(true);
	nAttr.setDefault(true);
	nAttr.setKeyable(true);
	nAttr.setChannelBox(true);
	addAttribute(primitiveGenerator::aDisableBaseMeshOverride);

	// Jiggle
	//


	primitiveGenerator::aJiggleEnabled = nAttr.create( "jiggleEnabled", "jiggleEnabled", MFnNumericData::kBoolean );
	nAttr.setStorable( true );
	nAttr.setDefault( false );
	nAttr.setKeyable( true );
	nAttr.setChannelBox( true );
	nAttr.setHidden(false);
	addAttribute( primitiveGenerator::aJiggleEnabled );

	primitiveGenerator::aTime = uAttr.create("time", "time", MFnUnitAttribute::kTime, 0.0);
	uAttr.setWritable(true);
	uAttr.setReadable(false);
	addAttribute(primitiveGenerator::aTime);

	primitiveGenerator::aStartFrame = nAttr.create("startFrame", "startFrame", MFnNumericData::kInt, 1.0, &status);
	nAttr.setKeyable(true);
	nAttr.setDefault(1.0);
	nAttr.setStorable( true );
	addAttribute(primitiveGenerator::aStartFrame);


	primitiveGenerator::aJiggleAmount = nAttr.create("jiggle", "jiggle", MFnNumericData::kFloat, 0.0);
	nAttr.setKeyable(true);
	nAttr.setDefault(1.0);
	nAttr.setMin(0.0);
	nAttr.setSoftMax(1.0);
	addAttribute(primitiveGenerator::aJiggleAmount);

	primitiveGenerator::aStiffness = nAttr.create("stiffness", "stiffness", MFnNumericData::kFloat, 1.0);
	nAttr.setKeyable(true);
	nAttr.setDefault(0.4);
	nAttr.setMin(0.0);
	nAttr.setMax(1.0);
	addAttribute(primitiveGenerator::aStiffness);

	primitiveGenerator::aDamping = nAttr.create("damping", "damping", MFnNumericData::kFloat, 1.0);
	nAttr.setKeyable(true);
	nAttr.setDefault(0.1);
	nAttr.setMin(0.0);
	nAttr.setMax(1.0);
	addAttribute(primitiveGenerator::aDamping);


	attributeAffects(primitiveGenerator::aTime, primitiveGenerator::aOutMesh);

	attributeAffects(primitiveGenerator::aStartFrame, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aJiggleEnabled, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aJiggleAmount, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aStiffness, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aDamping, primitiveGenerator::aOutMesh);


	attributeAffects(primitiveGenerator::aInCurve, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aRefCurve, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aRadius, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aWidth, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aHeight, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aRotate, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aTwist, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aSides, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aSegmentsLoop, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aSegments, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aUseInputCurve, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aSmoothNormals, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aCapTop, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aAlingToUpVector, primitiveGenerator::aOutMesh);

	attributeAffects(primitiveGenerator::aOnlyKnotSegmentsRes, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aAutoSegments, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aAutoSegmentsRes, primitiveGenerator::aOutMesh);


	attributeAffects(primitiveGenerator::aInLocAPos, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aInLocBPos, primitiveGenerator::aOutMesh);

	attributeAffects(primitiveGenerator::aFirstUpVec, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aCurveZOffset, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aSegmentRamp, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aProfilePresets, primitiveGenerator::aOutMesh);

	// Override
	attributeAffects(primitiveGenerator::aDisableBaseMeshOverride, primitiveGenerator::aOutMesh);

	// UV

	attributeAffects(primitiveGenerator::aCapUVSize, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aUOffset, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aVOffset, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aUOffsetCap, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aVOffsetCap, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aUWidth, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aVWidth, primitiveGenerator::aOutMesh);
	attributeAffects(primitiveGenerator::aUVRotate, primitiveGenerator::aOutMesh);

	//MGlobal::executeCommand( "makePaintable -attrType multiFloat -sm deformer deltaRelax weights;" );

	return MStatus::kSuccess;
}