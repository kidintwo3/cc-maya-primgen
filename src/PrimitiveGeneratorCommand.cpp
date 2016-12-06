#include "PrimitiveGeneratorCommand.h"
#include "DagUtils.h"

primitiveGeneratorCommand::primitiveGeneratorCommand()
{
}


void* primitiveGeneratorCommand::creator()
{
	return new primitiveGeneratorCommand;
}


bool primitiveGeneratorCommand::isUndoable() const
{
	return true;
}

MSyntax primitiveGeneratorCommand::newSyntax()
{
	MSyntax syntax;

	syntax.addFlag( "-r", "-remove", MSyntax::kString );
	syntax.addFlag( "-m", "-muscle", MSyntax::kBoolean );

	syntax.setObjectType( MSyntax::kSelectionList, 0, 1 );
	syntax.useSelectionAsDefault( true );

	syntax.enableEdit( false );
	syntax.enableQuery( false );

	return syntax;
}

MStatus primitiveGeneratorCommand::createPrimGenFromCurves(MDagPathArray p_currSelTrA, MDagPathArray p_currSelShapeA)
{
	MStatus status;



	for (unsigned int i = 0; i < p_currSelTrA.length(); i++)
	{
		p_currSelShapeA[i].extendToShape();

		if (p_currSelShapeA[i].node().apiType() != MFn::kNurbsCurve)
		{
			MGlobal::displayWarning(MString() + "[PrimGen] Only select curves, / " + p_currSelShapeA[i].partialPathName() + " / is not a curve");
			return MStatus::kFailure;
		}

	}

	o_primGenNodeA.clear();
	o_outputMeshA.clear();

	for (unsigned int i = 0; i < p_currSelShapeA.length(); i++)
	{

		MObject o_primGenNode = createNodeCustom(m_DAGMod, "primitiveGenerator" );
		MObject o_outputMesh = createNodeMaya(m_DEPNode, "mesh");

		o_primGenNodeA.append(o_primGenNode);
		o_outputMeshA.append(o_outputMesh);

		connectPlug(m_DAGMod,p_currSelShapeA[i].node(),o_primGenNode,"worldSpace","inCurve");
		connectPlug(m_DAGMod,o_primGenNode,o_outputMesh,"outMesh","inMesh");

		MFnNurbsCurve mFnC(p_currSelShapeA[i]);
		// MGlobal::displayInfo(MString() + mFnC.length() );

		//		int i_curveLength = int(mFnC.length());

		setPlugs(o_primGenNode, "segments", "50");
		setPlugs(o_primGenNode, "profilePresets", "0");
		setPlugs(o_primGenNode, "radius", "0.5");
		setPlugs(o_primGenNode, "sides", "8");

		MFnDependencyNode mfDgN(o_primGenNode);
		MPlug a_curveAttribute = mfDgN.findPlug("segmentsRamp", status);

		MRampAttribute a_Ramp(a_curveAttribute);
		a_Ramp.setValueAtIndex(1.0,0);

		assignInitialShadingGroup(o_outputMesh);


		// Depricated - Now done during runtime with optional plug

		//// Overrides on output mesh

		//MFnDependencyNode mfDgN_out(o_outputMesh);
		//MPlug p_overrideEnabled = mfDgN_out.findPlug("overrideEnabled", false, &status);
		//CHECK_MSTATUS_AND_RETURN_IT(status);
		//MPlug p_overrideDisplayType = mfDgN_out.findPlug("overrideDisplayType", false, &status);
		//CHECK_MSTATUS_AND_RETURN_IT(status);

		//p_overrideEnabled.setBool(true);
		//p_overrideDisplayType.setInt(2);



	}


	return MStatus::kSuccess;
}

MStatus primitiveGeneratorCommand::createPrimGenFromLocators()
{
	MStatus status;

	o_primGenNodeA.clear();
	o_outputMeshA.clear();

	MObject o_primGenNode = createNodeCustom(m_DAGMod, "primitiveGenerator" );
	MObject o_outputMesh = createNodeMaya(m_DEPNode, "mesh");
	o_locA = createNodeMaya(m_DEPNode, "PrimitiveGeneratorLoc" );
	o_locB = createNodeMaya(m_DEPNode, "PrimitiveGeneratorLoc" );

	o_primGenNodeA.append(o_primGenNode);
	o_outputMeshA.append(o_outputMesh);

	connectPlug(m_DAGMod,o_primGenNode,o_outputMesh,"outMesh","inMesh");

	connectPlug(m_DAGMod,o_locA,o_primGenNode,"worldMatrix","locatorAPos");
	connectPlug(m_DAGMod,o_locB,o_primGenNode,"worldMatrix","locatorBPos");


	MFnTransform mfnTrA(o_locA);
	MFnTransform mfnTrB(o_locB);

	mfnTrA.setTranslation(MPoint(-5.0,0.0,0.0), MSpace::kTransform);
	mfnTrB.setTranslation(MPoint(5.0,0.0,0.0), MSpace::kTransform);

	MFnDependencyNode mfDgN(o_primGenNode);
	MPlug a_curveAttribute = mfDgN.findPlug("segmentsRamp", &status);
	CHECK_MSTATUS_AND_RETURN_IT(status);
	MRampAttribute a_Ramp(a_curveAttribute);

	//// Overrides on output mesh

	//MFnDependencyNode mfDgN_out(o_outputMesh);
	//MPlug p_overrideEnabled = mfDgN_out.findPlug("overrideEnabled", false, &status);
	//CHECK_MSTATUS_AND_RETURN_IT(status);
	//MPlug p_overrideDisplayType = mfDgN_out.findPlug("overrideDisplayType", false, &status);
	//CHECK_MSTATUS_AND_RETURN_IT(status);

	//p_overrideEnabled.setBool(true);
	//p_overrideDisplayType.setInt(2);



	if (!m_muscle)
	{

		setPlugs(o_primGenNode, "segments", "20");
		setPlugs(o_primGenNode, "profilePresets", "0");
		setPlugs(o_primGenNode, "radius", "0.5");
		setPlugs(o_primGenNode, "sides", "20");
		setPlugs(o_primGenNode, "jiggleEnabled", "false");

		a_Ramp.setValueAtIndex(1.0,0);

	}




	if (m_muscle)
	{
		setPlugs(o_primGenNode, "jiggleEnabled", "true");
		setPlugs(o_primGenNode, "segments", "20");
		setPlugs(o_primGenNode, "profilePresets", "0");
		setPlugs(o_primGenNode, "radius", "1.0");
		setPlugs(o_primGenNode, "sides", "10");

		MIntArray m_curve_interps;
		MFloatArray	m_curve_positions;
		MFloatArray	m_curve_values;

		m_curve_interps.append(MRampAttribute::kSpline);
		m_curve_interps.append(MRampAttribute::kSpline);
		m_curve_interps.append(MRampAttribute::kSpline);

		m_curve_positions.append(0.0);
		m_curve_positions.append(0.5);
		m_curve_positions.append(1.0);

		m_curve_values.append(0.0);
		m_curve_values.append(1.0);
		m_curve_values.append(0.0);


		status = a_Ramp.setRamp(m_curve_values, m_curve_positions, m_curve_interps);
		CHECK_MSTATUS_AND_RETURN_IT(status);



		// Connect time
		MSelectionList selection;
		status = selection.add( "time1" );
		CHECK_MSTATUS_AND_RETURN_IT(status);
		MObject time1;
		status = selection.getDependNode( 0, time1 );
		CHECK_MSTATUS_AND_RETURN_IT(status);
		MFnDependencyNode fnTime1( time1 );

		MDGModifier dgMod;
		status = dgMod.connect( fnTime1.findPlug( "outTime" ), mfDgN.findPlug( "time" ) );
		CHECK_MSTATUS_AND_RETURN_IT(status);
		status = dgMod.doIt();
		CHECK_MSTATUS_AND_RETURN_IT(status);

	}

	assignInitialShadingGroup(o_outputMesh);

	return MStatus::kSuccess;
}


MStatus primitiveGeneratorCommand::doIt( const MArgList& argList )
{
	MStatus status;

	MArgDatabase argData( syntax(), argList, &status );
	CHECK_MSTATUS_AND_RETURN_IT(status);

	m_muscle = false;

	if ( argData.isFlagSet( "muscle" ) ) { m_muscle = argData.flagArgumentBool("muscle",0, &status); }
	CHECK_MSTATUS_AND_RETURN_IT(status);


	// Get Selected Object


	MDagPathArray p_currSelTrA = getCurrSelectionDAGArray();
	MDagPathArray p_currSelShapeA = p_currSelTrA;

	if (p_currSelTrA.length() == 0)
	{
		MGlobal::displayInfo(MString() + "[PrimGen] Nothing selected. Creating PrimGen with Locators");
		createPrimGenFromLocators();
	}

	else
	{
		status = createPrimGenFromCurves(p_currSelTrA, p_currSelShapeA);
		if (status == MStatus::kFailure)
		{
			MGlobal::displayWarning(MString() + "[PrimGen] Selection is not a curve" );
			return::MStatus::kFailure;
		}
		MGlobal::displayInfo(MString() + "[PrimGen] Attaching PrimGen to curves" );
	}

	if (o_primGenNodeA.length() != 0)
	{

		MStringArray strA_pg;
		MStringArray strA_om;

		for (int i = 0; i < o_primGenNodeA.length(); i++)
		{
			MFnDependencyNode mfDgN(o_primGenNodeA[i]);
			strA_pg.append(mfDgN.name());
		}



		for (int i = 0; i < o_outputMeshA.length(); i++)
		{
			MFnDependencyNode mfDgN(o_outputMeshA[i]);
			strA_om.append(mfDgN.name());
		}


		MPxCommand::clearResult();
		MPxCommand::appendToResult(strA_pg);
		MPxCommand::appendToResult(strA_om);
	}




	return redoIt();
}

MStatus primitiveGeneratorCommand::redoIt()
{
	MStatus status;







	return MS::kSuccess;
}

MStatus primitiveGeneratorCommand::undoIt()
{
	MStatus status;

	// Restore the initial state
	status = m_DGMod.undoIt();
	CHECK_MSTATUS_AND_RETURN_IT( status );

	status = m_DAGMod.undoIt();
	CHECK_MSTATUS_AND_RETURN_IT(status);

	for (unsigned int i = 0; i < o_primGenNodeA.length(); i++)
	{
		deleteNode(o_primGenNodeA[i]);
		deleteNode(o_outputMeshA[i]);
	}



	return MS::kSuccess;
}