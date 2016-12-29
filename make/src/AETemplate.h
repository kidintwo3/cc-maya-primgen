//
//  AETemplates.h
//  primGen
//
//  Created by Hunyadi J�nos on 2015. 01. 03..
//  Copyright (c) 2015. Janos Hunyadi. All rights reserved.
//

#ifndef MELSCRIPTS_H
#define MELSCRIPTS_H

#include <maya/MString.h>


MString mel_AETemplate()
{
	MString s_aeTemplate = MString() + "//deleteUI AttrEdprimitiveGeneratorFormLayout;\r\n"
		"global proc AEprimitiveGeneratorTemplate( string $nodeName )\r\n"
		"{\r\n"
		"    editorTemplate -beginScrollLayout;\r\n"
		"    editorTemplate -beginLayout \"Profile\" -collapse 0;\r\n"
		"        editorTemplate -addControl \"profilePresets\";\r\n"
		"        editorTemplate -addControl \"radius\";\r\n"
		"        editorTemplate -addControl \"width\";\r\n"
		"        editorTemplate -addControl \"height\";\r\n"
		"		editorTemplate -label \"Custom Sides\" -addControl \"sides\";\r\n"
		"		editorTemplate -addSeparator;\r\n"
		"		editorTemplate -addControl \"alingToUpVector\";\r\n"
		"		editorTemplate -addControl \"firstUpVec\";\r\n"
		"    editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Segments\" -collapse 0;\r\n"
		"		editorTemplate -label \"Segments\" -addControl \"segments\";\r\n"
		"		editorTemplate -addControl \"autoSegmentsMultiplier\";\r\n"
		"		editorTemplate -addSeparator;\r\n"
		"		editorTemplate -label \"Auto segments\" -addControl \"autoSegments\";\r\n"
		"		editorTemplate -label \"Segments at knots only\" -addControl \"autoSegmentsKnotsOnly\";\r\n"
		"    editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Strands\" -collapse 0;\r\n"
		"		editorTemplate -label \"Strands\" -addControl \"strands\";\r\n"
		"		editorTemplate -label \"Strand Offset\" -addControl \"strandOffset\";\r\n"
		"		editorTemplate -label \"Strand Thinning\" -addControl \"strandThinning\";\r\n"
		"		editorTemplate -label \"Strand Curl\" -addControl \"strandCurl\";\r\n"
		"		editorTemplate -label \"Strand Curl Wave\" -addControl \"strandCurlWave\";\r\n"
		"		AEaddRampControl( $nodeName + \".strandOffsetRamp\" );\r\n"
		"    editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Segments Translation\" -collapse 0;\r\n"
		"		editorTemplate -label \"Rotation\" -addControl \"rotate\";\r\n"
		"		editorTemplate -label \"Twist\" -addControl \"twist\";\r\n"
		"		editorTemplate -label \"Offset curve normal\" -addControl \"curveZOffset\";\r\n"
		"    editorTemplate -endLayout;\r\n"
		"    \r\n"
		"    editorTemplate -beginLayout \"Ramp Deformer Attributes\" -collapse 0;\r\n"
		"            AEaddRampControl( $nodeName + \".segmentsRamp\" );\r\n"
		"    editorTemplate -endLayout;\r\n"
		"    \r\n"
		"    editorTemplate -beginLayout \"UV settings\" -collapse 0;\r\n"
		"        editorTemplate -addControl \"capUvSize\";\r\n"
		"        editorTemplate -addControl \"uOffset\";\r\n"
		"        editorTemplate -addControl \"vOffset\";\r\n"
		"        editorTemplate -addControl \"uOffsetCap\";\r\n"
		"        editorTemplate -addControl \"vOffsetCap\";\r\n"
		"        editorTemplate -addControl \"uWidth\";\r\n"
		"        editorTemplate -addControl \"vWidth\";\r\n"
		"        editorTemplate -addControl \"uvRotate\";\r\n"
		"    editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Custom profile Curve\" -collapse 0;\r\n"
		"		editorTemplate -callCustom \"AE_referenceCurve_create\" \"AE_referenceCurve_edit\" \"\"; \r\n"
		"	editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Muscle Jiggle / A-B only\" -collapse 0;\r\n"
		"		editorTemplate -label \"Enabled\" -addControl \"jiggleEnabled\";\r\n"
		"		editorTemplate -addSeparator;\r\n"
		"		editorTemplate -label \"Jiggle Amount\" -addControl \"jiggle\";\r\n"
		"		editorTemplate -addControl \"stiffness\";\r\n"
		"		editorTemplate -addControl \"damping\";\r\n"
		"		editorTemplate -addControl \"startFrame\";\r\n"
		"	editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Global switches\" -collapse 0;\r\n"
		"        editorTemplate -addControl \"smoothNormals\";\r\n"
		"        editorTemplate -addControl \"capTop\";\r\n"
		"		editorTemplate -addControl \"baseMeshDisplayOverride\";\r\n"
		"    editorTemplate -endLayout;\r\n"
		"    \r\n"
		"	editorTemplate -beginLayout \"Bake Mesh\" -collapse 1;\r\n"
		"		editorTemplate -callCustom \"AE_bakeMesh_create\" \"AE_bakeMesh_edit\" \"\"; \r\n"
		"	editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Extra options\" -collapse 0;\r\n"
		"		editorTemplate -callCustom \"AE_selectMesh_create\" \"AE_selectMesh_edit\" \"\"; \r\n"
		"	editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	editorTemplate -beginLayout \"Plug-in Info\" -collapse 1;\r\n"
		"		editorTemplate -callCustom \"AE_primgen_website_create\" \"AE_primgen_website_edit\" \"\";\r\n"
		"	editorTemplate -endLayout;\r\n"
		"	\r\n"
		"	AEabstractBaseCreateTemplate $nodeName;\r\n"
		"	\r\n"
		"    editorTemplate -addExtraControls;\r\n"
		"	editorTemplate -endScrollLayout;\r\n"
		"	editorTemplate -suppress \"outMesh\";\r\n"
		"	editorTemplate -suppress \"inCurve\";\r\n"
		"	\r\n"
		"	editorTemplate -suppress \"locatorAPos\";\r\n"
		"	editorTemplate -suppress \"locatorBPos\";\r\n"
		"	editorTemplate -suppress \"firstUpVecX\";\r\n"
		"	editorTemplate -suppress \"firstUpVecY\";\r\n"
		"	editorTemplate -suppress \"firstUpVecZ\";\r\n"
		"	// Legacy\r\n"
		"	editorTemplate -suppress \"segmentsLoop\";\r\n"
		"}\r\n"
		"// ----------------------------\r\n"
		"global proc AE_primgen_launch_website()\r\n"
		"{\r\n"
		"    launch -web \"http://gumroad.com/creativecase\";\r\n"
		"}\r\n"
		"global proc AE_primgen_website_create(string $attrName)\r\n"
		"{\r\n"
		"	string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    rowLayout -numberOfColumns 3 -adjustableColumn 2 -bgc 0.2 0.2 0.2;\r\n"
		"    text -al \"left\" -label \"Support / Info\";\r\n"
		"    textField -ed false;\r\n"
		"    iconTextButton -al \"left\" -ann \"Gumroad Page\" -style \"iconOnly\" -image1 \"primitiveGenerator_CCLogo.png\" -c \"AE_primgen_launch_website()\";\r\n"
		"    setParent ..;\r\n"
		"}\r\n"
		"global proc AE_primgen_website_edit(string $attrName)\r\n"
		"{\r\n"
		"}\r\n"
		"global proc AE_referenceCurve_create(string $attrName)\r\n"
		"{\r\n"
		"	string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    string $sm_pr_inputCurve[] = eval(\"listConnections \" + $nodeName[0] + \".refCurve\");\r\n"
		"    rowLayout -numberOfColumns 4 -adjustableColumn 2 -bgc 0.2 0.2 0.2;\r\n"
		"    text -al \"left\" -label \"Profile Curve\";\r\n"
		"    textField -ed false -bgc 0.5 0.5 0.5 -tx $sm_pr_inputCurve[0] \"t_refName\";\r\n"
		"    button -label \"Set\" -h 20 -w 40 -bgc 0.2 0.8 0.6 -c  (\"AE_referenceCurve_add \" + $nodeName[0] ) \"pres_refcurve_set\";\r\n"
		"    button -label \"Remove\" -h 20 -w 50 -bgc 0.8 0.2 0.3 -c  (\"AE_referenceCurve_remove \" + $nodeName[0] ) \"pres_refcurve_remove\";\r\n"
		"    setParent ..;\r\n"
		"}\r\n"
		"global proc AE_referenceCurve_edit(string $attrName)\r\n"
		"{\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"	button -edit -c (\"AE_referenceCurve_add \" + $nodeName[0] ) \"pres_refcurve_set\";\r\n"
		"	button -edit -c (\"AE_referenceCurve_remove \" + $nodeName[0] ) \"pres_refcurve_remove\";\r\n"
		"	string $sm_pr_inputCurve[] = eval(\"listConnections \" + $nodeName[0] + \".refCurve\");\r\n"
		"	textField -edit -tx $sm_pr_inputCurve[0] \"t_refName\";\r\n"
		"}\r\n"
		"global proc AE_selectMesh_create(string $attrName)\r\n"
		"{\r\n"
		"	string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    string $sm_pr_inputCurve[] = eval(\"listConnections \" + $nodeName[0] + \".refCurve\");\r\n"
		"    rowLayout -numberOfColumns 3  -bgc 0.2 0.2 0.2;\r\n"
		"    button -label \"Select input curve / output mesh\" -h 20 -w 200 -bgc 0.8 0.4 0.4 -c  (\"AE_select_curve_mesh \" + $nodeName[0] ) \"select_curve_mesh\";\r\n"
		"    button -label \"Select output mesh\" -h 20 -w 100 -bgc 0.6 0.2 0.8 -c  (\"AE_select_mesh \" + $nodeName[0] ) \"select_mesh\";\r\n"
		"    button -label \"Select input curve\" -h 20 -w 100 -bgc 0.2 0.6 0.8 -c  (\"AE_select_curve \" + $nodeName[0] ) \"select_curve\";\r\n"
		"    setParent ..;\r\n"
		"    \r\n"
		"    \r\n"
		"}\r\n"
		"global proc AE_selectMesh_edit(string $attrName)\r\n"
		"{\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"	button -edit -c (\"AE_select_mesh \" + $nodeName[0] ) \"select_mesh\";\r\n"
		"	button -edit -c (\"AE_select_curve \" + $nodeName[0] ) \"select_curve\";\r\n"
		"	button -edit -c (\"AE_select_curve_mesh \" + $nodeName[0] ) \"select_curve_mesh\";\r\n"
		"}\r\n"
		"global proc AE_select_mesh(string $attrName)\r\n"
		"{\r\n"
		"	\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    string $pg_list_otputMeshes[] = eval(\"listConnections \" + $nodeName[0] + \".outMesh\");\r\n"
		"    if (size($pg_list_otputMeshes) > 0) \r\n"
		"    {\r\n"
		"        select -r $pg_list_otputMeshes[0];\r\n"
		"    }\r\n"
		"}\r\n"
		"global proc AE_select_curve(string $attrName)\r\n"
		"{\r\n"
		"	\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    string $pg_list_inputCurves[] = eval(\"listConnections \" + $nodeName[0] + \".inCurve\");\r\n"
		"    if (size($pg_list_inputCurves) > 0) \r\n"
		"    {\r\n"
		"        select -r $pg_list_inputCurves[0];\r\n"
		"    }\r\n"
		"}\r\n"
		"global proc AE_select_curve_mesh(string $attrName)\r\n"
		"{\r\n"
		"	\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    string $pg_list_inputCurves[] = eval(\"listConnections \" + $nodeName[0] + \".inCurve\");\r\n"
		"    string $pg_list_otputMeshes[] = eval(\"listConnections \" + $nodeName[0] + \".outMesh\");\r\n"
		"    if (size($pg_list_inputCurves) > 0) \r\n"
		"    {\r\n"
		"        if (size($pg_list_otputMeshes) > 0)\r\n"
		"        {\r\n"
		"            select -r $pg_list_inputCurves[0];\r\n"
		"            select -add $pg_list_otputMeshes[0];\r\n"
		"        }\r\n"
		"    }\r\n"
		"}\r\n"
		"global proc AE_bakeMesh_create(string $attrName)\r\n"
		"{\r\n"
		"	string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    string $sm_pr_inputCurve[] = eval(\"listConnections \" + $nodeName[0] + \".refCurve\");\r\n"
		"    rowLayout -numberOfColumns 2 -adjustableColumn 2 -bgc 0.2 0.2 0.2;\r\n"
		"    button -label \"Bake mesh / delete node\" -h 20 -w 200 -bgc 0.8 0.2 0.3 -c  (\"AE_bake_and_delete_mesh \" + $nodeName[0] ) \"bake_mesh_and_delete\";\r\n"
		"    button -label \"Bake mesh\" -h 20 -w 200 -bgc 0.2 0.8 0.6 -c  (\"AE_bake_mesh \" + $nodeName[0] ) \"bake_mesh\";\r\n"
		"    setParent ..;\r\n"
		"}\r\n"
		"global proc AE_bakeMesh_edit(string $attrName)\r\n"
		"{\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"	button -edit -c (\"AE_bake_and_delete_mesh \" + $nodeName[0] ) \"bake_mesh_and_delete\";\r\n"
		"	button -edit -c (\"AE_bake_mesh \" + $nodeName[0] ) \"bake_mesh\";\r\n"
		"}\r\n"
		"global proc AE_bake_and_delete_mesh(string $attrName)\r\n"
		"{\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    \r\n"
		"    string $pg_conn_outputMeshs[] = `listConnections -d on -s off ($nodeName[0] + \".outMesh\")`;\r\n"
		"    if (size($pg_conn_outputMeshs) > 0) \r\n"
		"    {\r\n"
		"        undoInfo -ock;\r\n"
		"        setAttr ($nodeName[0] + \".baseMeshDisplayOverride\") 0;\r\n"
		"        delete $nodeName[0];\r\n"
		"        select -r $pg_conn_outputMeshs[0]; \r\n"
		"        undoInfo -cck;\r\n"
		"    }\r\n"
		"}\r\n"
		"global proc AE_bake_mesh(string $attrName)\r\n"
		"{\r\n"
		"    string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    string $pg_list_inputCurves[] = eval(\"listConnections \" + $nodeName[0] + \".outMesh\");\r\n"
		"    if (size($pg_list_inputCurves) > 0) \r\n"
		"    {\r\n"
		"        string $new_node_Name[] = `duplicate $pg_list_inputCurves[0]`;\r\n"
		"        \r\n"
		"        if (size($new_node_Name) > 0)\r\n"
		"        {\r\n"
		"            string $new_node_shape[] = `listRelatives -s -path $new_node_Name[0]`;\r\n"
		"            undoInfo -ock;\r\n"
		"            setAttr ($new_node_shape[0] + \".overrideEnabled\") 0;\r\n"
		"            select -r $new_node_Name[0];\r\n"
		"            undoInfo -cck;\r\n"
		"        }\r\n"
		"    }\r\n"
		"}\r\n"
		"global proc AE_referenceCurve_add(string $attrName)\r\n"
		"{\r\n"
		"	string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"	\r\n"
		"    string $sel[] = `ls -sl`;\r\n"
		"    if (size($sel) > 0) \r\n"
		"    {\r\n"
		"        string $selShape[] = `listRelatives -shapes $sel[0]`;\r\n"
		"        if (size($selShape) > 0) \r\n"
		"        {\r\n"
		"            string $nodeType;\r\n"
		"            $nodeType = `nodeType $selShape[0]`;\r\n"
		"            if( $nodeType == \"nurbsCurve\")\r\n"
		"            {\r\n"
		"                connectAttr -force ($selShape[0]+\".local\") ($nodeName[0]+\".refCurve\") ;\r\n"
		"                textField -edit -bgc 0.5 0.5 0.5 -tx $selShape[0] \"t_refName\";\r\n"
		"            }\r\n"
		"        }\r\n"
		"    }\r\n"
		"	\r\n"
		"    \r\n"
		"}\r\n"
		"global proc AE_referenceCurve_remove(string $attrName)\r\n"
		"{\r\n"
		"	string $nodeName[];\r\n"
		"    tokenize($attrName, \".\", $nodeName);\r\n"
		"    \r\n"
		"    string $pg_list_inputCurves[] = eval(\"listConnections \" + $nodeName[0] + \".refCurve\");\r\n"
		"    if (size($pg_list_inputCurves) > 0) \r\n"
		"    {\r\n"
		"        eval(\"disconnectAttr \" + $pg_list_inputCurves[0] + \".local\" + \" \" + $nodeName[0] + \".refCurve\");\r\n"
		"        textField -edit -bgc 0.5 0.5 0.5 -tx \"\"  \"t_refName\";\r\n"
		"    }\r\n"
		"}\r\n";

	return s_aeTemplate;

}

MString mel_createShelf()
{

	MString s_aeTemplate = MString() + "int $cc_doesShelfExist = `shelfLayout -query -ex \"CreativeCase\"`;"
		"if ($cc_doesShelfExist == 1)"
		"{"
		"	string $shelfButtons[] = `shelfLayout -q -ca \"CreativeCase\"`;"
		"	int $ex_b01,$ex_b02 = 0;"
		"	for( $i=0; $i<size($shelfButtons); ++$i )"
		"	{"
		"		if( `control -exists $shelfButtons[$i]` == true)"
		"		{"
		"			if (`control -q -docTag $shelfButtons[$i]` == \"pg_createPgButton\") {$ex_b01 = 1;}"
		"			if (`control -q -docTag $shelfButtons[$i]` == \"pg_addMuscleButton\") {$ex_b02 = 1;}"
		"		}"
		"	}"
		"	if ($ex_b01 == 0) {shelfButton -p \"CreativeCase\" -dtg \"pg_createPgButton\" -annotation \"Add a PrimGen modifier on its own or to curves\" -image1 \"primitiveGenerator.png\" -command \"primitiveGeneratorCommand\";}"
		"	if ($ex_b02 == 0) {shelfButton -p \"CreativeCase\" -dtg \"pg_addMuscleButton\" -annotation \"Add a PrimGen muscle modifier to the scene\" -image1 \"primitiveGenerator_muscle.png\" -command \"primitiveGeneratorCommand -muscle 1\";}"
		"}"
		"	"
		"if ($cc_doesShelfExist == false)"
		"{"
		"		shelfLayout -cellWidth 33 -cellHeight 33 -p $gShelfTopLevel CreativeCase;"
		"		shelfButton -p \"CreativeCase\" -dtg \"pg_createPgButton\" -annotation \"Add a PrimGen modifier on its own or to curves\" -image1 \"primitiveGenerator.png\" -command \"primitiveGeneratorCommand\";"
		"		shelfButton -p \"CreativeCase\" -dtg \"pg_addMuscleButton\" -annotation \"Add a PrimGen muscle modifier to the scene\" -image1 \"primitiveGenerator_muscle.png\" -command \"primitiveGeneratorCommand -muscle 1\";"
		"}";


	return s_aeTemplate;
}

#endif