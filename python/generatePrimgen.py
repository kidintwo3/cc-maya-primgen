import maya.cmds as cmds
import random

def setHairCurve(tr_node):
    if tr_node:
        cmds.select(tr_node, r=True)
        primgen_node = cmds.primitiveGeneratorCommand()

        if primgen_node:
            setPrimgen(primgen_node)

def setPrimgen(node):
    if node:
        if len(node) == 2:

            # Basics
            cmds.setAttr(node[0] + ".radius", 1.2)
            cmds.setAttr(node[0] + ".sides", 2)
            cmds.setAttr(node[0] + ".segments", 23)
            cmds.setAttr(node[0] + ".strands", 1)
            cmds.setAttr(node[0] + ".strandOffset", 0)
			cmds.setAttr(node[0] + ".capTop", 0)

            # Rotation
            cmds.setAttr(node[0] + ".rotate", random.randint(120,180))
            cmds.setAttr(node[0] + ".rotate", 160)
            cmds.setAttr(node[0] + ".rotate", 160)
            cmds.setAttr(node[0] + ".twist", 2)

            # Ramp
            cmds.setAttr(node[0] + ".segmentsRamp[2].segmentsRamp_Position", 1)
            cmds.setAttr(node[0] + ".segmentsRamp[2].segmentsRamp_FloatValue", 1)

            # UV
            cmds.setAttr(node[0] + ".uOffset", 1.0)
            cmds.setAttr(node[0] + ".vOffset", 1.0)

            cmds.setAttr(node[0] + ".uWidth", 0.3)
            cmds.setAttr(node[0] + ".vWidth", 0.68)
			
            # Rotation
            cmds.setAttr(node[0] + ".uvRotate", 180)
            cmds.sets(node[1], edit=True, forceElement="SG_hair")

def generatePrimGen():
    curr_sel = cmds.ls(sl=True)

    if not curr_sel:
        cmds.warning("[PrimGen] No selected groups...")
        return

    sel_grp = cmds.listRelatives( allDescendents=True, type="transform" )

    for i in sel_grp:
        shapes = cmds.listRelatives(i)

        if shapes:
            if cmds.nodeType(shapes[0]) == "nurbsCurve":
                setHairCurve(i)


    cmds.select(curr_sel, r=True)


generatePrimGen()
