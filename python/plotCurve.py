
import maya.cmds as cmds
import maya.OpenMaya as OpenMaya

def getMDagPath(nodeName):
    selList = OpenMaya.MSelectionList()
    selList.add(nodeName)
    mDagPath = OpenMaya.MDagPath()
    selList.getDagPath(0, mDagPath)
    return mDagPath

def getSelection():
    selList = OpenMaya.MSelectionList()
    OpenMaya.MGlobal.getActiveSelectionList(selList)
    sel = []
    selList.getSelectionStrings(sel)
    return sel

curveFn = OpenMaya.MFnNurbsCurve(getMDagPath(getSelection()[0]))
array = OpenMaya.MPointArray()
curveFn.getCVs(array,OpenMaya.MSpace.kObject)

m_locPointsNum = 0

for i in xrange(array.length()):
    print "{ " + str(array[i].x) + "f ," + str(array[i].y) + "f ," + str(array[i].z) + "f },"
    m_locPointsNum += 1

print m_locPointsNum