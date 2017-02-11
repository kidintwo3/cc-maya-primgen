import maya.cmds as cmds
import maya.api.OpenMaya as om

target_obj = "trish_head_m_high_geo"
source_obj = om.MGlobal.getActiveSelectionList()

target_obj_sl = om.MSelectionList()
target_obj_sl.add(target_obj)
target_path = target_obj_sl.getDagPath(0)

up_vec = None

for i in xrange(0, source_obj.length()):

    source_path = source_obj.getDagPath(i)
    
    if source_path.apiType() == 110:
    
        source_path_shape = source_path
        source_path_shape.extendToShape()
        
        if source_path_shape.apiType() == 267:

            mFn_mesh = om.MFnMesh(target_path)
            mFn_curve = om.MFnNurbsCurve(source_path)

            root_cv_p = mFn_curve.cvPosition(0)
            
            if root_cv_p:
                closest_n = mFn_mesh.getClosestNormal(root_cv_p)
                
                if closest_n:
                    up_vec = om.MVector(closest_n[0])
                        
            if up_vec:
                
                conn =  cmds.listConnections( source_path_shape.partialPathName() + '.worldSpace', d=True, s=False )
                
                if conn:
                    cmds.setAttr(conn[0] + '.firstUpVecX', up_vec[0])
                    cmds.setAttr(conn[0] + '.firstUpVecY', up_vec[1])
                    cmds.setAttr(conn[0] + '.firstUpVecZ', up_vec[2])
                        
