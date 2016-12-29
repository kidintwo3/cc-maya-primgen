import maya.cmds as cmds

sel_obj = cmds.ls(sl=True, l=True)

if sel_obj:

    child_curves = cmds.listRelatives(sel_obj, ad=True, typ="nurbsCurve", f=True)
    
    if child_curves:
        
        for i in child_curves:
            
            if cmds.attributeQuery( 'worldSpace', node=i, exists=True ):
            
                conn = cmds.listConnections( i + '.worldSpace', s=False, d=True )
                if conn:
                    if cmds.objectType( conn[0], isType='primitiveGenerator' ):
                        cmds.setAttr(conn[0] + '.baseMeshDisplayOverride', False)
    else:
        cmds.warning("No curves selected...")
        
else:
    cmds.warning("Nothing selected...")