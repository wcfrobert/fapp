import numpy as np

class Node:
    """
    Node object definition.
        tag:               unique tag number defined by user
        coord:             1x3 array containing x,y,z coordinate
        dof:               1x6 list containing dof numbers
        
        load:              1x6 list of nodal loads (Fx,Fy,Fz,Mx,My,Mz)
        fixity:            1x6 list of nodal fixity (ux,uy,uz,rotx,roty,rotz)
        disp:              1x6 array of nodal displacement (ux,uy,uz,rotx,roty,rotz)
        coord_disp:        1x3 array of displaced x,y,z coordinate post analysis
        reaction:          1x6 array of the support reactions (if the node is fixed)
    """
    def __init__(self, tag, x, y, z):
        self.tag = tag
        self.coord = np.array([x, y, z])
        self.dof = self.assign_dof()
        
        self.load = None
        self.fixity = None
        self.disp = None
        self.coord_disp = None
        self.reaction = None
        
    def assign_dof(self):
        """
        Given node tag number, determine DOF numbering.
        For convenience, start DOF at 0 to match 0-based indexing in python.
        
        Node 1 has DOF 0,1,2,3,4,5
        Node 2 has DOF 6,7,8,9,10,11
        and so on...
        """
        start = 6*self.tag
        node_dof = [start,start+1,start+2,start+3,start+4,start+5]
        return node_dof
    
