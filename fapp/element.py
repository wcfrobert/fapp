import numpy as np
import math

class Element():
    """
    Element object definition.
        tag:                    unique tag number defined by user
        i_node:                 beginning node object
        j_node:                 end node object
        A:                      cross-sectional area
        Ayy:                    shear area along y-axis
        Azz:                    shear area along z-axis
        Iy:                     moment of inertia about minor axis
        Iz:                     moment of inertia about major axis
        J:                      torsion constant
        E:                      Young's modulus
        v:                      Poisson's ratio
        G:                      shear modulus
        beta:                   web rotation. Default value is 0. Input in degrees
        
        length:                 length of member
        dof:                    1x12 list containing the element's DOF
        load:                   1x3 list containing uniform load in local coordinate {w}
        
        T:                      12x12 coordinate transformation matrix [T]
        k_local:                12x12 local stiffness matrix [Ke]
        k_global:               12x12 global stiffness matrix [T]'[Ke][T]
        FEF_local:              12x1 local fixed-end force vector {FEF}
        FEF_global:             12x1 global fixed-end force vector [T]'{FEF}
        f_local:                12x1 local element force vector {F}
        f_global:               12x1 global element force vector [T]'{F}
        d_local:                12x1 local element displacement vector  [T]{u}
        d_global:               12x1 global element displacement vector {u}
    """
    def __init__(self, tag, i_node, j_node, A, Ayy, Azz, Iy, Iz, J, E, v, beta):
        self.tag = tag
        self.i_node = i_node
        self.j_node = j_node
        self.A = A
        self.Ayy = Ayy
        self.Azz = Azz
        self.Iy = Iy
        self.Iz = Iz
        self.J = J
        self.E = E
        self.v = v
        self.G = E/(2*(1+v))
        self.beta = beta
        self.load = [0,0,0]

        self.length = self.compute_length()
        self.dof = self.get_dof()
        self.T = self.compute_T()
        self.k_local, self.k_global = self.compute_K()
        self.FEF_local, self.FEF_global = self.compute_FEF()
        
        self.f_local = None
        self.f_global = None
        self.d_local = None
        self.d_global = None

    def compute_length(self):
        """return the element's length by querying its end node's x,y,z position"""
        i_coord = self.i_node.coord
        j_coord = self.j_node.coord
        length = np.linalg.norm(j_coord - i_coord)
        return length
    
    def get_dof(self):
        """returns a 1x12 list associated with the element's 12 DOF"""
        element_dof = self.i_node.dof + self.j_node.dof
        return element_dof

    def compute_T(self):
        """
        Approach described in MASTAN textbook. 
            1.) compute x direction cosine {gamma_x}. It points from i_node to j_node
            2.) cross {gamma_x} with {0,1,0} then normalize to get {gamma_z}
                * elements that are vertical are special. {gamma_z} points to {0,0,1} or {0,0,-1}
                * whichever respects the right-hand rule and ensures {gamma_z} is always {-1,0,0}
            3.) cross {gamma_z} with {gamma_x} to get web direction vector {gamma_y}
        
        By default, major axis (web direction vector) always points skyward. For vertical elements 
        like columns, major axis always points to -X. This default setting works well. Beams always 
        have major axis resisting gravity load, and columns will always have major axis properties 
        when modeling frames in the X-Y plane.
        """
        # Direction cosine in x-direction point from i_node to j_node
        i_coord = self.i_node.coord
        j_coord = self.j_node.coord
        position_vector = j_coord - i_coord
        cosinex = np.array([
            position_vector[0]/self.length,
            position_vector[1]/self.length,
            position_vector[2]/self.length])
        
        # Direction cosine in z
        if np.all(cosinex == [0,1,0]):
            cosinez = np.array([0,0,1])
        elif np.all(cosinex == [0,-1,0]):
            cosinez = np.array([0,0,-1])
        else:
            cosinez = np.cross(cosinex,np.array([0,1,0]))
            cosinez = cosinez/np.linalg.norm(cosinez)
            
        # Direction cosine in y (web-direction)
        cosiney = np.cross(cosinez,cosinex)
        cosiney = cosiney/np.linalg.norm(cosiney)
        
        # Beta angle rotation
        beta_rad = self.beta*math.pi/180
        gamma_b = np.array([
            [1,0,0],
            [0, math.cos(beta_rad), math.sin(beta_rad)],
            [0, -math.sin(beta_rad), math.cos(beta_rad)]])
        
        # Assemble gamma matrix
        gamma_p = np.zeros([3,3])
        gamma_p[0,:] = cosinex
        gamma_p[1,:] = cosiney
        gamma_p[2,:] = cosinez
        gamma = gamma_b@gamma_p
        
        # Assemble transformation matrix
        zero = np.zeros([3,3])
        T = np.block([
            [gamma,zero,zero,zero],
            [zero,gamma,zero,zero],
            [zero,zero,gamma,zero],
            [zero,zero,zero,gamma]])
        return T
    
    def compute_K(self):
        """
        Assemble element 12x12 stiffness matrix.
        
        Shear deformation is included. Format of matrix follows the CEE 421L 
        course note by Dr. Henri Gavin. Theta_y and Theta_z defined here:
        Link: http://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
        """
        # assign to local variable so I'm not writing "self" a million times
        A = self.A
        Ayy = self.Ayy
        Azz = self.Azz
        Iz = self.Iz
        Iy = self.Iy
        L = self.length
        J = self.J
        G = self.G
        E = self.E
        theta_y = 12*E*Iz/G/Ayy/L/L 
        theta_z = 12*E*Iy/G/Azz/L/L
        
        # decoupled stiffness matrices
        k_axial = np.array([
            [E*A/L,-E*A/L],
            [-E*A/L,E*A/L]])
        k_torsion = np.array([
            [G*J/L,-G*J/L],
            [-G*J/L,G*J/L]])
        k_flexureZ = np.array([
            [12*E*Iz/L**3/(1+theta_y) , 6*E*Iz/L**2/(1+theta_y)       , -12*E*Iz/L**3/(1+theta_y), 6*E*Iz/L**2/(1+theta_y)],
            [6*E*Iz/L**2/(1+theta_y)  , (4+theta_y)*E*Iz/L/(1+theta_y), -6*E*Iz/L**2/(1+theta_y) , (2-theta_y)*E*Iz/L/(1+theta_y)],
            [-12*E*Iz/L**3/(1+theta_y), -6*E*Iz/L**2/(1+theta_y)      , 12*E*Iz/L**3/(1+theta_y) , -6*E*Iz/L**2/(1+theta_y)],
            [6*E*Iz/L**2/(1+theta_y)  , (2-theta_y)*E*Iz/L/(1+theta_y), -6*E*Iz/L**2/(1+theta_y) , (4+theta_y)*E*Iz/L/(1+theta_y)]])
        k_flexureY = np.array([
            [12*E*Iy/L**3/(1+theta_z) , -6*E*Iy/L**2/(1+theta_z)      , -12*E*Iy/L**3/(1+theta_z), -6*E*Iy/L**2/(1+theta_z)],
            [-6*E*Iy/L**2/(1+theta_z) , (4+theta_z)*E*Iy/L/(1+theta_z), 6*E*Iy/L**2/(1+theta_z)  , (2-theta_z)*E*Iy/L/(1+theta_z)],
            [-12*E*Iy/L**3/(1+theta_z), 6*E*Iy/L**2/(1+theta_z)       , 12*E*Iy/L**3/(1+theta_z) , 6*E*Iy/L**2/(1+theta_z)],
            [-6*E*Iy/L**2/(1+theta_z) , (2-theta_z)*E*Iy/L/(1+theta_z), 6*E*Iy/L**2/(1+theta_z)  , (4+theta_z)*E*Iy/L/(1+theta_z)]])

        # Combine into global stiffness matrix
        k_local = np.zeros([12,12])
        k_local[np.ix_([0,6],[0,6])] = k_axial
        k_local[np.ix_([3,9],[3,9])] = k_torsion
        k_local[np.ix_([1,5,7,11],[1,5,7,11])] = k_flexureZ
        k_local[np.ix_([2,4,8,10],[2,4,8,10])] = k_flexureY
        
        # Convert to global coordinate
        k_global = self.T.transpose()@k_local@self.T
        return k_local,k_global
    
    def compute_FEF(self):
        """assemble fixed-end force vector"""
        L = self.length
        w = self.load
        FEF_local = np.array([
            -w[0]*L/2,
            -w[1]*L/2,
            -w[2]*L/2,
            0,
            w[2]*L**2/12,
            -w[1]*L**2/12,
            -w[0]*L/2,
            -w[1]*L/2,
            -w[2]*L/2,
            0,
            -w[2]*L**2/12,
            w[1]*L**2/12,])
        
        # Convert to global coordinate
        FEF_global = self.T.transpose() @ FEF_local
        return FEF_local,FEF_global
    
    def force_recovery(self, d_global):
        """given end node displacement, recover internal forces"""
        # save displacement
        self.d_global = d_global
        self.d_local = self.T @ d_global
        
        # save force
        self.f_local = self.k_local@self.d_local + self.FEF_local
        self.f_global = self.T.transpose() @ self.f_local
        
    
