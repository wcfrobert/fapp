import math
import numpy as np
import fapp.node
import fapp.element
import scipy.sparse
import scipy.sparse.linalg
import time


class Analysis():
    """
    Analysis object definition.
        node_list                       list storing all nodes in the structure
        element_list                    list storing all elements in the structure      
        
        freeDOF                         list storing free degree of freedoms
        fixedDOF                        list storing fixed degree of freedoms (support)
        dispDOF                         list storing prescribed degree of freedoms (imposed displacement)
        K_structure                     overall structure stiffness matrix
        P_structure                     overall structure external load vector
        FEF_structure                   overall structure fixed-end force vector
        d_structure                     overall structure displacement vector
        
        DEFL                            N_node x 6 matrix containing nodal displacements
        REACT                           N_node x 6 matrix containing reaction forces at fixed nodes
        ELE_FOR                         N_element x 12 matrix containing member-end forces in local coordinate
        ERROR                           back calculated P_structure vector to quantify error
        SOLUTION_FLAG                   0 = failed or not solved yet. 1 = success
        
        linear system of equations: {P} - {FEF} = [K]{d}. The matrix and vectors are partitioned prior to solving.
        subscript f = free DOF, s = fixed DOF, n = prescribed DOF
        
            kff, kfn, kfs, knf, knn, kns, ksf, ksn, kss     partitioned stiffness matrix
            pf, ps, pn                                      partitioned external load vector
            feff, fefs, fefn                                partitioned fixed-end force vector
            df, ds = 0, dn = assigned by user               partitioned displacement vector
            Rf = 0, Rn, Rs                                  partitioned reaction vector
            
        Note that ps != Rs and pn != Rn. We use ps and pn to denote loads applied directly on top of
        support, which is an unusual edge case because the applied load goes directly into the support
        without affecting the overall structure, but this is still possible condition.
    """
    def __init__(self):
        self.node_list=[]
        self.element_list=[]
        
        self.freeDOF = []
        self.fixedDOF = []
        self.dispDOF = []
        self.K_structure = None
        self.kff = None
        self.kfn = None
        self.kfs = None
        self.knf = None
        self.knn = None
        self.kns = None
        self.ksf = None
        self.ksn = None
        self.kss = None
        self.P_structure = None
        self.Pf = None
        self.Ps = None
        self.Pn = None
        self.FEF_structure = None
        self.feff = None
        self.fefs = None
        self.fefn = None
        self.d_structure = None
        self.df = None
        self.dn = None
        self.R_structure = None
        self.Rn = None
        self.Rs = None
        
        self.DEFL = None
        self.REACT = None
        self.ELE_FOR = None
        self.ERROR = None
        self.SOLUTION_FLAG = 0
    
    def add_node(self, node_tag, x, y, z):
        """Factory for node object. Create node object and append to node_list"""
        node_object = fapp.node.Node(node_tag, x, y, z)
        self.node_list.append(node_object)
    
    def add_element(self, ele_tag, i_tag, j_tag,A,Ayy,Azz,Iy,Iz,J,E,v=0.3,beta=0):
        """Factory for element object. Create element object and append to element_list"""
        i_node = self.node_list[i_tag]
        j_node = self.node_list[j_tag]
        element_object = fapp.element.Element(ele_tag,i_node,j_node,A,Ayy,Azz,Iy,Iz,J,E,v,beta)
        self.element_list.append(element_object)
    
    def add_fixity(self, node_tag, ux="nan", uy="nan", uz="nan", rx="nan", ry="nan", rz="nan"):
        """
        Add fixity information to structure.
        0 = fixed, "nan" = free, float = prescribed displacement
        """
        selected_node = self.node_list[node_tag]
        selected_node.fixity = [ux,uy,uz,rx,ry,rz]
    
    def add_load_nodal(self, node_tag, Fx=0, Fy=0, Fz=0, Mx=0, My=0, Mz=0):
        """add nodal load to structure"""
        selected_node = self.node_list[node_tag]
        selected_node.load = [Fx,Fy,Fz,Mx,My,Mz]
    
    def add_load_member(self, ele_tag, wx=0, wy=0, wz=0):
        """add member load to structure"""
        selected_member = self.element_list[ele_tag]
        selected_member.load = [wx,wy,wz]
        selected_member.FEF_local, selected_member.FEF_global = selected_member.compute_FEF()
    
    def solve(self, print_info = True):
        """start analysis routine"""
        if print_info:
            print("Starting analysis routine...")
            print("Total number of nodes: {}".format(len(self.node_list)))
            print("Total number of elements: {}".format(len(self.element_list)))
            print("Total number of DOF: {}".format(6*len(self.node_list)))
        time_start = time.time()
        self.classify_DOF()
        self.assemble_stiffness()
        self.assemble_load()
        self.check_condition()
        self.compute_displacement()
        self.compute_reaction()
        self.force_recovery()
        self.disp_recovery()
        self.compute_error()
        self.SOLUTION_FLAG = 1
        time_end = time.time()
        if print_info:
            print("Analysis completed! Elapsed time: {:.4f} seconds".format(time_end - time_start))

    def classify_DOF(self):
        """classify degree of freedoms into freeDOF, fixedDOF, dispDOF""" 
        # these list of indices are used to partition the stiffness matrix and load vector
        for node in self.node_list:
            if node.fixity == None:
                self.freeDOF = self.freeDOF + node.dof
            else:              
                for j in range(6):   
                    if node.fixity[j] == 0:
                        self.fixedDOF.append(node.dof[j])
                    elif type(node.fixity[j])==int or type(node.fixity[j])==float:
                        self.dispDOF.append(node.dof[j])
                    else:
                        self.freeDOF.append(node.dof[j])
                        
        # indices must be sorted before use
        self.freeDOF = np.sort(self.freeDOF).astype(int)
        self.fixedDOF = np.sort(self.fixedDOF).astype(int)
        self.dispDOF = np.sort(self.dispDOF).astype(int)
    
    def assemble_stiffness(self):
        """assemble and partition sparse global stiffness matrix"""
        # First create matrix in COO sparse format
        COO_i=[]
        COO_j=[]
        COO_v=[]
        for ele in self.element_list:
            non_zero_index = ele.k_global.nonzero() #i,j of all non-zero positions
            for i in range(len(non_zero_index[0])):
                # convert local element index to global structure matrix index (i.e. DOF)
                m = non_zero_index[0][i]
                n = non_zero_index[1][i]
                idx_i = ele.dof[m]
                idx_j = ele.dof[n]
                val = ele.k_global[m,n]
                COO_i.append(idx_i)
                COO_j.append(idx_j)
                COO_v.append(val)
        
        # Convert to CSR sparse format
        self.K_structure = scipy.sparse.csr_matrix( (COO_v,(COO_i,COO_j)))
        
        # Partition stiffness matrix
        f = self.freeDOF
        n = self.dispDOF
        s = self.fixedDOF
        self.kff = self.K_structure[np.ix_(f,f)]
        self.kfn = self.K_structure[np.ix_(f,n)]
        self.kfs = self.K_structure[np.ix_(f,s)]
        self.knf = self.K_structure[np.ix_(n,f)]
        self.knn = self.K_structure[np.ix_(n,n)]
        self.kns = self.K_structure[np.ix_(n,s)]
        self.ksf = self.K_structure[np.ix_(s,f)]
        self.ksn = self.K_structure[np.ix_(s,n)]
        self.kss = self.K_structure[np.ix_(s,s)]

    def assemble_load(self):
        """assemble and partition load vector"""
        # external load vector
        self.P_structure = []
        for node in self.node_list:
            if node.load == None:
                self.P_structure = self.P_structure + [0,0,0,0,0,0]
            else:
                self.P_structure = self.P_structure + node.load
        self.P_structure = np.array(self.P_structure)
        self.Pf = self.P_structure[self.freeDOF]
        self.Ps = self.P_structure[self.fixedDOF]
        self.Pn = self.P_structure[self.dispDOF]
        
        # fixed-end force vector
        self.FEF_structure = np.zeros(len(self.P_structure))
        for ele in self.element_list:
            self.FEF_structure[np.ix_(ele.dof)] = self.FEF_structure[np.ix_(ele.dof)] + ele.FEF_global
        self.feff = self.FEF_structure[self.freeDOF]
        self.fefs = self.FEF_structure[self.fixedDOF]
        self.fefn = self.FEF_structure[self.dispDOF]
    
    def check_condition(self):
        """calculate condition number"""
        cond_number = np.linalg.cond(self.kff.toarray())
        cond_log10 = math.log10(cond_number)
        if cond_log10 > 16:
            print("\n-------WARNING----------")
            print("Log 10 of condition Number = {:.2f}".format(cond_log10))
            print("Stiffness matrix may be unstable or ill-conditioned")
            print("------------------------")
    
    def compute_displacement(self):
        """solve linear equation and determine nodal displacements"""       
        # imposed disp (dn) already known and fixed disp (ds) = 0
        self.dn = np.array([])
        for node in self.node_list:
            if node.fixity != None:
                for fix in node.fixity:
                    if (type(fix) == int or type(fix) == float) and fix != 0:
                        self.dn = np.append(self.dn, fix)

        # calculate free disp (df) via scipy
        self.df = scipy.sparse.linalg.dsolve.spsolve(self.kff,self.Pf - self.feff - self.kfn @ self.dn, use_umfpack=True)
        
        # assemble DEFL matrix
        N_node = len(self.node_list)
        self.d_structure = np.zeros(N_node*6)
        self.d_structure[np.ix_(self.freeDOF)] = self.df
        self.d_structure[np.ix_(self.dispDOF)] = self.dn
        self.DEFL = self.d_structure.reshape(N_node,6)  
    
    def compute_reaction(self):
        """After solving for nodal displacement, compute reaction forces"""
        # calculate reactions
        self.Rn = self.knf @ self.df + self.knn @ self.dn + self.fefn - self.Pn
        self.Rs = self.ksf @ self.df + self.ksn @ self.dn + self.fefs - self.Ps
        
        # assemble REACT matrix
        N_node = len(self.node_list)
        self.R_structure = np.zeros(N_node*6)
        self.R_structure[np.ix_(self.fixedDOF)] = self.Rs
        self.R_structure[np.ix_(self.dispDOF)] = self.Rn
        self.REACT = self.R_structure.reshape(N_node,6)
    
    def force_recovery(self):
        """pass displacements to each element and recover element internal forces"""
        # calculate element internal forces
        for ele in self.element_list:
            d_ele = self.d_structure[ele.dof]
            ele.force_recovery(d_ele)
        
        # assemble ELE_FOR matrix
        self.ELE_FOR = np.zeros([len(self.element_list),12])
        for i in range(len(self.element_list)):
            self.ELE_FOR[i,:] = self.element_list[i].f_local
        
    def disp_recovery(self):
        """pass displacements to each node for storage"""        
        for i in range(len(self.node_list)):
            self.node_list[i].disp = self.DEFL[i,:]
            self.node_list[i].reaction = self.REACT[i,:]
            self.node_list[i].coord_disp = self.node_list[i].coord + self.DEFL[i,0:3]
    
    def compute_error(self):
        """back-calculate residual vector to quantify error"""
        # back-calculated load vector
        P_actual = self.K_structure @ self.d_structure
        
        # assembled load vector with FEF, reactions, and load on support
        P_original = self.P_structure - self.FEF_structure + self.R_structure
        P_original[np.ix_(self.fixedDOF)] = P_original[np.ix_(self.fixedDOF)] + self.Ps
        P_original[np.ix_(self.dispDOF)] = P_original[np.ix_(self.dispDOF)] + self.Pn
        
        # compute residual
        self.ERROR = P_actual - P_original
        

    
    
    
    
    