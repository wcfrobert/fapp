# Import
import fapp.analysis as uda
import fapp.plotter as udp
import numpy as np

# Options
N_SEG = 128

def main():
    ele_counter = 0
    node_counter = 7
    
    # initialize structure
    structure1 = uda.Analysis()
    
    # define nodes
    structure1.add_node(node_tag=0, x=0, y=0, z=0)
    structure1.add_node(node_tag=1, x=2000, y=0, z=0)
    structure1.add_node(node_tag=2, x=1000, y=2500, z=0)
    structure1.add_node(node_tag=3, x=1000, y=2500, z=1500)
    structure1.add_node(node_tag=4, x=0, y=0, z=3000)
    structure1.add_node(node_tag=5, x=2000, y=0, z=3000)
    structure1.add_node(node_tag=6, x=1000, y=2500, z=3000)
    
    # define elements
    ele_counter,node_counter = add_segment(structure1,0,2,ele_counter,node_counter, 
                                           A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
    ele_counter,node_counter = add_segment(structure1,1,2,ele_counter,node_counter, 
                                           A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
    ele_counter,node_counter = add_segment(structure1,2,3,ele_counter,node_counter, 
                                           A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
    ele_counter,node_counter = add_segment(structure1,3,6,ele_counter,node_counter, 
                                           A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
    ele_counter,node_counter = add_segment(structure1,4,6,ele_counter,node_counter, 
                                           A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
    ele_counter,node_counter = add_segment(structure1,5,6,ele_counter,node_counter, 
                                           A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
    
    
    # define fixity (0 = fixed, "nan" = free, float = prescribed displacement)
    structure1.add_fixity(node_tag=0, ux=0, uy=0, uz=0, rx="nan", ry="nan", rz="nan")
    structure1.add_fixity(node_tag=1, ux=0, uy=0, uz=0, rx="nan", ry="nan", rz="nan")
    structure1.add_fixity(node_tag=4, ux=0, uy=0, uz=0, rx="nan", ry="nan", rz="nan")
    structure1.add_fixity(node_tag=5, ux=0, uy=0, uz=0, rx="nan", ry="nan", rz="nan")
    
    # define loading
    structure1.add_load_nodal(node_tag=3, Fx=4.5, Fy=0, Fz=0, Mx=0, My=0, Mz=0)
    
    # start analysis
    structure1.solve()
    
    # visualization options
    #fig1 = udp.plot(structure1)
    #fig2 = udp.plot_results(structure1, display="d")
    #fig3 = udp.plot_diagrams(structure1, ele_tag=0)
    

    return structure1


    
def add_segment(structure, start_node, end_node, ele_counter, node_counter,
                A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200):
    if N_SEG == 1:
        structure.add_element(ele_tag=ele_counter, i_tag=start_node, j_tag=end_node,
                               A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
        ele_counter += 1
    else:
        # segments in local domain from 0 to 1
        du = np.linspace(0,1,N_SEG+1)
        du_int = du[1:-1]
        
        # extract end coordinates
        coord0 = structure.node_list[start_node].coord
        coord1 = structure.node_list[end_node].coord
        
        # connectivity list
        connect = []
        connect.append(start_node)
        
        for i in range(len(du_int)):
            u0 = du_int[i]
            x0 = coord0[0] * (1 - u0) + coord1[0] * (u0)
            y0 = coord0[1] * (1 - u0) + coord1[1] * (u0)
            z0 = coord0[2] * (1 - u0) + coord1[2] * (u0)
            
            connect.append(node_counter)
            structure.add_node(node_tag=node_counter, x=x0, y=y0, z=z0)
            node_counter += 1
            
        connect.append(end_node)
        
        # chain together nodes
        for i in range(len(connect)-1):
            node_i = connect[i]
            node_j = connect[i+1]
            structure.add_element(ele_tag=ele_counter, i_tag=node_i, j_tag=node_j,
                                   A=1430, Ayy=1430, Azz=1430, Iy=1.26e6, Iz=1.26e6, J=2.52e6,E=200)
            ele_counter += 1
        
    
    return ele_counter, node_counter



if __name__ == "__main__":
    structure = main()
else: 
    print('module imported') 