# Import
import fapp.analysis as uda
import fapp.plotter as udp

# Options
STIFFNESS = 1e16


def main():
    # initialize structure
    structure1 = uda.Analysis()
    
    # define nodes
    structure1.add_node(node_tag=0, x=0, y=0, z=0)
    structure1.add_node(node_tag=1, x=1, y=0, z=0)
    structure1.add_node(node_tag=2, x=2, y=0, z=0)
    structure1.add_node(node_tag=3, x=3, y=0, z=0)
    
    # define elements
    structure1.add_element(ele_tag=0, i_tag=0, j_tag=1, A=1, Ayy=1, Azz=1, Iy=1, Iz=1, J=1,E=1, v=0.3, beta=0)
    structure1.add_element(ele_tag=1, i_tag=1, j_tag=2, A=STIFFNESS, Ayy=1, Azz=1, Iy=1, Iz=1, J=1,E=1, v=0.3, beta=0)
    structure1.add_element(ele_tag=2, i_tag=2, j_tag=3, A=1, Ayy=1, Azz=1, Iy=1, Iz=1, J=1,E=1, v=0.3, beta=0)
    
    # define fixity (0 = fixed, "nan" = free, float = prescribed displacement)
    structure1.add_fixity(node_tag=0, ux=0, uy=0, uz=0, rx=0, ry=0, rz="nan")
    structure1.add_fixity(node_tag=1, ux="nan", uy=0, uz=0, rx=0, ry=0, rz="nan")
    structure1.add_fixity(node_tag=2, ux="nan", uy=0, uz=0, rx=0, ry=0, rz="nan")
    structure1.add_fixity(node_tag=3, ux="nan", uy=0, uz=0, rx=0, ry=0, rz="nan")
    
    # define loading
    structure1.add_load_nodal(node_tag=1, Fx=1)
    structure1.add_load_nodal(node_tag=3, Fx=1)
    
    # start analysis
    structure1.solve()
    
    # visualization options
    fig1 = udp.plot(structure1)
    fig2 = udp.plot_results(structure1, display="d")
    fig3 = udp.plot_diagrams(structure1, ele_tag=0)

    return structure1,fig1,fig2,fig3


    
    
















if __name__ == "__main__":
    structure,fig1,fig2,fig3 = main()
else: 
    print('module imported') 