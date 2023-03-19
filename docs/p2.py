# Import
import fapp.analysis as uda
import fapp.plotter as udp


def main():
    # initialize structure
    structure1 = uda.Analysis()
    
    # define nodes
    structure1.add_node(node_tag=0, x=0, y=0, z=0)
    structure1.add_node(node_tag=1, x=5000, y=0, z=0)
    structure1.add_node(node_tag=2, x=5000, y=0, z=-8000)
    structure1.add_node(node_tag=3, x=2000, y=0, z=-8000)
    
    # define elements
    structure1.add_element(ele_tag=0, i_tag=0, j_tag=1, A=1.1e4, Ayy=1.1e4, Azz=1.1e4, Iy=1.74e7, Iz=1.06e8, J=5.29e7,E=200, v=0.3, beta=0)
    structure1.add_element(ele_tag=0, i_tag=1, j_tag=2, A=1.1e4, Ayy=1.1e4, Azz=1.1e4, Iy=1.74e7, Iz=1.06e8, J=5.29e7,E=200, v=0.3, beta=0)
    structure1.add_element(ele_tag=0, i_tag=2, j_tag=3, A=1.1e4, Ayy=1.1e4, Azz=1.1e4, Iy=1.74e7, Iz=1.06e8, J=5.29e7,E=200, v=0.3, beta=0)
    
    # define fixity (0 = fixed, "nan" = free, float = prescribed displacement)
    structure1.add_fixity(node_tag=0, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
    structure1.add_fixity(node_tag=3, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
    
    # define loading
    structure1.add_load_member(ele_tag=1, wx=0, wy=-0.005, wz=0)
    
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