# Import
import fapp.analysis as uda
import fapp.plotter as udp

# Options
CONSIDER_SHEAR_DEFORMATION = True


def main():
    # initialize structure
    structure1 = uda.Analysis()
    
    # define nodes
    structure1.add_node(node_tag=0, x=0, y=0, z=0)
    structure1.add_node(node_tag=1, x=0, y=156, z=0)
    structure1.add_node(node_tag=2, x=0, y=312, z=0)
    structure1.add_node(node_tag=3, x=0, y=468, z=0)
    
    structure1.add_node(node_tag=4, x=360, y=0, z=0)
    structure1.add_node(node_tag=5, x=360, y=156, z=0)
    structure1.add_node(node_tag=6, x=360, y=312, z=0)
    structure1.add_node(node_tag=7, x=360, y=468, z=0)
    
    structure1.add_node(node_tag=8, x=720, y=0, z=0)
    structure1.add_node(node_tag=9, x=720, y=156, z=0)
    structure1.add_node(node_tag=10, x=720, y=312, z=0)
    structure1.add_node(node_tag=11, x=720, y=468, z=0)
    
    structure1.add_node(node_tag=12, x=1080, y=0, z=0)
    structure1.add_node(node_tag=13, x=1080, y=156, z=0)
    structure1.add_node(node_tag=14, x=1080, y=312, z=0)
    structure1.add_node(node_tag=15, x=1080, y=468, z=0)
    
    Ayy=1e9
    Azz=1e9
    # define elements
    # W14x211: A=62.0, Ayy=15.4, Azz=49.3, Iy=1030, Iz=2660, J=44.6, E=30000
    if CONSIDER_SHEAR_DEFORMATION:
        Ayy=15.4
        Azz=49.3
    structure1.add_element(ele_tag=0, i_tag=0, j_tag=1, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=1, i_tag=1, j_tag=2, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=2, i_tag=2, j_tag=3, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=3, i_tag=4, j_tag=5, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=4, i_tag=5, j_tag=6, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=5, i_tag=6, j_tag=7, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=6, i_tag=8, j_tag=9, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=7, i_tag=9, j_tag=10, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=8, i_tag=10, j_tag=11, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=9, i_tag=12, j_tag=13, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=10, i_tag=13, j_tag=14, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    structure1.add_element(ele_tag=11, i_tag=14, j_tag=15, A=62.0, Ayy=Ayy, Azz=Azz, Iy=1030, Iz=2660, J=44.6, E=30000)
    
    # W36x150: A=44.3, Ayy=22.4, Azz=22.5, Iy=270, Iz=9040, J=10.1, E=30000
    if CONSIDER_SHEAR_DEFORMATION:
        Ayy=22.4
        Azz=22.5
    structure1.add_element(ele_tag=12, i_tag=1, j_tag=5, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=13, i_tag=5, j_tag=9, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=14, i_tag=9, j_tag=13, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=15, i_tag=2, j_tag=6, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=16, i_tag=6, j_tag=10, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=17, i_tag=10, j_tag=14, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=18, i_tag=3, j_tag=7, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=19, i_tag=7, j_tag=11, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    structure1.add_element(ele_tag=20, i_tag=11, j_tag=15, A=44.3, Ayy=Ayy, Azz=Azz, Iy=270, Iz=9040, J=10.1, E=30000)
    
    # define fixity (0 = fixed, "nan" = free, float = prescribed displacement)
    structure1.add_fixity(node_tag=0, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
    structure1.add_fixity(node_tag=4, ux=0, uy=-1., uz=0, rx=0, ry=0, rz=0)
    structure1.add_fixity(node_tag=8, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
    structure1.add_fixity(node_tag=12, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
    
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