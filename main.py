# Import
import fapp.analysis as fpa
import fapp.plotter as fpp

# initialize structure
my_structure = fpa.Analysis()

# define nodes - node tag must be defined chronologically from 0 to N. Named arguments recommended for clarity
my_structure.add_node(node_tag=0, x=0, y=0, z=0)
my_structure.add_node(node_tag=1, x=0, y=4000, z=0)
my_structure.add_node(node_tag=2, x=12000, y=7000, z=0)
my_structure.add_node(node_tag=3, x=12000, y=0, z=0)

# define elements - element tag must be defined chronologically from 0 to N. Named arguments recommended for clarity
my_structure.add_element(ele_tag=0, i_tag=0, j_tag=1, A=2e4, Ayy=2e4, Azz=2e4, Iy=999, Iz=1.5e9, J=999, E=200, v=0.3, beta=0)
my_structure.add_element(ele_tag=1, i_tag=1, j_tag=2, A=4e4, Ayy=4e4, Azz=4e4, Iy=999, Iz=3.8e9, J=999, E=200, v=0.3, beta=0)
my_structure.add_element(ele_tag=2, i_tag=2, j_tag=3, A=2e4, Ayy=2e4, Azz=2e4, Iy=999, Iz=1.5e9, J=999, E=200, v=0.3, beta=0)

# define fixity - 0 = fixed, "nan" = free, float = prescribed displacement
my_structure.add_fixity(node_tag=0, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
my_structure.add_fixity(node_tag=3, ux=0, uy=0, uz=0)

# define loading - nodal or member loads
my_structure.add_load_member(ele_tag=1, wx=-0.00363, wy=-0.01455, wz=0)
my_structure.add_load_nodal(node_tag=1, Fx=5, Fy=0, Fz=0, Mx=0, My=0, Mz=0)

# start analysis
my_structure.solve()

# visualization options
fig1 = fpp.plot(my_structure)
fig2 = fpp.plot_results(my_structure, display="d")
fig3 = fpp.plot_diagrams(my_structure, ele_tag=1)

import plotly.io as pio
pio.renderers.default='browser'
fig1.show()
fig2.show()
fig3.show()



    
    







