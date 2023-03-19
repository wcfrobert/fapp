<img src="C:\Users\wcfro\OneDrive\Desktop\Personal Projects\fapp\docs\logo.png" alt="logo" style="zoom:50%;" />

[INSERT LOGO]

A lightweight and incredibly <u>easy-to-use</u> Python implementation of the Finite Element Method (FEM). Perform quick first-order elastic frame analyses right in Python.

[INSERT DEMO GIF]

[INSERT TABLE OF CONTENT]

## Introduction

fapp (pronounced "F - app"), stands for Frame Analysis Program in Python. As the name suggests, only line elements are supported at this time (no surface or solid elements). The entire code base is only about 500 lines (not counting comments and the visualization functionalities). I'm sure you'll probably find it satisfactory for most use-cases within the context of structural engineering. It was a joy to develop and I hope you'll love how simple it is to use!

The module has an object-oriented design that is quite intuitive. I made a conscious effort during development to keep things pedagogically easy to understand. I hope it will serve as excellent reference material for educators looking for a python implementation of the direct stiffness method. 

* Analyze any 3-D frame structures
* Timoshenko beam elements with shear deformation
* First-order elastic analyses
* Nodal and member loads
* Imposed displacements
* Object-oriented design. Easy to access and store results.
* Beautiful and fully interactive visualization in your web browser. Pan, zoom, orbit, and hover.

Disclaimer: this package is meant for personal or educational use only. As evident from the empty "\tests" folder, fapp is NOT robust enough to be used for commercial purposes of any kind! I'm sure the results are correct for typical use, but there are plenty of edge cases that I simply haven't had the chance to explore/debug.



## Quick Start

**Step 0:** Installation. See "Installation" section below for more info

* For casual users, simply install Anaconda Python, download this module, and open "main.py" in Spyder IDE. The following open source packages are used in this project:
  * Numpy
  * Scipy
  * Plotly 


**Step 1:** Import fapp and initialize structure

```python
# import
import fapp.analysis as fpa
import fapp.plotter as fpp

# initialize structure
my_structure = fpa.Analysis()
```

**Step 2:** Define nodes and elements

```python
# node tag must be defined chronologically from 0 to N.
my_structure.add_node(node_tag=0, x=0, y=0, z=0)
my_structure.add_node(node_tag=1, x=0, y=4000, z=0)
my_structure.add_node(node_tag=2, x=12000, y=7000, z=0)
my_structure.add_node(node_tag=3, x=12000, y=0, z=0)

# element tag must be defined chronologically from 0 to N
my_structure.add_element(ele_tag=0, i_tag=0, j_tag=1, A=2e4, Ayy=2e4, Azz=2e4, Iy=999, Iz=1.5e9, J=999, E=200, v=0.3, beta=0)
my_structure.add_element(ele_tag=1, i_tag=1, j_tag=2, A=4e4, Ayy=4e4, Azz=4e4, Iy=999, Iz=3.8e9, J=999, E=200, v=0.3, beta=0)
my_structure.add_element(ele_tag=2, i_tag=2, j_tag=3, A=2e4, Ayy=2e4, Azz=2e4, Iy=999, Iz=1.5e9, J=999, E=200, v=0.3, beta=0)
```

**Step 3:** Define fixity and loading

```python
# define fixity - 0 = fixed, "nan" = free, float = prescribed displacement
my_structure.add_fixity(node_tag=0, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)
my_structure.add_fixity(node_tag=3, ux=0, uy=0, uz=0)

# define loading - nodal or member loads
my_structure.add_load_member(ele_tag=1, wx=-0.00363, wy=-0.01455, wz=0)
my_structure.add_load_nodal(node_tag=1, Fx=5, Fy=0, Fz=0, Mx=0, My=0, Mz=0)
```

**Step 4:** Run analysis

```python
my_structure.solve()
```

**Step 5:** Visualize results

```python
fig1 = fpp.plot(my_structure)
fig2 = fpp.plot_results(my_structure, display="d")
fig3 = fpp.plot_diagrams(my_structure, ele_tag=0)
fig1.show()
fig2.show()
fig3.show()
```



## Installation

**<u>Option 1: Anaconda Python Distribution</u>**

For the casual users, Anaconda Python distribution is recommended. The bundle of packages that comes with the Anaconda is more than enough. This is by far the easiest method of installation. Users don't need to worry about dependency management and setting up virtual environments. The following open source packages are used in this project:

* Numpy
* Scipy
* Plotly 

Installation procedure:

1. Download Anaconda python
2. Download this package (click the green "Code" button and download zip file)
3. Open "main.py" in Anaconda's Spyder IDE. Make sure working directory is correctly configured.
4. Run

**<u>Option 2: Command Prompt + Plain Python</u>**

Pip install is available.

```
pip install fapp
```

Alternatively:

1. Download this project to a folder of your choosing
    ```
    git clone https://github.com/wcfrobert/wcfrobert.github.io
    ```
2. Create virtual environment
    ```
    py -m venv venv
    ```
3. Activate virtual environment
    ```
    venv/Scripts/activate
    ```
4. Install requirements
    ```
    pip install -r requirements.txt
    ```
5. fapp is ready for use
    ```
    py main.py
    ```




## Usage

### Step 0: Instantiate analysis object

fapp has an object-oriented design that is quite intuitive to understand and use. Start by instantiating an "Analysis" object. 

```python
import fapp.Analysis

my_structure = fapp.analysis.Analysis()
```



### Step 1: Add nodes

Add nodes to your structure with ".add_node" method. Please note **node tag must be defined chronologically** from 0 to N. This is because the node tag should correspond to its index in a list that stores all the nodes in the structure. Furthermore, degree of freedoms (DOFs) are numbered chronologically (also known as the plain numberer). It would have been simple enough to add a mapping step such that the node tag (used to determine DOF) is decoupled from a user-specified label, but that adds a layer of obfuscation that I was not particularly fond of.

`Analysis.add_node(node_tag, x, y, z)` adds a node to the structure. Parameters:

* node_tag: int
  * Node tag must be defined chronologically from 0 to N. Node 0 will have DOF [0,1,2,3,4,5], node 1 will have DOF [6,7,8,9,10,11], and so on...
* x: float
  * x coordinate
* y: float
  * y coordinate
* z: float
  * z coordinate

```python 
# named arguments recommended for clarity
my_structure.add_node(0, 10, 0, 0)
my_structure.add_node(node_tag=1, x=0, y=0, z=0)

# nodes must be defined chronologically starting at 0. THe following is NOT allowed:
my_structure.add_node(node_tag=1, x=0, y=0, z=0) # ERROR! Must start from 0!
my_structure.add_node(node_tag=999, x=0, y=0, z=0) # ERROR! Second node so tag = 1
```



### Step 2: Add elements

Once all the nodes have been defined, create elements by linking nodes together. In this step, we are defining both the connectivity as well as section/material properties. Again, **element tag must be defined chronologically** from 0 to N.

`Analysis.add_element(ele_tag,i_tag,j_tag,A,Ayy,Azz,Iy,Iz,J,E,v=0.3,beta=0)` adds an element to the structure. Parameters:

* ele_tag: int
  * Element tag must be defined chronologically from 0 to N
* i_tag: int
  * Start node tag
* j_tag: int
  * End node tag
* A: float
  * Section area
* Ayy: float
  * Shear area along y axis (major bending direction). Input very large value to neglect shear deformation
* Azz: float
  * Shear area along z axis (minor bending direction). Input very large value to neglect shear deformation
* Iy: float
  * Moment of inertia about y axis (minor bending direction)
* Iz: float
  * Moment of inertia about z axis (major bending direction)
* J: float
  * Torsion constant
* E: float
  * Young's modulus
* v: float (optional)
  * Poisson's ratio. Default = 0.3. Used to determine shear modulus $G = frac{E}{2(1+v)}$
* beta: float (optional)
  * Member rotation along its longitudinal axis in degrees. Default  = 0. Refer to "Notes and Assumptions" section for default geometric transformation

```python
# Named arguments recommended for clarity
my_structure.add_element(0, 0, 1, 2e4, 2e4, 2e4, 999, 1.5e9, 999, 200, 0.3, 0)
my_structure.add_element(ele_tag=1, i_tag=1, j_tag=2, A=4e4, Ayy=4e4, Azz=4e4,
                         Iy=999, Iz=3.8e9, J=999, E=200, v=0.3, beta=0)
```



### Step 3: Add fixity

Apply boundary conditions. You may wish to fix a node, or impose a prescribed displacement. Enter 0 for fixed, "nan" for free, float for prescribed displacement.

`Analysis.add_fixity(node_tag, ux="nan", uy="nan", uz="nan", rx="nan", ry="nan", rz="nan")` adds fixity to the structure. Parameters:

* node_tag: int
  * Node you wish to fix
* ux: int, float, string (optional)
  * Displacement in global X direction. Default = "nan"
* uy: int, float, string (optional)
  * Displacement in global Y direction. Default = "nan"
* uz: int, float, string (optional)
  * Displacement in global Z direction. Default = "nan"
* rx: int, float, string (optional)
  * Rotation about global X axis. Default = "nan"
* ry: int, float, string (optional)
  * Rotation about global Y axis. Default = "nan"
* rz: int, float, string (optional)
  * Rotation about global Z axis. Default = "nan"

```python
# fully fixed support
my_structure.add_fixity(node_tag=0, ux=0, uy=0, uz=0, rx=0, ry=0, rz=0)

# pin support. Note you only need to specify DOFs to fix. "nan" is default input
my_structure.add_fixity(node_tag=3, ux=0, uy=0, uz=0)

# fix displacement in X, Y, Z and rotation about X and Y. 2D pin if modeling in XY plane
my_structure.add_fixity(node_tag=4, ux=0, uy=0, uz=0, rx=0, ry=0)
```



### Step 4: Add Loading

Apply loading to member or node.

`Analysis.add_load_member(ele_tag, wx=0, wy=0, wz=0)` adds uniform load to a member. Defined with respect to its **local axis**. Parameters:

* ele_tag: int
  * Element tag to identify element you wish to load
* wx: float (optional)
  * uniform load along member local x axis. Axial direction. Default = 0
* wy: float (optional)
  * uniform load along member local y axis. Major bending direction. Don't forget to put negative number for gravity load. Default = 0
* wz: float (optional)
  * uniform load along member local z axis. Minor bending direction. Default = 0

```python
# apply uniform load of -4.5 to element 1 in its local y axis
my_structure.add_load_member(ele_tag=1, wy=-0.01455)
```



`Analysis.add_load_nodal(node_tag, Fx=0, Fy=0, Fz=0, Mx=0, My=0, Mz=0)` adds external nodal load. Defined with respect to global axis. Parameters:

* node_tag: int
  * Node tag to identify node you wish to load
* Fx: float (optional)
  * Nodal load in global X direction. Default = 0
* Fy: float (optional)
  * Nodal load in global Y direction. Default = 0
* Fz: float (optional)
  * Nodal load in global Z direction. Default = 0
* Mx: float (optional)
  * Concentrated moment about global X axis. Default = 0
* My: float (optional)
  * Concentrated moment about global Y axis. Default = 0
* Mz: float (optional)
  * Concentrated moment about global Z axis. Default = 0

```python
# Apply nodal load of 5.5 to node 1 in global X direction
my_structure.add_load_nodal(node_tag=1, Fx=5)
```



### Step 5: Solve And Post-Process

When ready, start analysis.

`Analysis.solve(print_info=True)` starts analysis routine.

* print_info: Boolean (optional)
  * Set to False if you do not want anything printed to terminal

The beauty of fapp's object-oriented design is that all the information related to the structure, from geometry, connectivity, to force/displacement results will be stored in one variable. The "Analysis" object doubles as a factory and stores a list of "Node" objects and "Element" objects. You can use dot notation to access everything. For example:

```python
# retrieve x,y,z coordinate of node 2
my_structure.node_list[2].coord

# view the stiffness matrix of element 0
my_structure.element_list[0].k_local

# view nodal displacement at node 6 (after solving)
my_structure.node_list[6].disp
```

For ease of access, results are collected and stored in matrices. 

```python
# N_node x 6 matrix containing nodal displacements
my_structure.DEFL
# N_node x 6 matrix containing reaction forces at fixed nodes
my_structure.REACT
# N_element x 12 matrix containing member-end forces in local coordinate
my_structure.ELE_FOR
# for example, return nodal displacement at node 3
my_structure.DEFL[3, :]
```

The use of Spyder's variable explorer is highly recommended!

[INSERT SCREENSHOT OF SPYDER]



### Step 6: Visualize

To help with debugging, I built a visualization engine in Plotly. I must say it turned out much better than expected. You can zoom, pan, orbit, and hover over nodes and elements for more information. I also added camera buttons to help user navigate between views.

There are currently three visualization options:

`fapp.plotter.plot(structure)` for visualizing geometry, node and element tags, and connectivity.

* structure: {fapp.Analysis}
  * Pass in the structure you want to visualize

[INSERT SAMPLE]

* Node hoverinfo shows tag, xyz coordinate, fixity, and loading
* Element hoverinfo shows tag, start and end node, length, member loading, and member rotation





`fapp.plotter.plot_results(structure, display="d", scale=15, show_value = False)` for results after analysis such as displacement and force diagrams.

* structure: {fapp.Analysis}
  * Pass in the structure you want to visualize
* display: string (optional)
  * Indicate what result you wish to see. Default = "d"
  * "d" - displacement
  * "N" - axial
  * "Vy" - shear in local y axis (major bending direction)
  * "Vz" - shear in local z axis (minor bending direction)
  * "Tx" - torsion
  * "My" - moment about local y axis (minor bending direction)
  * "Mz" - moment about local z axis (major bending direction)
* scale: float (optional)
  * Scale of displacement or force diagram. Adjust manually as needed. 15 seems to be a good number. Larger number = larger scale. Default = 15
* show_value: boolean (optional)
  * Show force diagram magnitude at the two ends. Gets quite cluttered so it is off by default. Default = False

[INSERT SAMPLE]

* Node hoverinfo shows tag, displacements, and reactions if applicable
* Element hoverinfo shows tag, and element end forces





`fapp.plotter.plot_diagrams(structure, ele_tag)` for displacement and force diagrams for a single element about its local coordinate system.

* structure: {fapp.Analysis}
  * Pass in the structure you want to visualize
* ele_tag: int
  * Element tag to identify which element you would like to plot

[INSERT SAMPLE]

* You may click on the legend to show/hide the diagrams you wish to see.





## Notes and Assumptions

[INSERT DEFAULT COORDINATE ASSUMPTION SNIP]

* Global Coordinate (X, Y, Z): 
  * Y is the vertical axis (Elevation)
  * X and Z are the axes within the horizontal plane (Plan)
  * Recommend modeling 2-D plane structures in the X-Y plane
  * Recommend modeling 2D floor grids in the X-Z plane
* Local Coordinate (x, y, z):
  * x-axis = element longitudinal axis
  * z-axis = element major bending axis (relevant section properties: Iz, Ayy)
  * y-axis = element minor bending axis (relevant section properties: Iy, Azz)
* Default geometric transformation (i.e. member orientation):
  * Element major axis (web direction vector) always points skyward. 
  * For elements that are completely vertical (columns), major axis always points towards -X. 
  * This default setting works well because beams will always have major axis resisting gravity load, and columns will always have major axis bending properties when modeling within the X-Y plane.
* fapp is agnostic when it comes to unit. Please ensure input is consistent




## License

MIT License

Copyright (c) 2023 Robert Wang