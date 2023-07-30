import numpy as np
import plotly.graph_objects as go


def plot(structure):
    """visualize geometry, connectivity, dofs of structure"""
    fig = go.Figure()
    
    # plot nodes
    xlist,ylist,zlist=[],[],[]
    for n in structure.node_list:
        draw_node(fig, n, hoverinfo="geometry", showID=True)
        xlist.append(n.coord[0])
        ylist.append(n.coord[1])
        zlist.append(n.coord[2])
    
    # plot elements
    for e in structure.element_list:
        draw_element(fig, e, hoverinfo="geometry", showID=True)
    
    # for scaling purposes
    dx = max(xlist) - min(xlist)
    dy = max(ylist) - min(ylist)
    dz = max(zlist) - min(zlist)
    
    # apply some style and axes scaling
    my_scene = dict(
        aspectmode="data",
        xaxis = dict(showticklabels=False,gridcolor='white',backgroundcolor='white',showspikes=False,visible=False),
        yaxis = dict(showticklabels=False,gridcolor='white',backgroundcolor='white',showspikes=False,visible=False),
        zaxis = dict(showticklabels=False,gridcolor='white',backgroundcolor='white',showspikes=False,visible=False))
    fig.update_layout(scene=my_scene,showlegend=False,hoverlabel=dict(bgcolor="white"))
    my_camera = dict(up=dict(x=0,y=1,z=0),center=dict(x=0,y=0,z=0),eye=dict(x=1.25,y=1.25,z=1.25))
    fig.update_layout(scene_camera=my_camera,scene_dragmode='pan',
                      margin=dict(l=50,r=50,b=50,t=10,pad=4))
    
    # plot origin axis marker
    draw_origin_marker(fig, dx, dy, dz)
    
    # add buttons for camera position
    draw_camera_buttons(fig)
    
    return fig




def plot_diagrams(structure, ele_tag):
    """plot displacements and forces for a specified element in its local coordinate"""
    # retrieve element info
    ele = structure.element_list[ele_tag]
    L = ele.length
    i_node = ele.i_node.tag
    j_node = ele.j_node.tag
    
    # interpolate
    plot_x, plot_y_disp, plot_y_force = interpolate_internal(ele, N_pts=50)
    
    # plot diagram
    fig = go.Figure()
    
    # draw black line along x-axis to show member
    fig.add_trace(go.Scatter(x=[0,L],y=[0,0],
                            mode='lines+markers+text',name="Member",
                            line = dict(color='black',width = 2),
                            marker = dict(size = 7),
                            hoverinfo = 'skip',
                            text=[f"Node {i_node}",f"Node {j_node}"],
                            textposition="top center",
                            textfont=dict(family="Arial",size=14,color="black")))

    # draw diagrams
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_force[:,0],mode='lines+markers',fill='tozeroy',name = 'N(x)',
                             visible = "legendonly", marker_size=4))
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_force[:,1],mode='lines+markers',fill='tozeroy',name = 'Vy(x)',
                             visible = "legendonly", marker_size=4))
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_force[:,2],mode='lines+markers',fill='tozeroy',name = 'Vz(x)',
                             visible = "legendonly", marker_size=4))
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_force[:,3],mode='lines+markers',fill='tozeroy',name = 'Tx(x)',
                             visible = "legendonly", marker_size=4))
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_force[:,4],mode='lines+markers',fill='tozeroy',name = 'My(x)',
                             visible = "legendonly", marker_size=4))
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_force[:,5],mode='lines+markers',fill='tozeroy',name = 'Mz(x)',
                             visible = "legendonly", marker_size=4))
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_disp[:,1],mode='lines+markers',name = 'y disp',
                             marker_size=4))
    fig.add_trace(go.Scatter(x=plot_x,y=plot_y_disp[:,2],mode='lines+markers',name = 'z disp',
                             visible = "legendonly", marker_size=4))

    # styling
    fig.update_layout(title="Element {} Diagrams".format(ele.tag),
                      xaxis_title="Along Member Local x",
                      yaxis_title="Force or displacement Magnitude",
                      title_xanchor="center",
                      title_x=0.5)
    return fig



def plot_results(structure, display="d", scale=15, show_value = False):
    """
    Same as plot() but change hoverinfo. Also draw lines for displacement shape or force diagram
    Display modes:
        "d"   -   displaced shape
        "N"   -   axial force
        "Vy"  -   shear force (major axis)
        "Vz"  -   shear force (minor axis)
        "Tx"  -   torsion
        "My"  -   moment (minor axis)
        "Mz"  -   moment (major axis)
    """
    fig = go.Figure()
    
    # plot nodes
    xlist,ylist,zlist=[],[],[]
    for n in structure.node_list:
        draw_node(fig, n, hoverinfo="results")
        xlist.append(n.coord[0])
        ylist.append(n.coord[1])
        zlist.append(n.coord[2])
    
    # plot elements
    for e in structure.element_list:
        draw_element(fig, e, hoverinfo="results")
    
    # for scaling purposes
    dx = max(xlist) - min(xlist)
    dy = max(ylist) - min(ylist)
    dz = max(zlist) - min(zlist)
    
    # apply some style and axes scaling
    my_scene = dict(
        aspectmode="data",
        xaxis = dict(showticklabels=False,gridcolor='white',backgroundcolor='white',showspikes=False,visible=False),
        yaxis = dict(showticklabels=False,gridcolor='white',backgroundcolor='white',showspikes=False,visible=False),
        zaxis = dict(showticklabels=False,gridcolor='white',backgroundcolor='white',showspikes=False,visible=False))
    fig.update_layout(scene=my_scene,showlegend=False,hoverlabel=dict(bgcolor="white"))
    my_camera = dict(up=dict(x=0,y=1,z=0),center=dict(x=0,y=0,z=0),eye=dict(x=1.25,y=1.25,z=1.25))
    fig.update_layout(scene_camera=my_camera,scene_dragmode='pan',
                      margin=dict(l=50,r=50,b=50,t=10,pad=4))
    
    # plot origin axis marker
    draw_origin_marker(fig, dx, dy, dz)
    
    # add buttons for camera position
    draw_camera_buttons(fig)
    
    
    # ---------------- draw line work for force or displacement -----------------
    # add title
    title_dict = {"N": "Axial (N)",
                  "Vy": "Major Shear (Vy)",
                  "Vz": "Minor Shear (Vz)",
                  "Tx": "Torsion (Tx)",
                  "My": "Minor Bending Moment (My)",
                  "Mz": "Major Bending Moment (Mz)",
                  "d": "Displacement (d)"}
    appropriate_title = title_dict[display]
    fig.update_layout(title="{} Diagrams".format(appropriate_title),
                      title_xanchor="center",title_x=0.5, title_y=0.95)
    
    # Need to store max magnitude across all element for consistent scaling
    max_ordinate_all = [[],[],[],[],[],[],[],[],[]] 
    
    # extract interpolated quantities from all elements
    N_pts = 20
    x_local_all = []
    y_disp_all = []
    y_force_all = []
    for ele in structure.element_list:
        x_local, y_disp, y_force = interpolate_internal(ele, N_pts)
        x_local_all.append(x_local)
        y_disp_all.append(y_disp)
        y_force_all.append(y_force)
        max_ordinate_all[0].append(  max(abs(y_disp[:,0]))  )
        max_ordinate_all[1].append(  max(abs(y_disp[:,1]))  )
        max_ordinate_all[2].append(  max(abs(y_disp[:,2]))  )
        max_ordinate_all[3].append(  max(abs(y_force[:,0]))  )
        max_ordinate_all[4].append(  max(abs(y_force[:,1]))  )
        max_ordinate_all[5].append(  max(abs(y_force[:,2]))  )
        max_ordinate_all[6].append(  max(abs(y_force[:,3]))  )
        max_ordinate_all[7].append(  max(abs(y_force[:,4]))  )
        max_ordinate_all[8].append(  max(abs(y_force[:,5]))  )

    # plot results onto figure
    if display == "d":
        Qmax = max(dx,dy,dz) / (100/scale)
        Ymax = max(max(max_ordinate_all[0]),max(max_ordinate_all[1]),max(max_ordinate_all[2]))
        SF = Qmax/Ymax
        
        for i in range(len(x_local_all)):
            # get coordinates to plot in local domain
            ele = structure.element_list[i]
            dx_loc = [a*SF for a in y_disp_all[i][:,0]]
            x_loc = x_local_all[i] + np.array(dx_loc)
            y_loc = [a*SF for a in y_disp_all[i][:,1]]
            z_loc = [a*SF for a in y_disp_all[i][:,2]]
            
            # obtain rotation matrix
            T_rotate = ele.T[0:3,0:3]
            T_rotate = np.transpose(T_rotate)
            
            # obtain translation
            T_translate = ele.i_node.coord
            
            # rotate and translate coordinates
            x_global, y_global, z_global = coord_transform(x_loc,y_loc,z_loc,T_rotate,T_translate)
            
            # add to plot
            plot_data = go.Scatter3d(x=x_global, y=y_global, z=z_global,
                mode='lines',line_color="red", hoverinfo='skip')
            fig.add_trace(plot_data) 
    else:
        # for force diagrams, different index depending on "display"
        display_dict = {"N":0,"Vy":1,"Vz":2,"Tx":3,"My":4,"Mz":5}
        index = display_dict[display]
        
        # scaling same as above
        Qmax = max(dx,dy,dz) / (100/scale)
        Ymax = max(max_ordinate_all[index+3]) #first 3 index x,y,z disp
        SF = 1 if Ymax==0 else Qmax/Ymax
        
        for i in range(len(x_local_all)):
            ele = structure.element_list[i]
            x_loc = x_local_all[i]
            
            # major and minor axis force diagrams are plotted differently
            if display == "Vz" or display == "My":
                z_loc = [a*SF for a in y_force_all[i][:,index]]
                y_loc = np.zeros(len(z_loc))
            else:
                y_loc = [a*SF for a in y_force_all[i][:,index]]
                z_loc = np.zeros(len(y_loc))
                
            # append and prepend
            x_loc = np.insert(x_loc,0,0)
            x_loc = np.append(x_loc,ele.length)
            y_loc = np.insert(y_loc,0,0)
            y_loc = np.append(y_loc,0)
            z_loc = np.insert(z_loc,0,0)
            z_loc = np.append(z_loc,0)
            
            # obtain rotation matrix
            T_rotate = ele.T[0:3,0:3]
            T_rotate = np.transpose(T_rotate)
            
            # obtain translation
            T_translate = ele.i_node.coord
            
            # rotate and translate coordinates
            x_global, y_global, z_global = coord_transform(x_loc,y_loc,z_loc,T_rotate,T_translate)
            
            # add to plot
            plot_data = go.Scatter3d(x=x_global, y=y_global, z=z_global,
                mode='lines',line_color="blue", hoverinfo='skip', line_width=3)
            fig.add_trace(plot_data)
            
            # add data annotation. Get coordinate of end points
            if show_value:
                x_anno = x_global[[1,-2]]
                y_anno = y_global[[1,-2]]
                z_anno = z_global[[1,-2]]
                text = []
                text.append("{:.2f}".format(y_force_all[i][:,index][0]))
                text.append("{:.2f}".format(y_force_all[i][:,index][-1]))
                
                plot_data = go.Scatter3d(x=x_anno,y=y_anno,z=z_anno,
                    mode='markers+text', hoverinfo='skip', text=text, 
                    textfont_color="blue",
                    marker = dict(size = 2,color = "blue"))
                fig.add_trace(plot_data)
            
    return fig





def coord_transform(x_loc,y_loc,z_loc,T_rot,T_trans):
    """takes in local coordinates, transform to global"""
    x_glob, y_glob, z_glob = [],[],[]
    for i in range(len(x_loc)):
        u_loc = np.array([x_loc[i], y_loc[i], z_loc[i]])
        u_glob = T_rot @ u_loc + T_trans
        x_glob.append(u_glob[0])
        y_glob.append(u_glob[1])
        z_glob.append(u_glob[2])
    return np.array(x_glob), np.array(y_glob), np.array(z_glob)




def interpolate_internal(ele, N_pts):
    """
    Use shape functions to interpolate displacements between nodes. See this link:
        https://doc/comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_beam.11.13.html
    Use statics and equilibrium equations to interpolate forces.
    """
    # element information
    L = ele.length
    vec_d = ele.d_local
    vec_f = ele.f_local[0:6]
    wx,wy,wz = ele.load
    
    # set up lists for plotting
    plot_x = []
    plot_y_disp = []
    plot_y_force = []
    
    for i in range(N_pts):
        # global and local domain
        x = i / (N_pts-1) * L
        u = x/L
    
        # shape functions for displacement interpolation (cubic)
        N1 = 1 - u
        N2 = u
        N3 = 1 - 3*u*u + 2*u*u*u
        N4 = 3*u*u - 2*u*u*u
        N5 = L* (u - 2*u*u + u*u*u)
        N6 = L* (-u*u + u*u*u)
        M1 = 1 - u
        M2 = u
        M3 = -6/L * (u - u*u)
        M4 = 6/L * (u - u*u)
        M5 = 1 - 4*u + 3*u*u
        M6 = -2*u + 3*u*u
        N = np.array([
            [N1, 0  , 0  , 0 , 0  , 0 , N2, 0  , 0 , 0 , 0  , 0 ],
            [0 , N3 , 0  , 0 , 0  , N5, 0 , N4 , 0 , 0 , 0  , N6],
            [0 , 0  , N3 , 0 , -N5, 0 , 0 , 0  , N4, 0 , -N6, 0 ],
            [0 , 0  , 0  , M1, 0  , 0 , 0 , 0  , 0 , M2, 0  , 0 ],
            [0 , 0  , M3 , 0 , M5 , 0 , 0 , 0  , M4, 0 , M6 , 0 ],
            [0 , -M3, 0  , 0 , 0  , M5, 0 , -M4, 0 , 0 , 0  , M6]])
        interp_d = N @ vec_d
        
        # cubic interpolation of v => quadratic theta => linear phi => linear M => constant V
        # given our limited scope, use statics to get a better interpolation of internal forces
        Nf = np.array([
            [-1  , 0  , 0   , 0 , 0  , 0 ],
            [0   , 1  , 0   , 0 , 0  , 0 ],
            [0   , 0  , 1   , 0 , 0  , 0 ],
            [0   , 0  , 0   , 1 , 0  , 0 ],
            [0   , 0  , x   , 0 , 1  , 0 ],
            [0   , x  , 0   , 0 , 0  , -1]])
        load_vec = np.array([
            -wx*x,
            wy*x,
            wz*x,
            0,
            wz*x*x/2,
            wy*x*x/2])
        interp_f = Nf @ vec_f + load_vec
        
        # store results
        plot_x.append(x)
        plot_y_disp.append(interp_d)
        plot_y_force.append(interp_f)
        
    plot_y_disp = np.array(plot_y_disp)
    plot_y_force = np.array(plot_y_force)
    
    return plot_x, plot_y_disp, plot_y_force


def draw_node(fig, node, hoverinfo, showID=False):
    """draw node on figure and add hover info"""
    # node coordinates
    x = [node.coord[0]]
    y = [node.coord[1]]
    z = [node.coord[2]]
    
    # hover information
    hovertemplate ='<b>%{text}</b><extra></extra>'
    if hoverinfo == "geometry":
        text = 'Node number: {}'.format(node.tag) +\
                '<br>X: {:.1f}'.format(x[0])+\
                '<br>Y: {:.1f}'.format(y[0])+\
                '<br>Z: {:.1f}'.format(z[0])+\
                '<br>Fixity: {}'.format(node.fixity)+\
                '<br>Loading: {}'.format(node.load)
    elif hoverinfo == "results":
        disp_hover = [round(x,2) for x in node.disp]
        react_hover = [round(x,2) for x in node.reaction]
        text = 'Node number: {}'.format(node.tag) +\
                '<br>Displacement: {}'.format(disp_hover)+\
                '<br>Reaction: {}'.format(react_hover)
    
    # change colors depending on fixity, loading
    if type(node.load) == type(None) and type(node.fixity) == type(None):
        node_color = "black"
    else:
        node_color = "red"
    
    # add trace to figure
    if showID:
        plot_data = go.Scatter3d(x=x,y=y,z=z,
            mode='markers+text', marker = dict(size = 3,color = node_color),
            hovertemplate = hovertemplate, text=[text], 
            texttemplate = f"N{node.tag}", textfont_color="red")
    else:
        plot_data = go.Scatter3d(x=x,y=y,z=z,
            mode='markers', marker = dict(size = 3,color = node_color),
            hovertemplate = hovertemplate, text=[text])
    
    fig.add_trace(plot_data)



def draw_element(fig, element, hoverinfo, showID=False):
    """draw element on figure and add hover info"""
    # member color depending on if member loading exists
    if element.load==[0,0,0]:
    	member_color = "black"
    else:
    	member_color = "black"
        
    # get end node coordinates
    x = [element.i_node.coord[0], element.j_node.coord[0]]
    y = [element.i_node.coord[1], element.j_node.coord[1]]
    z = [element.i_node.coord[2], element.j_node.coord[2]]
    
    # hover information
    hovertemplate ='<b>%{text}</b><extra></extra>'
    if hoverinfo == "geometry":
        text = 'Element number: {}'.format(element.tag) +\
                '<br>Length: {:.1f}'.format(element.length)+\
                '<br>i node: {}'.format(element.i_node.tag)+\
                '<br>j node: {}'.format(element.j_node.tag)+\
                '<br>Beta rotation: {}'.format(element.beta)+\
                '<br>Member load: {}'.format(element.load)
    elif hoverinfo == "results":
        f_local_hover = [round(x,2) for x in element.f_local]
        text = 'Element number: {}'.format(element.tag) +\
                '<br>i node: {}'.format(element.i_node.tag)+\
                '<br>j node: {}'.format(element.j_node.tag)+\
                '<br>Fi: {}'.format(f_local_hover[0:6])+\
                '<br>Fj: {}'.format(f_local_hover[6:])
    
    # plot member
    plot_data = go.Scatter3d(x=x, y=y, z=z,
        mode='lines',line = dict(color=member_color),hoverinfo='skip')
    fig.add_trace(plot_data)
    
    # plot an invisible node at midpoint to display element information
    x_mid = [(x[0] + x[1])/2]
    y_mid = [(y[0] + y[1])/2]
    z_mid = [(z[0] + z[1])/2]
    
    if showID:
        plot_data = go.Scatter3d(x=x_mid, y=y_mid, z=z_mid,
            mode='markers+text', marker = dict(size = 25, color="rgba(255, 255, 255, 0)"),
            hovertemplate = hovertemplate, text = [text],
            texttemplate = f"E{element.tag}")
    else:
        plot_data = go.Scatter3d(x=x_mid, y=y_mid, z=z_mid,
            mode='markers', opacity=0, marker = dict(size = 25),
            hovertemplate = hovertemplate, text = [text])
    
    fig.add_trace(plot_data)



def draw_origin_marker(fig, dx, dy, dz):
    """plot the cartesian coordinate origin marker"""
    # dmax used to scale the marker appropriately
    dmax = max(dx,dy,dz)
    
    # plot hidden node to stretch out axes
    fig.add_trace(go.Scatter3d(x=[dmax],y=[dmax],z=[dmax],mode="markers", opacity=0, hoverinfo="skip"))
    
    # Plot X (in blue)
    X = go.Scatter3d(
        x=[0,dmax/14],
        y=[0,0],
        z=[0,0],
        mode='lines+text',
        hoverinfo = 'skip',
        line=dict(color='blue', width=4),
        text=["","X"],
        textposition="middle right",
        textfont=dict(
            family="Arial",
            size=8,
            color="blue"))
    fig.add_trace(X)
    # Plot Y (in red)
    Y = go.Scatter3d(
        x=[0,0],
        y=[0,dmax/14],
        z=[0,0],
        mode='lines+text',
        hoverinfo = 'skip',
        line=dict(color='red', width=4),
        text=["","Y"],
        textposition="top center",
        textfont=dict(
            family="Arial",
            size=8,
            color="red"))
    fig.add_trace(Y)
    # Plot Z (in green)
    Z = go.Scatter3d(
        x=[0,0],
        y=[0,0],
        z=[0,dmax/14],
        mode='lines+text',
        hoverinfo = 'skip',
        line=dict(color='green', width=4),
        text=["","Z"],
        textposition="middle center",
        textfont=dict(
            family="Arial",
            size=8,
            color="green"))
    fig.add_trace(Z)


def draw_camera_buttons(fig):
    """add some camera buttons to help users navigate between views"""
    button1 = dict(
        method = "relayout",
        args=[{"scene.camera.up": {'x':0,'y':1,'z':0},
                 "scene.camera.eye":{'x':1.25,'y':1.25,'z':1.25},
                 "scene.camera.center":{'x':0,'y':0,'z':0}}], 
        label = "Axonometric")
    button2 = dict(
        method = "relayout",
        args=[{"scene.camera.up":{'x':0, 'y':1,'z':0},
                 "scene.camera.eye":{'x':0,'y':2,'z':0},
                 "scene.camera.center":{'x':0,'y':0,'z':0}}], 
        label="Top View XZ")
    button3 = dict(
        method = "relayout",
        args=[{"scene.camera.up":{'x':0, 'y':1,'z':0}, 
                "scene.camera.eye":{'x':0,'y':0,'z':2},
                "scene.camera.center":{'x':0,'y':0,'z':0}}], 
        label="Plane View XY")
    button4 = dict(
        method = "relayout",
        args=[{"scene.camera.up":{'x':0, 'y':1,'z':0}, 
                "scene.camera.eye":{'x':2,'y':0,'z':0},
                "scene.camera.center":{'x':0,'y':0,'z':0}}], 
        label="Plane View ZY")
    fig.update_layout(
        updatemenus=[
            dict(buttons=[button1, button2, button3, button4], 
                 direction="right",
                 pad={"r": 10, "t": 10},
                 showactive=True,
                 x=0.1,
                 xanchor="center",
                 y=0.1,
                 yanchor="top")
            ]
        )





