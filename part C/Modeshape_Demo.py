import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

from matplotlib.backend_bases import MouseButton
from matplotlib.patches import Rectangle

from functools import partial



#--------------------------------------#
## System settings
m1 = 5 # Mass 1
m2 = 5 # Mass 2

k1 = 50  # Spring 1
k2 = 50  # Spring 2
k3 = 10  # Spring 3

## Animation settings
sim_time = 10  # Time of simulation
fps = 30   # Target FPS

## Plot settings
disp_limit = 2

box_length = 3
box_one_centre = 0

min_spring_length = 2
num_spring_sections = 3

## Plot style
marker_colour = "b"
marker_size = 6
#--------------------------------------#
## Equation of motion
M = np.diag([m1,m2])
K = np.array([[k1+k2,-k2],[-k2,k2+k3]])
(eigenvalues,eigenvectors) = np.linalg.eig(np.linalg.inv(M)@K)
natural_frequency = np.sqrt(eigenvalues)
num_time_points = fps*sim_time
t = np.linspace(0,sim_time,num_time_points)

def modal_displacement(t,initial_x):
    #assume initial velocity is always zero
    Q = np.linalg.inv(eigenvectors)@initial_x
    q = np.array([Q[0]*np.cos(natural_frequency[0]*t),
         Q[1]*np.cos(natural_frequency[1]*t)])
    return q


def modal_to_physical(q):
    return eigenvectors@q

#--------------------------------------#
## Calculate system geometry
box_distance = box_length + 2*disp_limit + min_spring_length
max_disp_left = box_one_centre - disp_limit - min_spring_length
max_disp_right = box_one_centre + box_length + box_distance + disp_limit + min_spring_length
system_disp_lims = [max_disp_left,max_disp_right]
#--------------------------------------#
## Setup figures
plt.close("all")
fig,(ax_system,ax_coords) = plt.subplots(1,2)

ax_system.get_xaxis().set_visible(False)
ax_system.get_yaxis().set_visible(False)
ax_system.set_aspect(1)


ax_system.set_xlim(system_disp_lims)
ax_system.set_ylim([0,box_length])


ax_coords.set_label("coords")

ax_coords.set_xlabel("$x_1$")
ax_coords.set_ylabel("$x_2$")
ax_coords.set_xticks([])
ax_coords.set_yticks([])
ax_coords.set_aspect(1)

coord_lim = [-disp_limit,disp_limit]
ax_coords.set_xlim(coord_lim)
ax_coords.set_ylim(coord_lim)
ax_coords.plot(coord_lim,[0,0],"k")
ax_coords.plot([0,0],coord_lim,"k")


#--------------------------------------#
## Plot initial data
marker_style = dict(markerfacecolor = marker_colour,
                    markeredgecolor = marker_colour,
                    marker = "o",
                    linestyle = "none",
                    markersize = marker_size)

[coord_marker] = ax_coords.plot(0,0,zorder=10,**marker_style)


lumped_mass_shapes = [Rectangle((box_one_centre+offset,0),box_length,box_length)
                     for offset in (0,box_distance)]

system_box_one = ax_system.add_patch(lumped_mass_shapes[0])
system_box_two = ax_system.add_patch(lumped_mass_shapes[1])


#--------------------------------------#
## Springs
SPRING_SEQUENCE = [+1,-1,-1,+1] #defines direction of components in spring section 

num_components = len(SPRING_SEQUENCE) #num components (lines) per section
num_spring_points = 3+num_components*num_spring_sections
box_mid_height = [box_length/2]*2
component_length = min_spring_length/2

def get_spring_points(spring_ends_x,spring_ends_y,length):
    angle = np.arccos((length/component_length -2)/(num_components*num_spring_sections))
    spring_coords = np.zeros([2,num_spring_points])
    spring_coords[:,0] = [spring_ends_x[0],spring_ends_y[0]]
    spring_coords[:,1] = [spring_ends_x[0] + component_length,spring_ends_y[0]]
    
    spring_coords[:,-1] = [spring_ends_x[1],spring_ends_y[1]]
    spring_coords[:,-2] = [spring_ends_x[1] - component_length,spring_ends_y[1]]
    
    x_diff = component_length*np.cos(angle)
    y_diff = component_length*np.sin(angle)
    
    for iPoint in range(num_spring_points-3):
        direction = SPRING_SEQUENCE[iPoint%num_components]
        next_point = spring_coords[:,iPoint+1] + [x_diff,direction*y_diff]
        spring_coords[:,iPoint+2] = next_point
    
    return spring_coords


end_left = max_disp_left
box_one_left = system_box_one.get_xy()[0]
box_one_right = box_one_left + box_length
box_two_left = system_box_two.get_xy()[0]
box_two_right = box_two_left + box_length
end_right = max_disp_right

spring_one_length = disp_limit + min_spring_length
spring_two_length = box_distance-box_length
spring_three_length = spring_one_length

spring_one_points = get_spring_points([end_left,box_one_left],box_mid_height,spring_one_length)
spring_two_points = get_spring_points([box_one_right,box_two_left],box_mid_height,spring_two_length)
spring_three_points = get_spring_points([box_two_right,end_right],box_mid_height,spring_three_length)

[spring_one] = ax_system.plot(spring_one_points[0],spring_one_points[1],"k")
[spring_two] = ax_system.plot(spring_two_points[0],spring_two_points[1],"k")
[spring_three] = ax_system.plot(spring_three_points[0],spring_three_points[1],"k")


def update_spring_points(spring,x):
    x_data = spring.get_xdata()
    y_data = spring.get_ydata()
    
    component_length = x_data[1] - x_data[0]
    x_data[0] = x[0]
    x_data[-1] = x[1]
    
    x_data[1:-1] = np.linspace(x[0] + component_length,x[1] - component_length,num_spring_points - 2)
    spring.set_xdata(x_data)
    
    component_dx = abs(x_data[2] - x_data[1])
    if component_dx < component_length:
        component_dy = np.sqrt(component_length**2 - component_dx**2)
    else:
        component_dy = 0.01 #stops the spring getting stuck fully extended

    y0 = y_data[0]
    for iPoint in range(2,num_spring_points-2):
        if iPoint%2 == 1:
            continue
        if y_data[iPoint] < y0:
            y_data[iPoint] = y0 - component_dy
        elif y_data[iPoint] > y0:
            y_data[iPoint] = y0+component_dy
    spring.set_ydata(y_data)

#--------------------------------------#
## Setup animation

    
def update_animation(frame_id,x,orbit):
    x_1 = x[0,frame_id]
    x_2 = x[1,frame_id]
    
    x_data = orbit.get_xdata()
    y_data = orbit.get_ydata()
    
    x_data = np.append(x_data,x_1)
    y_data = np.append(y_data,x_2)
    
    orbit.set_data(x_data,y_data)
    
    update_system(x_1,x_2)
    return orbit

def update_system(x_1,x_2):
    coord_marker.set_ydata([x_2])
    coord_marker.set_xdata([x_1])
    
    system_box_one.set_xy([x_1,0])
    system_box_two.set_xy([x_2+box_distance,0])
    
    x_spring_one = [end_left,x_1]
    x_spring_two = [x_1+box_length,x_2+box_distance]
    x_spring_three = [x_2+box_length+box_distance,end_right]
    
    update_spring_points(spring_one,x_spring_one)
    update_spring_points(spring_two,x_spring_two)
    update_spring_points(spring_three,x_spring_three)
    
def pause_animation(animation):
    animation.pause()
    animation.paused = True    
    return animation
    
def resume_animation(animation):
    animation.resume()
    animation.paused = False
    return animation
    
def toggle_animation(animation):
    if animation.paused:
        animation = resume_animation(animation)
    else:
        animation = pause_animation(animation)
    return animation
#--------------------------------------#

## Define interaction
def on_mouse_move(event):
    
    if "animation" in globals():
        #checks if an animation is currently playing
        return
    ax = event.inaxes
    if not ax:
        return
    if ax.get_label() != "coords":
        return
    
    x_1 = event.xdata
    x_2 = event.ydata

    update_system(x_1,x_2)
    event.canvas.draw()
    

    
def on_mouse_click(event):
    global animation
    
    if event.button is MouseButton.LEFT:
        ax = event.inaxes
        if not ax or ax.get_label() != "coords":
            if "animation" in globals():
                #checks if an animation exists
                animation = toggle_animation(animation)
            return
        
        
        #animate dynamics from selected initial condition
        if "animation" in globals():
            animation = pause_animation(animation)
        x_1 = event.xdata
        x_2 = event.ydata
        
        q = modal_displacement(t, [x_1,x_2])
        x = modal_to_physical(q)
        
        [orbit] = ax_coords.plot([],[],label="orbit")

        
        animation = ani.FuncAnimation(fig, partial(update_animation,orbit = orbit,x=x),
                                      frames=range(num_time_points),interval = 1000/fps)
        animation.paused = False
        event.canvas.draw()
        
    elif event.button is MouseButton.RIGHT:
        #clear old animations
        if "animation" in globals():
            #checks if an animation exists
            animation = pause_animation(animation)
            del animation
        
        lines = ax_coords.get_children()
        for line in lines:
            if line.get_label() == "orbit":
                line.remove()
        event.canvas.draw()
        
       
     
mouse_move_binding_id = plt.connect('motion_notify_event', on_mouse_move)
mouse_click_binding_id = plt.connect('button_press_event', on_mouse_click)