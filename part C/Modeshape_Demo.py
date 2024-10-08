import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backend_bases import MouseButton
from matplotlib.patches import Rectangle




#--------------------------------------#
## System settings
m1 = 5 # Mass 1
m2 = 5 # Mass 2

k1 = 50  # Spring 1
k2 = 50  # Spring 2
k3 = 10  # Spring 3

## Animation settings
T = 10  # Time of simulation
FPS = 30   # Target FPS

## Plot settings
disp_limit = 2

box_length = 3
box_one_centre = 0

min_spring_length = 1
num_spring_sections = 3

## Plot style
marker_colour = "b"
marker_size = 6

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

[coord_marker] = ax_coords.plot(0,0,**marker_style)


lumped_mass_shapes = [Rectangle((box_one_centre+offset,0),box_length,box_length)
                     for offset in (0,box_distance)]

system_box_one = ax_system.add_patch(lumped_mass_shapes[0])
system_box_two = ax_system.add_patch(lumped_mass_shapes[1])


#--------------------------------------#
## Springs
REST_ANGLE = 2*np.pi/5
SPRING_SEQUENCE = [+1,-1,-1,+1]

num_spring_points = 3+4*num_spring_sections
box_mid_height = [box_length/2]*2


def get_spring_points(spring_ends_x,spring_ends_y,angle,length):
    component_length = length/(2+4*num_spring_sections*np.cos(angle))
    
    spring_coords = np.zeros([2,num_spring_points])
    spring_coords[:,0] = [spring_ends_x[0],spring_ends_y[0]]
    spring_coords[:,1] = [spring_ends_x[0] + component_length,spring_ends_y[0]]
    
    spring_coords[:,-1] = [spring_ends_x[1],spring_ends_y[1]]
    spring_coords[:,-2] = [spring_ends_x[1] - component_length,spring_ends_y[1]]
    
    x_diff = component_length*np.cos(angle)
    y_diff = component_length*np.sin(angle)
    
    for iPoint in range(num_spring_points-3):
        direction = SPRING_SEQUENCE[iPoint%4]
        next_point = spring_coords[:,iPoint+1] + [x_diff,direction*y_diff]
        spring_coords[:,iPoint+2] = next_point
    
    return spring_coords

def plot_springs():
    end_left = max_disp_left
    box_one_left = system_box_one.get_xy()[0]
    box_one_right = box_one_left + box_length
    box_two_left = system_box_two.get_xy()[0]
    box_two_right = box_two_left + box_length
    end_right = max_disp_right
    
    spring_one_length = disp_limit + min_spring_length
    spring_two_length = box_distance-box_length
    spring_three_length = spring_one_length

    
    spring_one_points = get_spring_points([end_left,box_one_left],box_mid_height,REST_ANGLE,spring_one_length)
    spring_two_points = get_spring_points([box_one_right,box_two_left],box_mid_height,REST_ANGLE,spring_two_length)
    spring_three_points = get_spring_points([box_two_right,end_right],box_mid_height,REST_ANGLE,spring_three_length)
    
    [spring_one] = ax_system.plot(spring_one_points[0],spring_one_points[1],"k")
    [spring_two] = ax_system.plot(spring_two_points[0],spring_two_points[1],"k")
    [spring_three] = ax_system.plot(spring_three_points[0],spring_three_points[1],"k")
        
    return spring_one,spring_two,spring_three


def update_spring_points(spring,x):
    x_data = spring.get_xdata()
    y_data = spring.get_ydata()
    
    component_length = x_data[1] - x_data[0]
    x_data[0] = x[0]
    x_data[-1] = x[1]
    
    num_free_points = num_spring_points - 4
    x_data[1:-1] = np.linspace(x[0] + component_length,x[1] - component_length,num_free_points+2)
    spring.set_xdata(x_data)
    
    component_dx = abs(x_data[2] - x_data[1])
    if component_dx < component_length:
        component_dy = np.sqrt(component_length**2 - component_dx**2)
    else:
        component_dy = 0.01

    y0 = y_data[0]
    for iPoint in range(2,num_spring_points-2):
        if iPoint%2 == 1:
            continue
        if y_data[iPoint] < y0:
            y_data[iPoint] = y0 - component_dy
        elif y_data[iPoint] > y0:
            y_data[iPoint] = y0+component_dy
    spring.set_ydata(y_data)

spring_one,spring_two,spring_three = plot_springs()
#--------------------------------------#
## Define interaction
def on_mouse_move(event):
    ax = event.inaxes
    if not ax:
        return
    if ax.get_label() != "coords":
        return
    
    x_1 = event.xdata
    x_2 = event.ydata

    coord_marker.set_ydata([x_2])
    coord_marker.set_xdata([x_1])
    
    system_box_one.set_xy([x_1,0])
    system_box_two.set_xy([x_2+box_distance,0])
    
    x_spring_one = [max_disp_left,x_1]
    x_spring_two = [x_1+box_length,x_2+box_distance]
    x_spring_three = [x_2+box_length+box_distance,max_disp_right]
    
    update_spring_points(spring_one,x_spring_one)
    update_spring_points(spring_two,x_spring_two)
    update_spring_points(spring_three,x_spring_three)
    event.canvas.draw()
 

mouse_move_binding_id = plt.connect('motion_notify_event', on_mouse_move)
