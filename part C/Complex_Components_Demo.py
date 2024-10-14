#--------------------------------------#
############### CONTROLS ############### 
# Left click to toggle animation
#--------------------------------------#
## User settings
period = 1/2
amplitude = 1
phase_deg =  0 #in degrees

## Plot settings
num_periods = 2
num_time_points = 1000
complex_ax_view = dict(elev = 15, azim = -140, roll = 0)

## Animation settings
animation_time = 10 #number of seconds for complete animation

## Plot style
marker_size = 6
line_width = 1

signal_colour = "k"
plus_colour = "b"
minus_colour = "r" 

y_plus_label = "$1/2 e^{+i(\omega_2t + \phi)}$"
y_minus_label = "$1/2 e^{-i(\omega_2t + \phi)}$"

#--------------------------------------#
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
from functools import partial
#--------------------------------------#
## Governing equations
frequency = 2*np.pi/period
t = np.linspace(0,period*num_periods,num_time_points)

phase = phase_deg*np.pi/180
angle = frequency*t + phase
y = amplitude*np.cos(angle)

y_plus = amplitude/2 * np.exp(1j*angle)
y_plus_real = np.real(y_plus)
y_plus_imag = np.imag(y_plus)

y_minus = amplitude/2 * np.exp(-1j*angle)
y_minus_real = np.real(y_minus)
y_minus_imag = np.imag(y_minus)
#--------------------------------------#

## Setup figures
fig = plt.figure(1)
ax_signal = plt.subplot(121)
ax_complex = plt.subplot(122,projection = "3d")

time_span = [0,period*num_periods]
amplitude_span = [-amplitude,amplitude]

ax_signal.set_xlabel("t")
ax_signal.set_ylabel("y")
ax_signal.set_xlim([0,period*num_periods])
ax_signal.set_ylim([-amplitude,amplitude])


amplitude_half_span = [-amplitude/2,amplitude/2]

ax_complex.view_init(**complex_ax_view)
ax_complex.set_xlabel("t")
ax_complex.set_ylabel("real(y)")
ax_complex.set_zlabel("imag(y)")
ax_complex.set_xlim(time_span)
ax_complex.set_ylim(amplitude_half_span)
ax_complex.set_zlim(amplitude_half_span)



#--------------------------------------#
## Plot initial data
line_style = dict(linewidth = line_width)
marker_style = dict(markersize = marker_size,
                    marker = "o",
                    linewidth = line_width)

signal_style = dict(color = signal_colour)
plus_style = dict(color = plus_colour)
minus_style = dict(color = minus_colour)


ax_signal.plot(t,y,**signal_style,**line_style)
[signal_marker] = ax_signal.plot(t[0],y[0],**signal_style,**marker_style)

ax_complex.plot(t,y_minus_real,y_minus_imag,**minus_style,**line_style,label=y_minus_label)
[minus_marker] = ax_complex.plot(t[0],y_minus_real[0],y_minus_imag[0],**minus_style,**marker_style)

ax_complex.plot(t,y_plus_real,y_plus_imag,**plus_style,**line_style,label=y_plus_label)
[plus_marker] = ax_complex.plot(t[0],y_plus_real[0],y_plus_imag[0],**plus_style,**marker_style)

plt.legend()

SIGNAL_INDEX, PLUS_INDEX, MINUS_INDEX = [0,1,2]
markers = [signal_marker,plus_marker,minus_marker]
#--------------------------------------#
## Setup animation
frame_interval = animation_time/num_time_points

def update_animation(frame_id,markers):
    t_frame = t[frame_id]
    y_frame = y[frame_id]
    y_plus_real_frame = y_plus_real[frame_id]
    y_plus_imag_frame = y_plus_imag[frame_id]
    y_minus_real_frame = y_minus_real[frame_id]
    y_minus_imag_frame = y_minus_imag[frame_id]
    
    markers[SIGNAL_INDEX].set_xdata([t_frame])
    markers[SIGNAL_INDEX].set_ydata([y_frame])
    
    markers[PLUS_INDEX].set_data_3d([t_frame],[y_plus_real_frame],[y_plus_imag_frame])
    markers[MINUS_INDEX].set_data_3d([t_frame],[y_minus_real_frame],[y_minus_imag_frame])
    return markers




#--------------------------------------#
## Define interaction
class Animation_interaction:
    animation = []
    animation_paused = []
    def animation_exists(self):
        return self.animation and self.animation.event_source

    def pause_animation(self):
        if not self.animation_exists(): return

        self.animation.pause()
        self.animation_paused = True    
        
    def resume_animation(self):
        if not self.animation_exists(): return
        
        self.animation.resume()
        self.animation_paused = False
        
    def toggle_animation(self):
        if not self.animation_exists(): return
        
        if self.animation_paused: self.resume_animation()
        else: self.pause_animation()
    
        
    def on_mouse_click(self,event):
        if self.animation_exists():
            self.toggle_animation()
            return
        
        
        self.animation = ani.FuncAnimation(fig, partial(update_animation,markers = markers),
                                      frames=range(num_time_points),interval = 1000*frame_interval,
                                      repeat = True,blit="True")
        self.animation_paused = False
        event.canvas.draw()



def on_fig_close(event):
    plt.disconnect(mouse_click_binding_id)

mouse_callback = Animation_interaction()
mouse_click_binding_id = plt.connect('button_press_event', mouse_callback.on_mouse_click)
plt.connect('close_event', on_fig_close)