# 2 dof mass-spring system |--k--(m)--γk--(m)--k--|
#--------------------------------------#
γ = 1       #stiffness multiplier
β = 0.01    #stiffness proportional damping
m = 2       #mass
k = 2       #spring stiffness
forced_dof = 1  #dof where a sinusoidal force is applied

## Plot settings
plotted_dof = 1
num_frequency_points = 1e4
frequency_limit = [0.1,100]
log_frequency_axis = True

## Plot style
mode_one_colour = "r"
mode_two_colour = "g"
physical_colour = "b"

mode_line_style = "--"
physical_line_style = "-"
modal_over_physical = True #toggles which set of lines are on top

cursor_line_colour = "k"
cursor_line_width = 0.75

marker_size = 5
#--------------------------------------#
import matplotlib.pyplot as plt
import numpy as np
#--------------------------------------#
## Define system
NUM_DOF = 2
forced_index = forced_dof - 1

M = np.diag([m,m])
K = k*np.array([[γ+1,-γ],[-γ,γ+1]])


(eigenvalues,eigenvectors) = np.linalg.eig(np.linalg.inv(M)@K)
mode_order = eigenvalues.argsort()
eigenvalues = eigenvalues[mode_order]
eigenvectors = eigenvectors[:,mode_order]

modal_contribution = eigenvectors*eigenvectors[:,forced_index]

modal_mass = np.diag(eigenvectors.T*M*eigenvectors)
modal_stiffness = eigenvalues
modal_damping = β*modal_stiffness;

def get_mode_receptance(ω):
    s = 1j*ω
    Q_receptance = 1/(modal_mass*s**2 + modal_damping*s + eigenvalues)
    return Q_receptance


def modal_to_physical(q):
    y = eigenvectors@q
    return y

modal_to_physical = np.vectorize(modal_to_physical,signature="(q) -> (y)")
get_mode_receptance = np.vectorize(get_mode_receptance,signature="()->(2)")
#--------------------------------------#
## Setup figures
plt.close("all")

fig = plt.figure()
ax_nyquist = plt.subplot(1,2,2,label="nyquist")
ax_magnitude = plt.subplot(2,2,1,label="bode")
ax_phase = plt.subplot(2,2,3,label="bode")


# Nyquist plot
nyquist_title = lambda ω : f"ω = {ω:.2f}"

ax_nyquist.set_aspect(1)
ax_nyquist.set_xlabel(f"Real($\dot Y_{plotted_dof}/F)$")
ax_nyquist.set_ylabel(f"Imaginary($\dot Y_{plotted_dof}/F)$")
ax_nyquist.set_title(nyquist_title(0))

# Magnitude plot
ax_magnitude.set_ylabel("Magnitude (dB)")

ax_magnitude.set_xlabel("Frequency (rad/s)")
ax_magnitude.set_xlim(frequency_limit)
if log_frequency_axis: ax_magnitude.set_xscale("log")



# Phase plot
ax_phase.set_ylabel("Phase (deg)")

ax_phase.set_xlabel("Frequency (rad/s)")
ax_phase.set_xlim(frequency_limit)
if log_frequency_axis: ax_phase.set_xscale("log")
#--------------------------------------#
## Plot intial data
CURSOR_Z = 1
MODAL_Z = 2
PHYSICAL_Z = 3
if modal_over_physical: MODAL_Z,PHYSICAL_Z = PHYSICAL_Z,MODAL_Z

mode_style = dict(linestyle = mode_line_style,
                  zorder = MODAL_Z)
mode_one_style = dict(color = mode_one_colour,
                      **mode_style)
mode_two_style = dict(color = mode_two_colour,
                      **mode_style)

physical_style = dict(color = physical_colour,
                      linestyle = physical_line_style,
                      zorder = PHYSICAL_Z)

cursor_line_style = dict(color = cursor_line_colour,
                         linewidth = cursor_line_width,
                         zorder = CURSOR_Z)

marker_style = dict(marker = "o",
                    markersize = marker_size,
                    markeredgecolor = cursor_line_colour,
                    linewidth = cursor_line_width)


mode_one_label = "$Q_{1}$"
mode_two_label = "$Q_{2}$"
physical_label = f"$Y_{plotted_dof}$"

dof_index = plotted_dof - 1

num_frequency_points = round(num_frequency_points)
if log_frequency_axis:
    frequency = 10**np.linspace(*np.log10(frequency_limit),num_frequency_points)
else:
    frequency = np.linspace(*frequency_limit,num_frequency_points)
    
mode_receptance = get_mode_receptance(frequency)
receptance = modal_to_physical(mode_receptance)
receptance = receptance[:,dof_index]

def get_amplitude(z): return 20*np.log10(abs(z))
def get_angle(z):
    DEG_RANGE = 180
    theta = np.angle(z,deg = True)
    unwrap_index = theta > 0
    theta[unwrap_index] = theta[unwrap_index] - DEG_RANGE
    return theta

physical_amplitude = get_amplitude(receptance)
modal_amplitude = get_amplitude(mode_receptance)

physical_phase = get_angle(receptance)
modal_phase =  get_angle(mode_receptance)

# Nyquist plot
physical_mobility = 1j*frequency*receptance
real_mobility = np.real(physical_mobility)
imaginary_mobility = np.imag(physical_mobility)
ax_nyquist.plot(real_mobility,imaginary_mobility,**physical_style)

#modal_mobility = (1j*frequency*mode_receptance.T).T
#real_modal_receptance = np.real(modal_mobility)
#imaginary_modal_receptance = np.imag(modal_mobility)
#ax_nyquist.plot(real_modal_receptance[:,0],imaginary_modal_receptance[:,0],**mode_one_style)
#ax_nyquist.plot(real_modal_receptance[:,1],imaginary_modal_receptance[:,1],**mode_two_style)

[nyquist_marker] = ax_nyquist.plot([real_mobility[0],0],[imaginary_mobility[0],0],
                                   markevery = 2,**physical_style,**marker_style)

# Magnitude plot
ax_magnitude.plot(frequency,modal_amplitude[:,0],**mode_one_style, label = mode_one_label)
ax_magnitude.plot(frequency,modal_amplitude[:,1],**mode_two_style, label = mode_two_label)
ax_magnitude.plot(frequency,physical_amplitude,**physical_style, label = physical_label)

ax_magnitude.legend(loc="upper right")

magnitude_ylim = ax_magnitude.get_ylim()
ax_magnitude.set_ylim(magnitude_ylim)
[magnitude_cursor_line] = ax_magnitude.plot(2*[frequency[0]],magnitude_ylim,**cursor_line_style)

[magnitude_mode_one_marker] = ax_magnitude.plot(frequency[0],modal_amplitude[0,0],**mode_one_style,**marker_style)
[magnitude_mode_two_marker] = ax_magnitude.plot(frequency[0],modal_amplitude[0,1],**mode_two_style,**marker_style)
[magnitude_physical_marker] = ax_magnitude.plot(frequency[0],physical_amplitude[0],**physical_style,**marker_style)




# Phase plot
ax_phase.plot(frequency,modal_phase[:,0],**mode_one_style)
ax_phase.plot(frequency,modal_phase[:,1],**mode_two_style)
ax_phase.plot(frequency,physical_phase,**physical_style)

phase_ylim = ax_phase.get_ylim()
ax_phase.set_ylim(phase_ylim)
[phase_cursor_line] = ax_phase.plot(2*[frequency[0]],phase_ylim,**cursor_line_style)

[phase_mode_one_marker] = ax_phase.plot(frequency[0],modal_phase[0,0],**mode_one_style,**marker_style)
[phase_mode_two_marker] = ax_phase.plot(frequency[0],modal_phase[0,1],**mode_two_style,**marker_style)
[phase_physical_marker] = ax_phase.plot(frequency[0],physical_phase[0],**physical_style,**marker_style)
#--------------------------------------#
def update_markers(ω):
    mode_receptance = get_mode_receptance(ω)
    receptance = modal_to_physical(mode_receptance)
    receptance = receptance[dof_index]
    
    # Nyquist plot
    mobility = 1j*ω*receptance
    real_mobility = np.real(mobility)
    imaginary_mobility = np.imag(mobility)
    
    nyquist_marker.set_data([[real_mobility,0],[imaginary_mobility,0]])
    ax_nyquist.set_title(nyquist_title(ω))
    
    # Magnitude plot
    physical_amplitude = get_amplitude(receptance)
    modal_amplitude = get_amplitude(mode_receptance)

    magnitude_mode_one_marker.set_data([[ω],[modal_amplitude[0]]])
    magnitude_mode_two_marker.set_data([[ω],[modal_amplitude[1]]])
    magnitude_physical_marker.set_data([[ω]],[physical_amplitude])
    magnitude_cursor_line.set_xdata([ω,ω])
    
    # Phase plot
    [physical_phase] = get_angle([receptance])
    [modal_phase] =  get_angle([mode_receptance])

    
    phase_mode_one_marker.set_data([[ω],[modal_phase[0]]])
    phase_mode_two_marker.set_data([[ω],[modal_phase[1]]])
    phase_physical_marker.set_data([[ω]],[physical_phase])
    phase_cursor_line.set_xdata([ω,ω])

    
    
#--------------------------------------#
## Define interaction
def on_mouse_move(event):
    ax = event.inaxes
    if not ax: return
    if ax.get_label() != "bode": return
    
    mouse_x = event.xdata
    update_markers(mouse_x)
    event.canvas.draw()
    
def on_fig_close(event):
    plt.disconnect(mouse_move_binding_id)

mouse_move_binding_id = plt.connect('motion_notify_event', on_mouse_move)
plt.connect('close_event', on_fig_close)