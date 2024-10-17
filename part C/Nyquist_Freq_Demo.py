## system settings
response_freq = 1;
sampling_freq = 1;

## plot settings
num_periods = 10
num_time_points = 100

signal_colour = "b"
sample_colour = "r"

guide_line_thickess = 0.4
#--------------------------------------#
import matplotlib.pyplot as plt
import numpy as np
#--------------------------------------#
## system
AMPLITUDE = 1
signal_period = 2*np.pi/response_freq
t_signal = np.linspace(0,signal_period*num_periods,num_time_points*num_periods)

sample_period = 2*np.pi/sampling_freq
t_sample = np.arange(signal_period*num_periods+sample_period,step=sample_period,)


N = np.ceil(2*response_freq/sampling_freq - 1)
x = AMPLITUDE*np.cos(response_freq*t_signal)
x_sample = AMPLITUDE*np.cos(response_freq*t_sample)
y = AMPLITUDE*np.cos((response_freq - N*sampling_freq)*t_signal)
#--------------------------------------#
## Setup figures
AMP_LIM_BUFFER = 1.1
fig = plt.figure(1)
ax = plt.subplot(111)

ax.set_xlabel("t")
ax.set_ylabel("x(t)")

ax.set_xlim([0,sample_period*num_periods])
y_limits = [-AMPLITUDE*AMP_LIM_BUFFER,AMPLITUDE*AMP_LIM_BUFFER]
ax.set_ylim(y_limits)

#--------------------------------------#
## Plot data
guide_line_style = dict(linewidth = guide_line_thickess,
                        color = "k")
num_sample = len(t_sample)
for t in t_sample: ax.plot(2*[t],y_limits,**guide_line_style)
ax.plot([0,sample_period*num_periods],[0,0],**guide_line_style)

marker_style = dict(marker = "o",
                    markeredgecolor = sample_colour,
                    markerfacecolor = signal_colour,
                    linestyle = "none")

ax.plot(t_signal,y,"--",color = sample_colour)
ax.plot(t_signal,x,color = signal_colour)
ax.plot(t_sample,x_sample,**marker_style)

