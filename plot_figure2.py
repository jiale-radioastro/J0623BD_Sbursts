import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

class burst:
    def __init__(self,time,time_sigma,freq,freq_span,drift,drift_sigma,snr,Iflux,Vflux):
        self.time=time
        self.time_sigma=time_sigma
        self.freq=freq
        self.freq_span=freq_span
        self.drift=drift
        self.drift_sigma=drift_sigma
        self.snr=snr
        self.Iflux=Iflux
        self.Vflux=Vflux

def trebin(data, ttbin):
    """
    Rebin along time axis.
    data: 1D (time) or 2D (time, freq)
    """
    data = np.asarray(data)
    if ttbin == 1:
        return data

    n = data.shape[0] // ttbin  # number of complete bins
    if data.ndim == 2:
        trimmed = data[:n * ttbin, :]
        reshaped = trimmed.reshape(n, ttbin, data.shape[1])
        return np.mean(reshaped, axis=1)
    elif data.ndim == 1:
        trimmed = data[:n * ttbin]
        reshaped = trimmed.reshape(n, ttbin)
        return np.mean(reshaped, axis=1)
    else:
        raise ValueError("Input must be 1D or 2D")


def frebin(data, ffbin):
    """
    Rebin along frequency axis.
    data: 1D (freq) or 2D (time, freq)
    """
    data = np.asarray(data)
    if ffbin == 1:
        return data
    
    n = data.shape[-1] // ffbin  # number of complete bins
    if data.ndim == 2:
        trimmed = data[:, :n * ffbin]
        reshaped = trimmed.reshape(data.shape[0], n, ffbin)
        return np.mean(reshaped, axis=2)
    elif data.ndim == 1:
        trimmed = data[:n * ffbin]
        reshaped = trimmed.reshape(n, ffbin)
        return np.mean(reshaped, axis=1)
    else:
        raise ValueError("Input must be 1D or 2D")


file_dir='data/'

bursts1_list=np.load(file_dir+'Sbursts20241031_list.npy',allow_pickle=True)
burst1_time=np.array([burst0.time for burst0 in bursts1_list])
burst1_tsigma=np.abs(np.array([burst0.time_sigma for burst0 in bursts1_list]))
burst1_twidth=np.abs(np.array([burst0.time_sigma*2*np.sqrt(2*np.log(2)) for burst0 in bursts1_list]))
burst1_freq=np.array([burst0.freq for burst0 in bursts1_list])
burst1_fwidth=np.array([np.diff(burst0.freq_span)[0] for burst0 in bursts1_list])
burst1_drift=np.array([burst0.drift for burst0 in bursts1_list])
burst1_drift_err=np.array([burst0.drift_sigma for burst0 in bursts1_list])
burst1_snr=np.array([burst0.snr for burst0 in bursts1_list])
burst1_Iflux=np.array([burst0.Iflux for burst0 in bursts1_list])
burst1_Vflux=np.array([burst0.Vflux for burst0 in bursts1_list])
burst1_Lflux=(burst1_Iflux+burst1_Vflux)/2
burst1_dur=burst1_fwidth/burst1_drift

bursts2_list=np.concatenate([np.load(file_dir+'Sbursts20251018_list1.npy',allow_pickle=True),\
                             np.load(file_dir+'Sbursts20251018_list2.npy',allow_pickle=True),\
                             np.load(file_dir+'Sbursts20251018_list3.npy',allow_pickle=True)])
burst2_time=np.array([burst0.time for burst0 in bursts2_list])
burst2_tsigma=np.abs(np.array([burst0.time_sigma for burst0 in bursts2_list]))
burst2_twidth=np.abs(np.array([burst0.time_sigma*2*np.sqrt(2*np.log(2)) for burst0 in bursts2_list]))
burst2_freq=np.array([burst0.freq for burst0 in bursts2_list])
burst2_fwidth=np.array([np.diff(burst0.freq_span)[0] for burst0 in bursts2_list])
burst2_drift=np.array([burst0.drift for burst0 in bursts2_list])
burst2_drift_err=np.array([burst0.drift_sigma for burst0 in bursts2_list])
burst2_snr=np.array([burst0.snr for burst0 in bursts2_list])
burst2_Iflux=np.array([burst0.Iflux for burst0 in bursts2_list])
burst2_Vflux=np.array([burst0.Vflux for burst0 in bursts2_list])
burst2_Lflux=(burst2_Iflux+burst2_Vflux)/2
burst2_dur=burst2_fwidth/burst2_drift

npy_file = file_dir +'Sbursts20241031_dspec.npy'

data=np.load(npy_file,allow_pickle=True).item()
V_dspec1=data['V_dspec']
I_dspec1=data['I_dspec']
L_dspec1=(I_dspec1+V_dspec1)/2
time1=data['timelist']
freq1=data['freqlist']
mjd1=data['mjd0']

L_std=[]
for i in range(np.shape(V_dspec1)[1]):
    L_std.append(np.nanstd(L_dspec1[:,i]))
L_std=np.array(L_std)

mask1=np.ones(np.shape(V_dspec1)[1])
for index in range(len(mask1)):
    if L_std[index]>170:
        mask1[index]=np.nan
freq1 = trebin(freq1,2)
dfreq1 = freq1[1]-freq1[0]
L_dspec1 = frebin(trebin(L_dspec1,4),2)
mask1 = trebin(mask1,2)

npy_file = file_dir +'Sbursts20251018_dspec1.npy'

data=np.load(npy_file,allow_pickle=True).item()
V_dspec2=data['V_dspec']
I_dspec2=data['I_dspec']
L_dspec2=(I_dspec2+V_dspec2)/2
time2=data['timelist']
freq2=data['freqlist']
mjd2=data['mjd0']

L_std=[]
for i in range(np.shape(V_dspec2)[1]):
    L_std.append(np.nanstd(L_dspec2[:,i]))
L_std=np.array(L_std)

mask2=np.ones(np.shape(V_dspec2)[1])
for index in range(len(mask2)):
    if L_std[index]>170:
        mask2[index]=np.nan
freq2 = trebin(freq2,2)
dfreq2 = freq2[1]-freq2[0]
L_dspec2 = frebin(trebin(L_dspec2,4),2)
mask2 = trebin(mask2,2)

fig = plt.figure(figsize=(10,8),dpi=200)
colormap=cm.viridis
norm = Normalize(vmin=0, vmax=100)

fig = plt.figure(figsize=(8,8),dpi=200)
ax1 = fig.add_axes([0.1,0.82,0.85,0.18])
ax2 = fig.add_axes([0.1,0.6,0.85,0.18])
ax3 = fig.add_axes([0.1,0.32,0.85,0.18])
ax4 = fig.add_axes([0.1,0.1,0.85,0.18])

for i in range(len(burst1_time)):
    line_time0=burst1_time[i]-burst1_dur[i]/2
    line_time1=burst1_time[i]+burst1_dur[i]/2
    line_freq0=bursts1_list[i].freq_span[0]
    line_freq1=bursts1_list[i].freq_span[1]
    color=colormap(norm(burst1_Lflux[i]))
    ax1.plot([line_time0,line_time1],[line_freq0,line_freq1],color=color, lw=0.7)
for freqi in range(len(freq1)):
    if np.isnan(mask1[freqi]):
        ax1.fill_between([-0.2,0.05],[freq1[freqi]-dfreq2/2,freq1[freqi]-dfreq1/2],\
                         [freq1[freqi]+dfreq2/2,freq1[freqi]+dfreq1/2],facecolor='red',alpha=1)
ax1.set_xlim([-0.2,10.2])
ax1.set_ylim([1000,1500])
ax1.set_ylabel("Frequency [MHz]")
ax1.fill_between([3,4],[1280,1280],[1480,1480],facecolor='darkorange',alpha=0.4,zorder=10)
ax1.plot([3,-0.2],[1280,900],linestyle=':',color='darkorange',clip_on=False)
ax1.plot([4,10.2],[1280,900],linestyle=':',color='darkorange',clip_on=False)
ax1.text(0.02,0.07,'a',color='k',transform=ax1.transAxes,fontsize=12,fontweight='bold')

ax2.imshow((L_dspec1*mask1).T,aspect='auto',origin='lower',cmap='viridis',\
           extent=[time1[0],time1[-1],freq1[0]-dfreq1/2,freq1[-1]+dfreq1/2],vmin=0,vmax=100,interpolation="none")
ax2.set_xlim([3,4])
ax2.set_ylim([1280,1480])
ax2.set_ylabel('Frequency [MHz]')
ax2.set_xlabel('Time [s] since MJD '+str(mjd1))
for freqi in range(len(freq1)):
    if np.isnan(mask1[freqi]):
        ax2.fill_between([3,3.025],[freq1[freqi]-dfreq1/2,freq1[freqi]-dfreq1/2],\
                         [freq1[freqi]+dfreq1/2,freq1[freqi]+dfreq1/2],facecolor='red',alpha=1)
ax2.text(0.02,0.07,'b',color='white',transform=ax2.transAxes,fontsize=12,fontweight='bold')

windows = [(-5, 35),(75, 95),(305, 335)]
gap = 6
offsets = []
gap_edges = []
cursor = 0
for j, (lo, hi) in enumerate(windows):
    offsets.append(cursor - lo)
    cursor += hi - lo

    if j < len(windows) - 1:
        gap_start = cursor
        gap_end = cursor + gap
        gap_edges.append((gap_start, gap_end))
        cursor += gap

def compressed_time(x):
    x = np.asarray(x, dtype=float)
    y = np.full_like(x, np.nan)

    for (lo, hi), offset in zip(windows, offsets):
        mask = (x >= lo) & (x <= hi)
        y[mask] = x[mask] + offset

    return y
for i in range(len(burst2_time)):
    line_time0=burst2_time[i]-burst2_dur[i]/2
    line_time1=burst2_time[i]+burst2_dur[i]/2
    line_freq0=bursts2_list[i].freq_span[0]
    line_freq1=bursts2_list[i].freq_span[1]
    color = colormap(norm(burst2_Lflux[i]))
    x = compressed_time([line_time0, line_time1])
    ax3.plot(x, [line_freq0, line_freq1], color=color, lw=0.7)
tick_values = [0, 10, 20, 30, 80, 90, 310, 320, 330]
tick_positions = compressed_time(tick_values)
ax3.set_xticks(tick_positions)
ax3.set_xticklabels(tick_values)
xmin = compressed_time([windows[0][0]])[0]
xmax = compressed_time([windows[-1][1]])[0]
ax3.set_xlim(xmin, xmax)
ax3.set_ylabel("Frequency [MHz]")
trans = ax3.get_xaxis_transform()
for gap_start, gap_end in gap_edges:
    gap_mid = 0.5 * (gap_start + gap_end)
    mark_dx = 0.5
    mark_sep = 0.6
    left_mark_x = gap_mid - mark_sep / 2
    right_mark_x = gap_mid + mark_sep / 2
    for yaxis in [0, 1]:
        ax3.plot(
            [left_mark_x, right_mark_x],
            [yaxis, yaxis],
            transform=trans,
            color=ax3.get_facecolor(),
            lw=4,
            clip_on=False,
            zorder=5
        )
        ax3.plot(
            [left_mark_x - mark_dx / 2, left_mark_x + mark_dx / 2],
            [yaxis - 0.035, yaxis + 0.035],
            transform=trans,
            color="black",
            lw=1,
            clip_on=False,
            zorder=6
        )
        ax3.plot(
            [right_mark_x - mark_dx / 2, right_mark_x + mark_dx / 2],
            [yaxis - 0.035, yaxis + 0.035],
            transform=trans,
            color="black",
            lw=1,
            clip_on=False,
            zorder=6
        )
for freqi in range(len(freq2)):
    if np.isnan(mask2[freqi]):
        ax3.fill_between(compressed_time([-5,-3]),[freq2[freqi]-dfreq2/2,freq2[freqi]-dfreq2/2],\
                         [freq2[freqi]+dfreq2/2,freq2[freqi]+dfreq2/2],facecolor='red',alpha=1)
ax3.set_ylim([1000,1500])
ax3.fill_between(compressed_time([13.1,14.1]),[1000,1000],[1200,1200],facecolor='darkorange',alpha=0.4,zorder=10)
ax3.plot(compressed_time([13.1,-5]),[1000,900],linestyle=':',color='darkorange',clip_on=False)
ax3.plot(compressed_time([14.1,335]),[1000,900],linestyle=':',color='darkorange',clip_on=False)
ax3.text(0.02,0.07,'c',color='k',transform=ax3.transAxes,fontsize=12,fontweight='bold')

im2=ax4.imshow((L_dspec2*mask2).T,aspect='auto',origin='lower',cmap='viridis',\
           extent=[time2[0],time2[-1],freq2[0]-dfreq2/2,freq2[-1]+dfreq2/2],vmin=0,vmax=100,interpolation="none")
ax4.set_xlim([13.1,14.1])
ax4.set_ylim([1000,1200])
ax4.set_ylabel('Frequency [MHz]')
ax4.set_xlabel('Time [s] since MJD '+str(mjd2))
for freqi in range(len(freq2)):
    if np.isnan(mask2[freqi]):
        ax4.fill_between([13.1,13.125],[freq2[freqi]-dfreq2/2,freq2[freqi]-dfreq2/2],\
                         [freq2[freqi]+dfreq2/2,freq2[freqi]+dfreq2/2],facecolor='red',alpha=1)
ax4.text(0.02,0.07,'d',color='white',transform=ax4.transAxes,fontsize=12,fontweight='bold')

cax = fig.add_axes([0.82, 0.88, 0.11, 0.01])  
cbar = fig.colorbar(im2, cax=cax, orientation='horizontal',label='LCP [mJy]')

cax = fig.add_axes([0.82, 0.38, 0.11, 0.01])  
cbar = fig.colorbar(im2, cax=cax, orientation='horizontal',label='LCP [mJy]')

plt.show()
fig.savefig('figures/figure2.pdf',dpi=200,bbox_inches='tight')