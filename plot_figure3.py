import numpy as np
import matplotlib.pyplot as plt
import string

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

file_dir='/Volumes/jiale_disk1/projects/PT2025_0037/20251018/'
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

fig=plt.figure(figsize=(9,5.5),dpi=200)
labels = list(string.ascii_lowercase)  # ['a','b','c',...]
for i in range(6):
    ax = plt.subplot(2, 3, i+1)

    ax.text(
        0.05, 0.95, f'{labels[i]}',
        transform=ax.transAxes,
        fontsize=12,
        fontweight='bold',
        va='top'
    )

    if i == 0:
        counts, bins = np.histogram(burst1_drift/1e3, bins=15, range=[0,10])
        ax.step(bins[:-1], counts, where='post', color='darkorange')
        ax.plot([np.median(burst1_drift/1e3),np.median(burst1_drift/1e3)],[0,120],':',color='darkorange')
        counts, bins = np.histogram(burst2_drift/1e3, bins=15, range=[0,10])
        ax.step(bins[:-1], counts, where='post', color='grey')
        ax.plot([np.median(burst2_drift/1e3),np.median(burst2_drift/1e3)],[0,120],':',color='grey')
        ax.set_ylim([0,120])
        ax.set_xticks([0,2,4,6,8])
        ax.set_xlabel(r'Drift rate [GHz s$^{-1}$]')
        ax.set_ylabel('Counts')

    elif i == 1:
        counts, bins = np.histogram(burst1_twidth*burst1_drift,bins=15,range=[-2,28])
        ax.step(bins[:-1], counts, where='post', color='darkorange')
        ax.plot([np.median(burst1_twidth*burst1_drift),np.median(burst1_twidth*burst1_drift)],[0,80],':',color='darkorange')
        counts, bins = np.histogram(burst2_twidth*burst2_drift,bins=15,range=[-2,28])
        ax.step(bins[:-1], counts, where='post', color='grey')
        ax.plot([np.median(burst2_twidth*burst2_drift),np.median(burst2_twidth*burst2_drift)],[0,80],':',color='grey')
        ax.set_xlabel('Instantaneous bandwidth [MHz]')
        ax.set_ylim([0,80])

    elif i == 2:
        counts, bins = np.histogram(burst1_dur*1e3,bins=15,range=[0,30])
        ax.step(bins[:-1], counts, where='post', color='darkorange')
        ax.plot([np.median(burst1_dur*1e3),np.median(burst1_dur*1e3)],[0,80],':',color='darkorange')
        counts, bins = np.histogram(burst2_dur*1e3,bins=15,range=[0,30])
        ax.step(bins[:-1], counts, where='post', color='grey')
        ax.plot([np.median(burst2_dur*1e3),np.median(burst2_dur*1e3)],[0,80],':',color='grey')
        ax.set_xlabel('Duration [ms]')
        ax.set_ylim([0,80])

    elif i == 3:
        counts, bins = np.histogram(burst1_snr,bins=15,range=[2,20])
        ax.step(bins[:-1], counts, where='post', color='darkorange')
        ax.plot([np.median(burst1_snr),np.median(burst1_snr)],[0,150],':',color='darkorange')
        counts, bins = np.histogram(burst2_snr,bins=15,range=[2,20])
        ax.step(bins[:-1], counts, where='post', color='grey')
        ax.plot([np.median(burst2_snr),np.median(burst2_snr)],[0,150],':',color='grey')
        ax.set_xlabel('SNR')
        ax.set_ylim([0,150])
        ax.set_ylabel('Counts')

    elif i == 4:
        counts, bins = np.histogram(burst1_Lflux,bins=15,range=[0,200])
        ax.step(bins[:-1], counts, where='post', color='darkorange')
        ax.plot([np.median(burst1_Lflux),np.median(burst1_Lflux)],[0,100],':',color='darkorange')
        counts, bins = np.histogram(burst2_Lflux,bins=15,range=[0,200])
        ax.step(bins[:-1], counts, where='post', color='grey')
        ax.plot([np.median(burst2_Lflux),np.median(burst2_Lflux)],[0,100],':',color='grey')
        ax.set_xlabel('LCP flux density [mJy]')
        ax.set_ylim([0,100])

    elif i == 5:
        counts, bins = np.histogram(burst1_Vflux/burst1_Iflux,bins=15,range=[0,3.3])
        ax.step(bins[:-1], counts, where='post', color='darkorange', label='MJD 60613.8')
        ax.plot([np.median(burst1_Vflux/burst1_Iflux),np.median(burst1_Vflux/burst1_Iflux)],[0,100],':',color='darkorange')
        counts, bins = np.histogram(burst2_Vflux/burst2_Iflux,bins=15,range=[0,3.3])
        ax.step(bins[:-1], counts, where='post', color='grey', label='MJD 60965.8')
        ax.plot([np.median(burst2_Vflux/burst2_Iflux),np.median(burst2_Vflux/burst2_Iflux)],[0,100],':',color='grey')
        ax.set_xlabel('Circular polarization')
        ax.set_ylim([0,100])
        ax.legend()

plt.tight_layout()
plt.show()

fig.savefig('figures/figure3.pdf',dpi=200,bbox_inches='tight')