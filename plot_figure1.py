import numpy as np
import matplotlib.pyplot as plt

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
        return np.nanmean(reshaped, axis=1)
    elif data.ndim == 1:
        trimmed = data[:n * ttbin]
        reshaped = trimmed.reshape(n, ttbin)
        return np.nanmean(reshaped, axis=1)
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
        return np.nanmean(reshaped, axis=2)
    elif data.ndim == 1:
        trimmed = data[:n * ffbin]
        reshaped = trimmed.reshape(n, ffbin)
        return np.nanmean(reshaped, axis=1)
    else:
        raise ValueError("Input must be 1D or 2D")

target='WISE_J0623-0456'
project='PT2025_0037'
beam=1
date_list=['20241027','20241029','20241031','20241102','20250506','20250508','20250509',\
           '20251018','20251022','20251023','20251026']
file_dir='data/'

mjd0_list = []
bjd0_list = []
bjds_list = []
times_list = []
freqs_list = []
dspec_list = []


for date in date_list:
    npy_name = file_dir+'dspec_filter_'+date+'_beam_'+str(beam).zfill(2)+'.npy'
    data=np.load(npy_name,allow_pickle=True).item()
    mjd0_list.append(data['mjd0'])
    bjd0_list.append(data['bjd0'])
    times_list.append(data['timelist'])
    bjds_list.append(data['bjd0']+data['timelist']/24)
    freqs_list.append(data['freqlist'])
    dspec_list.append(data['V_filter'])


period = 1.919475/24
bjd0_ref = 2460031.047
offset_cycle = [1,1,1,1,1,1,1,1,0,1,0]
bjds_sbursts = [2460614.328798336,2460966.4073807923]

fig = plt.figure(figsize=(10,6),dpi=200)

gs = fig.add_gridspec(6, 2, hspace=0.1, wspace=0.05)

for datei in range(11):
    mjd0 = mjd0_list[datei]
    bjd0 = bjd0_list[datei]
    bjds = bjds_list[datei]
    freqs = freqs_list[datei]
    times = times_list[datei]
    dspec = dspec_list[datei]
    cyclei = (bjd0-bjd0_ref)//period + offset_cycle[datei]
    phases = (bjds - bjd0_ref - cyclei*period)/period

    row = datei % 6
    col = datei // 6
    ax = fig.add_subplot(gs[row, col])
    im1 = ax.imshow(
        trebin(dspec,10).T,
        origin='lower',
        vmin=0, vmax=6,
        cmap='viridis',
        aspect='auto',
        extent=[phases[0], phases[-1], freqs[0], freqs[-1]]
    )
    ax.set_xlim([-0.8, 1.8])
    ax.set_ylim([freqs[0], freqs[-1]])
    if datei not in [5,10]:
        ax.set_xticks([])
    if datei>5:
        ax.set_yticks([])
    if datei in [5,10]:
        ax.set_xlabel('Phase')
    if datei==0:
        ax.set_ylabel('Freq [MHz]')
    if datei in [2,3,6,7,9]:
        ax.text(1.2, 1425, 'MJD ' + str(round(mjd0,1)),color='brown',fontsize=9)
    else:
        ax.text(-0.75, 1425, 'MJD ' + str(round(mjd0,1)),color='brown',fontsize=9)

    if datei == 2:
        ax.arrow((bjds_sbursts[0]- bjd0_ref - cyclei*period)/period,freqs[0],0,30,head_width=0.02,head_length=15,color='orange')

    if datei == 7:
        ax.arrow((bjds_sbursts[1]- bjd0_ref - cyclei*period)/period,freqs[0],0,30,head_width=0.02,head_length=15,color='orange')
    

cax2 = fig.add_axes([0.62, 0.15, 0.2, 0.017]) 
cbar2 = fig.colorbar(im1, cax=cax2, orientation='horizontal',label='Stokes V [mJy]')

fig.savefig('figures/figure1.pdf',dpi=200,bbox_inches='tight')