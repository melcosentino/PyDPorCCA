import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyhydrophone as pyhy
import soundfile
from pyporcc import click_detector
from pyporcc import porcc
from scipy import signal

model = 'ST300HF'
name = 'SoundTrap'
serial_number = 738496579
soundtrap = pyhy.soundtrap.SoundTrap(name=name, model=model, serial_number=serial_number)

MaxLenFile = 200000
LongFilt = 0.00001
LongFilt2 = 0.000001
ShortFilt = 0.1
DetThres = 10
PreSam = 40
PostSam = 40
MaxLenClick = 1024
MinSep = 100
Fs = 576000
MinFrq = 100000 / Fs
MaxFrq = 150000 / Fs
MainFolder = 'D:/JOANNA/SoundTrap'
PathDetSave = 'C:/Mel/CPODvsDPorCCA/TruncatedData'
models_config_path = 'C:/Mel/PythonEnvironments/DPorCCA/pyporcc/models/log_models.ini'
LQ = 0.6
HQ = 0.999999
Fs = 576000
# Signal, Fs = soundfile.read(FileToOpen, start=int(Start), stop=int(End))
pfilter = click_detector.Filter(filter_name='butter', filter_type='bandpass', order=4, frequencies=[MinFrq, MaxFrq])
dfilter = click_detector.Filter(filter_name='butter', filter_type='high', order=4, frequencies=20000)
classifier = porcc.PorCC(load_type='manual', config_file='default')
# update the thresholds
classifier.th1 = HQ
classifier.th2 = LQ
detector = click_detector.ClickDetector(hydrophone=soundtrap, long_filt=LongFilt, long_filt2=LongFilt2,
                                        short_filt=ShortFilt, threshold=DetThres, min_separation=MinSep,
                                        max_length=MaxLenClick, pre_samples=PreSam, post_samples=PostSam,
                                        prefilter=pfilter, dfilter=dfilter, save_max=MaxLenFile,
                                        save_folder=MainFolder, convert=True, click_model_path=None,
                                        classifier=classifier)
Clicks = classifier.classify_matrix(Clicks)
Clicks = Clicks.drop(Clicks[Clicks.pyPorCC == 3][:].index)
Clicks.to_csv('C:/Mel/CPODvsDPorCCA/CorrecteddB/Clicks.csv')
ClicksWOwave = Clicks.drop(['wave'], axis=1)
ClicksWOwave.to_csv('C:/Mel/CPODvsDPorCCA/CorrecteddB/ClicksWOwave.csv')
WavFile = 'C:/Mel/CPODvsDPorCCA/CorrecteddB/738496579.150812000428.wav'
Clicks = pd.read_pickle('C:/Mel/CPODvsDPorCCA/CorrecteddB/Detected_Clips_110815_230428.pkl')
# clicks = detector.detect_click_clips_file(WavFile, blocksize=34560000)
# clicks = classifier.classify_matrix(clicks)

Start = 0  # 4731473
End = Start + 1024
click, Fs = soundfile.read(WavFile, start=int(Start), stop=int(End))

plt.plot(Clicks.wave[300])
plt.show()

nfft = 512
window = signal.get_window('boxcar', nfft)
freq, psd = signal.periodogram(x=click, window=window, nfft=nfft, fs=Fs, scaling='spectrum')
plt.plot(freq, psd)
plt.show()

# Normalize spectrum
psd = psd / np.max(psd)
cf = np.sum(freq * (psd ** 2)) / np.sum(psd ** 2)
pf = freq[psd.argmax()]

# Calculate RMSBW
# BW = (sqrt(sum((f-CF).^2.*PSD.^2 ) / sum(PSD.^2)))/1000;
rmsbw = (np.sqrt((np.sum(((freq - cf) ** 2) * (psd ** 2))) / np.sum(psd ** 2))) / 1000.0

# Calculate click duration based on Madsen & Walhberg 2007 - 80#
ener = np.cumsum(click ** 2)
istart = np.where(ener <= (ener[-1] * 0.1))[0]  # index of where the 1.5% is
iend = np.where(ener <= (ener[-1] * 0.9))[0]  # index of where the 98.5% is
if len(istart) > 0:
    istart = istart[-1]
else:
    istart = 0
if len(iend) > 0:
    iend = iend[-1]
else:
    iend = len(ener)
duration = ((iend - istart) * 1e6) / Fs  # duration in microseconds

# Parameters according to Mhl & Andersen, 1973
q = (cf / rmsbw) / 1000.0
ratio = pf / cf

# Calculate -3dB bandwith: Consecutive frequencies of the psd that have more than half of the maximum freq power
half = np.max(psd) / (10 ** (3 / 10.0))
max_freq_i = psd.argmax()
i_left, i_right = 0, 0
for i in np.arange(0, max_freq_i):
    if psd[max_freq_i - i] < half:
        break
    else:
        i_left = max_freq_i - i

for i in np.arange(0, psd.size - max_freq_i):
    if psd[max_freq_i + i] < half:
        break
    else:
        i_right = max_freq_i + i

bw = (freq[i_right] - freq[i_left]) / 1000.0

Clicks = pd.read_pickle('C:/Mel/CPODvsDPorCCA/TruncatedData/Detected_Clips_240815_170131.pkl')
