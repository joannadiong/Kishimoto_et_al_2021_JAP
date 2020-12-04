import os, shutil
import numpy as np
import scipy.signal
from scipy import interpolate
import matplotlib.pyplot as plt

import trials_key

def find_trial_MVC_by_time(sub_info, sub_data):
    # get activation levels and trials key-value pairs
    sub_key = trials_key.gen(sub_info.sub)

    # make dict of trials sorted by activation level
    levels = ['01', '05', '10', '15', '25', '50', '75']
    bysort_level = {}
    for level in levels:
        trial = sub_key[level]
        bysort_level.update({level: trial})

    # make dict of trials sorted by test order over time
    bysort_time = {k: v for k, v in sorted(bysort_level.items(), key=lambda item: item[1])}
    levels_by_time = list(bysort_time.keys())

    # get MVC torques sorted by tested order over time
    mvcs_by_time = []
    levels_by_time_ = [] # get activation level only if MVC exists
    for level in levels_by_time:
        mvc_torque = _find_trial_MVC_normalize_torque_signals(sub_info, sub_data, level)
        if mvc_torque:
            mvcs_by_time.append(mvc_torque[0])
            levels_by_time_.append(level)

    return mvcs_by_time, levels_by_time_


def find_mmax_amp(sub_info, sub_data):
    # index the last maximal stimulation
    idx = int(sub_data['max_curr'].sig['trig'].times[-1] * sub_info.freq)
    # calculate the peak to peak within a 50 ms window, 5 ms after the stimulus
    ptp_start, ptp_stop = 0.005, 0.055 # in sec
    idx1 = int(idx + ptp_start * sub_info.freq) 
    idx2 = int(idx1 + ptp_stop * sub_info.freq)

    emg = sub_data['max_curr'].sig['emgSO'].raw
    time = sub_data['max_curr'].sig['emgSO'].times

    mmax_amp = np.ptp(emg[idx1: idx2])
    buffer = int(0.040 * sub_info.freq) # in ms

    # plot window and check that it is correct
    fig = plt.figure(figsize=(11, 7))
    plt.subplot(1, 3, 1)
    plt.grid()
    plt.plot((time[idx1: idx2] - time[idx1]) * 1000, emg[idx1: idx2], 'r')
    plt.xlabel('Time (ms)')
    plt.ylabel('EMG SO (mV)')
    plt.subplot(1, 3, (2, 3))
    # set stim instance to zero, plot time in ms
    plt.plot(0, emg[idx] + max(emg), 'r|', markersize=6) # idx / sub_info.freq
    plt.plot((time[idx - buffer: idx2 + buffer] - time[idx]) * 1000,
             emg[idx - buffer: idx2 + buffer], 'k')
    plt.plot((time[idx1: idx2] - time[idx1] + ptp_start) * 1000, emg[idx1: idx2], 'r')
    plt.xlabel('Time (ms)')
    plt.ylabel('EMG SO (mV)')
    plt.tight_layout()
    plt.savefig('mmax.png', dpi=300)
    shutil.move('mmax.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'mmax.png'))
    plt.close()

    return mmax_amp, idx1, idx2


def find_mmax_rms(sub_info, sub_data, idx1, idx2):
    # identify samples over which soleus M wave occurs; only soleus M waves were checked
    emg = sub_data['max_curr'].sig['emgSO'].raw
    time = sub_data['max_curr'].sig['emgSO'].times
    mmax = emg[idx1: idx2]

    # interpolate over the M wave
    xaxis = list(range(0, len(mmax)))
    f = interpolate.interp1d(xaxis, mmax)
    xaxis_new = np.arange(0, len(mmax) - 1, 0.1)
    mmax_new = f(xaxis_new)

    # identify the sample indexes where the first phase of the M wave crosses 0 volts
    # similarly to Thomas C (1997) Fatigue in human thenar muscles paralysed by spinal cord injury
    min_val = abs(min(mmax_new))
    max_val = max(mmax_new)
    height = np.mean([min_val, max_val]) * .7
    indexes, _ = scipy.signal.find_peaks(abs(mmax_new), height=height, distance=5)
    plt.plot(mmax_new,'.-')
    plt.plot(indexes, mmax_new[indexes], 'ro', label='min, max')

    peak_index = indexes[0]
    if mmax_new[peak_index] < 0:
        mmax_new *= -1
    for i in range(peak_index, 0, -1):
        if mmax_new[i] > 0 and mmax_new[i - 1] < 0:
            idx_start_p1_mmax = i-1
            break
        else:
            idx_start_p1_mmax = 0
    for i in range(peak_index, len(mmax_new)):
        if mmax_new[i] > 0 and mmax_new[i + 1] < 0:
            idx_stop_p1_mmax = i + 1
            break

    # print(idx_start_p1_mmax, idx_stop_p1_mmax)
    plt.plot([idx_start_p1_mmax, idx_stop_p1_mmax],
             [mmax_new[idx_start_p1_mmax], mmax_new[idx_stop_p1_mmax]],
             'bo', label='cross 0V')
    plt.legend()
    plt.xlabel('Samples')
    plt.ylabel('EMG SO (mV)')
    plt.tight_layout()
    plt.savefig('mmax-p1.png', dpi=300)
    shutil.move('mmax-p1.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'mmax-p1.png'))
    plt.close()

    # calculate the root-mean-square of the first phase of the M wave
    # from sklearn.metrics import mean_squared_error
    # i = np.zeros(len(mmax_new[idx_start_p1_mmax: idx_stop_p1_mmax]))
    # np.sqrt(mean_squared_error(i, mmax_new[idx_start_p1_mmax: idx_stop_p1_mmax])) # gets same answer
    mmax_p1_rms = np.sqrt(np.sum(mmax_new[idx_start_p1_mmax: idx_stop_p1_mmax] ** 2) / len(mmax_new[idx_start_p1_mmax: idx_stop_p1_mmax]))

    return mmax_p1_rms


def find_mvc_emg_rms(sub_info, max_vals_and_indexes, signals_above_threshold):
    # use index of MVC torque to find MVC EMG
    mvc_torque_idx = max_vals_and_indexes.mvc_torque[1]
    mvc_torque = max_vals_and_indexes.mvc_torque[0]
    torque = signals_above_threshold.torques
    emgSO = signals_above_threshold.emgSO # rectified EMG, not enveloped

    plt.subplot(3,1,1)
    plt.plot(torque, 'k')
    plt.plot(mvc_torque_idx, mvc_torque, 'ro')
    plt.ylabel('Torque (Nm)')

    plt.subplot(3,1,2)
    plt.plot(emgSO, 'k')
    plt.plot(mvc_torque_idx, emgSO[mvc_torque_idx], 'ro')
    plt.ylabel('EMG SO (mV)')

    # get root mean square soleus EMG over 50 ms window over the MVC index
    half_win = int(sub_info.freq * 0.05 / 2)
    mvc_indexes = list(range(mvc_torque_idx - half_win, mvc_torque_idx + half_win))
    mvc_emg = emgSO[mvc_torque_idx - half_win: mvc_torque_idx + half_win]
    mvc_emg_rms = np.sqrt(np.sum(mvc_emg ** 2) / len(mvc_emg))
    # from sklearn.metrics import mean_squared_error
    # i = np.zeros(len(mvc_emg))
    # np.sqrt(mean_squared_error(i, mvc_emg)) # gets same answer

    plt.subplot(3,1,3)
    plt.plot(list(range(mvc_torque_idx - int(sub_info.freq / 6), mvc_torque_idx + int(sub_info.freq / 6))),
             emgSO[mvc_torque_idx - int(sub_info.freq / 6): mvc_torque_idx + int(sub_info.freq / 6)], 'k')
    for i, j in zip(mvc_indexes, mvc_emg):
        plt.plot(i, j, 'go', markersize=2)
    plt.xlabel('Samples')
    plt.ylabel('EMG SO (mV)')
    plt.tight_layout()
    plt.savefig('mvc_torq_emg.png', dpi=300)
    shutil.move('mvc_torq_emg.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'mvc_torq_emg.png'))
    plt.close()

    return mvc_emg_rms


def find_trial_emg_rms(sub_info, sub_data):
    nsamples_before_trig = int(sub_info.freq * 0.05)  # get EMG over 50 ms window
    emgs_rect = dict()

    i = 1
    j = len(list(sub_data.keys())[4:])

    def _determine_sit_rest_indexes(sub_info, sub_data, key):

        nsamples_before_trig = int(sub_info.freq * 0.5)
        idx1 = int(sub_data[key].sig['trig'].times[0] * sub_info.freq)
        idx2 = int(sub_data[key].sig['trig'].times[1] * sub_info.freq)
        if np.mean(sub_data[key].sig['torque'].proc[idx1 - nsamples_before_trig: idx1]) > \
                np.mean(sub_data[key].sig['torque'].proc[idx2 - nsamples_before_trig: idx2]):
            index_sit = idx1
            index_rest = idx2
        else:
            index_sit = idx2
            index_rest = idx1

        return index_rest, index_sit

    for key in list(sub_data.keys())[4:]:
        index_rest, index_sit = _determine_sit_rest_indexes(sub_info, sub_data, key)
        # shift indexed EMG region away from filter artefact close to stimulus artefact
        filter_artefact_length = int(sub_info.freq * 0.05)
        index_start, index_stop = index_sit - (filter_artefact_length + nsamples_before_trig), index_sit - filter_artefact_length

        emgSO = sub_data[key].sig['emgSO'].rect[index_start: index_stop]
        emgMG = sub_data[key].sig['emgMG'].rect[index_start: index_stop]
        emgLG = sub_data[key].sig['emgLG'].rect[index_start: index_stop]
        emgs_rect_ = {key: {'emgSO': emgSO, 'emgMG': emgMG, 'emgLG': emgLG}}
        emgs_rect.update(emgs_rect_)

        plt.subplot(j, 1, i)
        plt.plot(emgSO, 'k', label='SO')
        plt.plot(emgMG, 'r', label='MG')
        plt.plot(emgLG, 'b', label='LG')
        plt.ylim(0, 0.5)
        plt.yticks(ticks=[], labels=[])
        if i == 2:
            plt.legend()
        if i == 6:
            plt.ylabel('EMG (ylim 0-0.2 mV)')
        i += 1
    plt.xlabel('Samples')
    plt.tight_layout()
    plt.savefig('emg_rect.png', dpi=300)
    shutil.move('emg_rect.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'emg_rect.png'))
    plt.close()

    return emgs_rect


def normalise_emg(sub_data, mvc_emg_rms, mmax_p1_rms, emgs_rect):
    emg_norm_mvc, emg_norm_mmax = dict(), dict()
    for key in list(sub_data.keys())[4:]:
        emg = emgs_rect[key]['emgSO']
        emg_rms = np.sqrt(np.sum(emg ** 2) / len(emg))
        # from sklearn.metrics import mean_squared_error
        # i = np.zeros(len(emg))
        # np.sqrt(mean_squared_error(i, emg)) # gets same answer

        emg_mvc = emg_rms / mvc_emg_rms * 100
        emg_norm_mvc_ = {key: {'norm_mvc': emg_mvc}}
        emg_norm_mvc.update(emg_norm_mvc_)

        emg_mmax = emg_rms / mmax_p1_rms * 100
        emg_norm_mmax_ = {key: {'norm_mmax': emg_mmax}}
        emg_norm_mmax.update(emg_norm_mmax_)

    return emg_norm_mvc, emg_norm_mmax


def _find_trial_MVC_normalize_torque_signals(sub_info, sub_data, level):

    # get calibrated max torque during the MVC for activation levels 1-75% trials
    torque = sub_data[level].sig['torque'].proc
    index_above_threshold = list(torque > 30) # set torque threshold at 30 Nm
    count = 0
    indexes = []

    # Extract indexes of torque data during the MVC preconditioning
    for i in range(0, len(index_above_threshold), 1):
        if index_above_threshold[i]:
            indexes.append(i)
            if not index_above_threshold[i-1]:
                count += 1
                if count == 5:
                    break
    if not indexes: # if MVC was not recorded
        pass
    else:
        mvc_torque = max(torque[indexes])
        mvc_torque_index = np.argmax(torque[indexes])
        mvc_torque = (mvc_torque, mvc_torque_index) # only torques above threshold are indexed

        fig = plt.figure(figsize=(11, 7))
        plt.subplot(1, 3, 1)
        plt.grid()
        plt.plot(torque[indexes], 'k')
        plt.plot(mvc_torque[1], mvc_torque[0] + 2, 'ro')
        plt.ylabel('Torque (Nm)')
        plt.subplot(1, 3, (2, 3))
        plt.plot(torque, 'k')
        plt.ylabel('Torque (Nm)')
        plt.tight_layout()
        plt.savefig('mvc_' + level + '.png', dpi=300)
        shutil.move('mvc_' + level + '.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'mvc_' + level + '.png'))
        plt.close()

        return mvc_torque

