import os, shutil
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
import json
import pandas as pd
import utilities as utils
import trials_key
import spike2py as spk2
import warnings

import process_R1

warnings.filterwarnings('ignore')

LENOVO = '/home/joanna/Dropbox/Projects/activation'
EXTDRV = '/media/joanna/Elements/Projects/activation'
REPO = '/home/joanna/Dropbox/Projects/activation_repo'
os.chdir(REPO); print(os.getcwd()) # set as REPO when finalised

# To test in Console:
# subjects = ['sub01']

def import_process_signals(subjects):
    subjects_data = dict()
    for subject in subjects:
        _mkdir_proc(subject)
        sub_info, sub_data, sub_info_short = _import_signals(subject)
        sub_data = _calibrate_loadcell_signals(sub_info, sub_data)
        sub_data = _remove_loadcell_offset_start_each_trial(sub_info, sub_data)
        sub_data, max_vals_and_indexes, signals_above_threshold = _find_MVC_normalize_torque_signals(sub_info, sub_data)
        activations = _calculate_activations(sub_info, sub_data, max_vals_and_indexes)
        torques_emgs = _calculate_torque_EMG_at_activations(sub_info, sub_data)

        # R1 functions
        mvcs_by_time, levels_by_time = process_R1.find_trial_MVC_by_time(sub_info, sub_data)
        mmax_amp, idx1, idx2 = process_R1.find_mmax_amp(sub_info, sub_data)
        mmax_p1_rms = process_R1.find_mmax_rms(sub_info, sub_data, idx1, idx2)
        mvc_emg_rms = process_R1.find_mvc_emg_rms(sub_info, max_vals_and_indexes, signals_above_threshold)
        emgs_rect = process_R1.find_trial_emg_rms(sub_info, sub_data)
        emg_norm_mvc, emg_norm_mmax = process_R1.normalise_emg(sub_data, mvc_emg_rms, mmax_p1_rms, emgs_rect)

        subjects_data.update({subject: {'sub_info': sub_info_short,
                                        'mvc_torque': max_vals_and_indexes.mvc_torque[0],
                                        'activations': activations, # Access activations and torques with trial key
                                        'torques_emgs': torques_emgs,
                                        'mvcs_levels_by_time': [mvcs_by_time, levels_by_time],
                                        'mmax_amp': mmax_amp,
                                        'emg_norm_mvc': emg_norm_mvc,
                                        'emg_norm_mmax': emg_norm_mmax}})
    return subjects_data

def write_signals_to_json(subjects_data):
    path = os.path.join('.', 'data', 'proc')
    with open(os.path.join(path, 'subjects_data.json'), 'w') as file:
        json.dump(subjects_data, file)

def write_sub_info_to_csv(subjects, subjects_data):
    path = os.path.join('.', 'data', 'proc')
    df = pd.DataFrame()
    for subject in subjects:
        age, sex, height, weight, act_base, mvc_torque, mmax_amp = [[] for i in range(7)]
        age.append(subjects_data[subject]['sub_info']['age'])
        sex.append(subjects_data[subject]['sub_info']['sex'])
        height.append(subjects_data[subject]['sub_info']['height'])
        weight.append(subjects_data[subject]['sub_info']['weight'])
        activation = np.array(subjects_data[subject]['sub_info']['activations_baseline']).mean() # mean of 2 activation MVCs
        act_base.append(activation)
        mvc_torque.append(subjects_data[subject]['mvc_torque'])
        mmax_amp.append(subjects_data[subject]['mmax_amp'])
        df_ = pd.DataFrame({'subject': subject, 'age': age, 'sex': sex, 'height': height, 'weight': weight,
                            'act_base': act_base, 'mvc_torque': mvc_torque, 'mmax_amp': mmax_amp})
        df = df.append(df_, ignore_index=False)
    df.to_csv(os.path.join(path,'subjects_info.csv'))
    df.describe().to_csv(os.path.join(path, 'subjects_describe.csv'))

def write_signals_to_csv(subjects, subjects_data):
    path = os.path.join('.', 'data', 'proc')
    keys = ['01', '05', '10', '15', '25', '50', '75', '90', '95', '100']
    df = pd.DataFrame()
    for subject in subjects:
        sub = list((subject,) * len(keys))
        activations, torques, emgSO, emgMG, emgLG, emg_norm_mvc, emg_norm_mmax = [[] for i in range(7)]
        for key in keys:
            activations.append(subjects_data[subject]['activations'][key]['activation'])
            torques.append(subjects_data[subject]['torques_emgs'][key]['torque'])
            emgSO.append(subjects_data[subject]['torques_emgs'][key]['emgSO'])
            emgMG.append(subjects_data[subject]['torques_emgs'][key]['emgMG'])
            emgLG.append(subjects_data[subject]['torques_emgs'][key]['emgLG'])
            emg_norm_mvc.append(subjects_data[subject]['emg_norm_mvc'][key]['norm_mvc'])
            emg_norm_mmax.append(subjects_data[subject]['emg_norm_mmax'][key]['norm_mmax'])
        df_ = pd.DataFrame({'subject': sub, 'trials': keys, 'activations': activations, 'torques': torques,
                           'emgSO': emgSO, 'emgMG': emgMG, 'emgLG': emgLG,
                            'lnemgSO': np.log(emgSO), 'lnemgMG': np.log(emgMG), 'lnemgLG': np.log(emgLG),
                            'emg_norm_mvc': emg_norm_mvc, 'emg_norm_mmax': emg_norm_mmax})
        df = df.append(df_, ignore_index=False)
    df.to_csv(os.path.join(path,'subjects_data.csv'))

def plot_rest_twitches(subjects):
    path = os.path.join('.', 'data', 'proc')
    with open(os.path.join(path, 'subjects_data.json'), 'r') as file:
        data = json.load(file)
    # check REST twitches are representative (across trials) for each subject
    keys = ['01', '05', '10', '15', '25', '50', '75', '90', '95', '100']
    fig = plt.figure(figsize=(8, 6))
    plt.subplot(1, 1, 1)
    for subject in subjects:
        activations, rests = [], []
        for key in keys:
            activations.append(data[subject]['activations'][key]['activation'])
            rests.append(data[subject]['activations'][key]['rest_ptp'])
        plt.plot(activations, rests, label=subject)
    plt.legend()
    plt.ylabel('Torque (%MVC)')
    plt.xlabel('Activation (%)')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'rest_twitch.png'), dpi=300)
    plt.close()

def plot_signals(subjects):
    path = os.path.join('.', 'data', 'proc')
    with open(os.path.join(path, 'subjects_data.json'), 'r') as file:
        data = json.load(file)
    keys = ['01', '05', '10', '15', '25', '50', '75', '90', '95', '100']

    # plot EMG-activation
    fig = plt.figure(figsize=(9, 11))
    ax1 = fig.add_subplot(3, 2, 1)
    ax2 = fig.add_subplot(3, 2, 2)
    ax3 = fig.add_subplot(3, 2, 3)
    ax4 = fig.add_subplot(3, 2, 4)
    ax5 = fig.add_subplot(3, 2, 5)
    ax6 = fig.add_subplot(3, 2, 6)
    for subject in subjects:
        activation, emgSO, emgMG, emgLG = [[] for i in range(4)]
        for key in keys:
            activation.append(data[subject]['activations'][key]['activation'])
            emgSO.append(data[subject]['torques_emgs'][key]['emgSO'])
            emgMG.append(data[subject]['torques_emgs'][key]['emgMG'])
            emgLG.append(data[subject]['torques_emgs'][key]['emgLG'])
        # log the EMG signals
        emgSO_ln = np.log(emgSO)
        emgMG_ln = np.log(emgMG)
        emgLG_ln = np.log(emgLG)
        # plot
        ax1.plot(activation, emgSO)
        ax2.plot(activation, emgSO_ln)
        ax3.plot(activation, emgMG)
        ax4.plot(activation, emgMG_ln)
        ax5.plot(activation, emgLG)
        ax6.plot(activation, emgLG_ln)

    ax1.set_ylabel('EMG SO [mV]')
    ax2.set_ylabel('EMG SO [ln(mV)]')
    ax3.set_ylabel('EMG MG [mV]')
    ax4.set_ylabel('EMG MG [ln(mV)]')
    ax5.set_ylabel('EMG LG [mV]')
    ax5.set_xlabel('Activation [%]')
    ax6.set_ylabel('EMG LG [ln(mV)]')
    ax6.set_xlabel('Activation [%]')
    ax1.autoscale(enable=True, axis='x', tight=True)
    ax2.autoscale(enable=True, axis='x', tight=True)
    ax3.autoscale(enable=True, axis='x', tight=True)
    ax4.autoscale(enable=True, axis='x', tight=True)
    ax5.autoscale(enable=True, axis='x', tight=True)
    ax6.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'emg_activation.png'), dpi=300)
    plt.close()

    # plot activation-torque
    fig = plt.figure(figsize=(11, 7))
    for subject in subjects:
        activation, torque = [[] for i in range(2)]
        for key in keys:
            activation.append(data[subject]['activations'][key]['activation'])
            torque.append(data[subject]['torques_emgs'][key]['torque'])
        plt.plot(torque, activation, '0.5')
    plt.ylabel('Activation [%]')
    plt.xlabel('Torque [%MVC]')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'activation_torque.png'), dpi=300)
    plt.close()

    # plot emg-torque
    fig = plt.figure(figsize=(7, 11))
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)
    for subject in subjects:
        torque, emgSO, emgMG, emgLG = [[] for i in range(4)]
        for key in keys:
            torque.append(data[subject]['torques_emgs'][key]['torque'])
            emgSO.append(data[subject]['torques_emgs'][key]['emgSO'])
            emgMG.append(data[subject]['torques_emgs'][key]['emgMG'])
            emgLG.append(data[subject]['torques_emgs'][key]['emgLG'])
        # plot
        ax1.plot(torque, emgSO, '0.5')
        ax2.plot(torque, emgMG, '0.5')
        ax3.plot(torque, emgLG, '0.5')

    ax1.set_ylabel('EMG SO [mV]')
    ax2.set_ylabel('EMG MG [mV]')
    ax3.set_ylabel('EMG LG [mV]')
    ax3.set_xlabel('Torque [%MVC]')
    ax1.autoscale(enable=True, axis='x', tight=True)
    ax2.autoscale(enable=True, axis='x', tight=True)
    ax3.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'emg_torque.png'), dpi=300)
    plt.close()

def plot_mvcs_over_time(subjects):
    path = os.path.join('.', 'data', 'proc')
    with open(os.path.join(path, 'subjects_data.json'), 'r') as file:
        data = json.load(file)

    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(1, 1, 1)
    for subject in subjects:
        # print (subject)
        mvc_start = data[subject]['mvc_torque']
        mvcs_by_time = data[subject]['mvcs_levels_by_time'][0]
        order = [1, 2, 3, 4, 5, 6, 7]
        if len(mvcs_by_time) < len(order):
            order = order[: len(mvcs_by_time)]
        # plt.plot([1] + order, [mvc_start] + mvcs_by_time, 'k')
        # plt.plot(1, mvc_start, 'ro')
        plt.plot(order, mvcs_by_time, 'k')
        plt.plot(order, mvcs_by_time, 'ko')
    ax.set_ylim(0, 160)
    plt.ylabel('Torque [%MVC]')
    plt.xlabel('Test order [a.u.]')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'mvcs_by_time.png'), dpi=300)
    plt.savefig(os.path.join(path, 'mvcs_by_time.svg'), dpi=300)
    plt.close()

def plot_normalised_emg(subjects):
    path = os.path.join('.', 'data', 'proc')
    with open(os.path.join(path, 'subjects_data.json'), 'r') as file:
        data = json.load(file)
    keys = ['01', '05', '10', '15', '25', '50', '75', '90', '95', '100']

    # plot figure
    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)

    slopes = []
    for subject in subjects:
        # extract activation, EMG data from subject dict
        activations, emg_norm_mvc, emg_norm_mmax = [], [], []
        for key in keys:
            activations.append(data[subject]['activations'][key]['activation'])
            emg_norm_mvc.append(data[subject]['emg_norm_mvc'][key]['norm_mvc'])
            emg_norm_mmax.append(data[subject]['emg_norm_mmax'][key]['norm_mmax'])

        # get slopes
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(emg_norm_mmax, emg_norm_mvc)
        slopes.append(slope)

        # plot
        ax1.plot(activations, emg_norm_mvc, 'k')
        ax2.plot(activations, emg_norm_mmax, 'k')
        ax3.plot(emg_norm_mmax, emg_norm_mvc, 'k')

    # analyse all slopes
    mean, sd = np.mean(slopes), np.std(slopes)
    se = sd / np.sqrt(len(slopes))
    ll, ul = stats.norm.interval(0.95, loc=mean, scale=se)
    # ll, ul = mean - 1.96 * se, mean + 1.96 * se # checked -- correct

    file = os.path.join(path, 'results.txt')
    open(file, 'w').close()
    with open(file, 'a') as file:
        file.write('\nMean {:.2f}, SD {:.2f}, 95% CI {:.2f} to {:.2f}'.format(mean, sd, ll, ul))

    # ax1.set_ylim(0, 140)
    ax1.set_ylabel('EMG (%MVC)')
    ax1.set_xlabel('Activation (%)')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.autoscale(enable=True, axis='x', tight=True)

    ax2.set_ylabel('EMG (%Mmax)')
    ax2.set_xlabel('Activation (%)')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.autoscale(enable=True, axis='x', tight=True)

    ax3.set_ylabel('EMG (%MVC)')
    ax3.set_xlabel('EMG (%Mmax)')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.autoscale(enable=True, axis='x', tight=True)

    plt.tight_layout()
    plt.savefig(os.path.join(path, 'emg_normalised.png'), dpi=300)
    plt.savefig(os.path.join(path, 'emg_normalised.svg'), dpi=300)
    plt.close()

def copy_figs(subjects):
    if not os.path.exists(os.path.join('.', 'data', 'proc', 'figs')):
        os.chdir(os.path.join('.', 'data', 'proc'))
        os.mkdir('figs')
        os.chdir(os.path.join('..', '..'))  # return to base directory: data

    for subject in subjects:
        shutil.copy(os.path.join('.', 'data', 'proc', subject, 'mmax-p1.png'),
                    os.path.join('.', 'data', 'proc', 'figs', 'mmax-p1-' + subject + '.png'))
        shutil.copy(os.path.join('.', 'data', 'proc', subject, 'mvc_torq_emg.png'),
                    os.path.join('.', 'data', 'proc', 'figs', 'mvc_torq_emg-' + subject + '.png'))
        shutil.copy(os.path.join('.', 'data', 'proc', subject, 'emg_rect.png'),
                    os.path.join('.', 'data', 'proc', 'figs', 'emg_rect-' + subject + '.png'))

def _mkdir_proc(subject):

    if not os.path.exists(os.path.join('.', 'data', 'proc', subject)):
        os.chdir(os.path.join('.', 'data', 'proc'))
        os.mkdir(subject)
        os.chdir(os.path.join('..', '..'))  # return to base directory: data

def _import_signals(subject):

    path = os.path.join('.', 'data', 'raw', subject)
    sub_info = utils.read_subject_log(path, subject)
    sub_key = trials_key.gen(sub_info.sub)
    signals = utils.generate_spike2py_signalInfos(sub_info)
    sub_info_short = {'age': sub_info.age, 'sex': sub_info.sex, 'height':sub_info.height, 'weight': sub_info.weight,
                      'activations_baseline': sub_info.activations_baseline}

    trial_data = dict()
    for trial, trialname in sub_key.items():
        filename = trialname + '.mat'
        trial_info = spk2.TrialInfo(cond=trial, path=path, filename=filename, signals=signals)
        trial_data[trial] = spk2.Trial(trial_info)

    return sub_info, trial_data, sub_info_short

def _calibrate_loadcell_signals(sub_info, sub_data):

    loadcell_offset_value = np.mean(sub_data['baseline'].sig['torque'].raw) * sub_info.scale_MVC_loadcell
    for key in sub_data.keys():
        sub_data[key].sig['torque'].calibrate(slope=sub_info.scale_MVC_loadcell, offset=loadcell_offset_value)

    return sub_data

def _remove_loadcell_offset_start_each_trial(sub_info, sub_data):

    for key in sub_data.keys():
        sub_data[key].sig['torque'].remove_offset(type_='start', val=int(sub_info.freq))

    return sub_data

def _find_MVC_normalize_torque_signals(sub_info, sub_data):

    torque = sub_data['mvc_vol'].sig['torque'].proc
    index_above_threshold = list(torque > 30) # set torque threshold at 30 Nm
    count = 0
    indexes = []
    # Extract indexes of torque data during the 5 MVC attempts (last 5 MVCs)
    for i in range(len(index_above_threshold)-1, 1, -1):
        if index_above_threshold[i]:
            indexes.append(i)
            if not index_above_threshold[i-1]:
                count += 1
                if count == 5:
                    break
    mvc_torque = max(torque[indexes])
    mvc_torque_index = np.argmax(torque[indexes])
    mvc_torque = (mvc_torque, mvc_torque_index) # only torques above threshold are indexed, not time series torques
    torques_above_threshold = torque[indexes]

    fig = plt.figure(figsize=(11, 7))
    plt.subplot(1, 1, 1)
    plt.grid()
    plt.plot(torque[indexes], 'k')
    plt.plot(mvc_torque[1], mvc_torque[0] + 2, 'ro')
    plt.ylabel('Torque (Nm)')
    plt.tight_layout()
    plt.savefig('mvc_vol.png', dpi=300)
    shutil.move('mvc_vol.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'mvc_vol.png'))
    plt.close()

    def _find_max_EMG(emg_signal, indexes):
        mvc_emg = max(emg_signal[indexes])
        mvc_emg_index = np.argmax(emg_signal[indexes])
        return mvc_emg, mvc_emg_index

    SO = sub_data['mvc_vol'].sig['emgSO'].envel
    MG = sub_data['mvc_vol'].sig['emgMG'].envel
    LG = sub_data['mvc_vol'].sig['emgLG'].envel
    mvc_SO = _find_max_EMG(SO, indexes)
    mvc_MG = _find_max_EMG(MG, indexes)
    mvc_LG = _find_max_EMG(LG, indexes)

    emgSO_above_threshold = sub_data['mvc_vol'].sig['emgSO'].rect[indexes]
    emgMG_above_threshold = sub_data['mvc_vol'].sig['emgMG'].rect[indexes]
    emgLG_above_threshold = sub_data['mvc_vol'].sig['emgLG'].rect[indexes]

    for i, val in enumerate(torque):
        if val > 40:            # Could break if torque goes above 40 during DF MVC
            last_index_for_TA = i
            break

    TA = sub_data['mvc_vol'].sig['emgTA'].envel
    TA_thresholded = list(TA[0:last_index_for_TA] > 0.05)
    indexes_for_TA_MVC = []
    current_run_True = []
    for i, val in enumerate(TA_thresholded):
        if val:
            current_run_True.append(i)
        else:
            if len(current_run_True) >= 500:
                indexes_for_TA_MVC.extend(current_run_True)
            current_run_True = []
    mvc_TA = _find_max_EMG(TA, indexes)

    Maximum_values = namedtuple('Maximum_values', 'mvc_torque mvc_SO mvc_MG mvc_LG mvc_TA')
    max_vals_and_indexes = Maximum_values(mvc_torque=mvc_torque,
                                        mvc_SO=mvc_SO,
                                        mvc_MG=mvc_MG,
                                        mvc_LG=mvc_LG,
                                        mvc_TA=mvc_TA)
    Above_threshold = namedtuple('Above_threshold', 'torques emgSO emgMG emgLG')
    signals_above_threshold = Above_threshold(torques=torques_above_threshold,
                                              emgSO=emgSO_above_threshold,
                                              emgMG=emgMG_above_threshold,
                                              emgLG=emgLG_above_threshold)

    for key in sub_data.keys():
        sub_data[key].sig['torque'].normalize(type_='value', value=max_vals_and_indexes.mvc_torque[0], signal_version='proc')
        # sub_data[key].sig['emgSO'].normalize(type_='value', value=max_vals_and_indexes.mvc_SO[0], signal_version='rect')
        # sub_data[key].sig['emgMG'].normalize(type_='value', value=max_vals_and_indexes.mvc_MG[0], signal_version='rect')
        # sub_data[key].sig['emgLG'].normalize(type_='value', value=max_vals_and_indexes.mvc_LG[0], signal_version='rect')
        # sub_data[key].sig['emgTA'].normalize(type_='value', value=max_vals_and_indexes.mvc_TA[0], signal_version='rect')

    return sub_data, max_vals_and_indexes, signals_above_threshold

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

def _calculate_activations(sub_info, sub_data, max_vals_and_indexes):

    def _calculate_peak_to_peak_amplitude(sub, key, freq, index, type):
        high_force_trials = ['90', '95', '100']
        if key in high_force_trials and type=='sit':
            nsamples = int(freq * 0.150) # Finds twitch peak at high force within 150 ms (could customise to 100 ms)
        else:
            nsamples = int(freq * 0.150) # Finds twitch peak at low-mod force within 150 ms
        index1, index2 = index, index + nsamples
        # find min and max torque, accounting for 15 ms electromechanical delay
        proc = sub_data[key].sig['torque'].proc # 50 Hz lowpass filtered
        time = sub_data[key].sig['torque'].times
        sig = proc[index1: index2]
        delay = int(0.010 * freq) # default: 0.015
        sig_after_delay = sig[delay:]
        time_after_delay = time[index1: index2][delay:]
        # find index and value of max torque
        index_max_nsamples, sig_max = np.where(sig == sig_after_delay.max())[0], sig_after_delay.max()
        # find index and value of min force in signal preceding the max force
        # sig_before_max = sig[:int(index_max_nsamples)] # signal is between stimulus and max torque
        # time_before_max = time[:int(index_max_nsamples)]
        sig_before_max = sig[delay: int(index_max_nsamples)] # signal is between EMD and max torque
        time_before_max = time[delay: int(index_max_nsamples)]
        if int(index_max_nsamples) == delay:
            index_min_nsamples, sig_min = index_max_nsamples, sig_max
        else:
            index_min_nsamples, sig_min = np.where(sig == sig_before_max.min())[0], sig_before_max.min()
        # calculate twitch amplitude
        signal_ptp = sig_max - sig_min
        if signal_ptp < 0:
            signal_ptp = 0
        index_min = index + int(index_min_nsamples)
        index_max = index + int(index_max_nsamples)
        # plot and check indexing of torque
        raw = sub_data[key].sig['torque'].raw  # torque in V
        loadcell_offset_value = np.mean(sub_data['baseline'].sig['torque'].raw) * sub_info.scale_MVC_loadcell
        raw = raw * sub_info.scale_MVC_loadcell - loadcell_offset_value  # torque in Nm
        raw = raw / max_vals_and_indexes.mvc_torque[0] * 100  # torque normalised to MVC
        raw = raw - np.mean(raw[:2000])  # unfiltered
        plt.figure()
        plt.plot(time[index1: index2], raw[index1: index2], label='raw')
        # plt.plot(time[index1: index2], sig, label='filtered')
        plt.plot(time_after_delay, sig_after_delay, label='filt, after EMD')
        plt.plot(delay / freq + index1 / freq, sig[delay], 'ko', label='EMD={}ms'.format(int(delay / freq * 1000)))
        plt.plot(index_max_nsamples / freq + index1 / freq, sig_max, 'ro', label='max')
        plt.plot(time_before_max + index1 / freq, sig_before_max, label='filt, before max')
        plt.plot(index_min_nsamples / freq + index1 / freq, sig_min, 'go', label='min')
        plt.legend()
        plt.ylabel('Torque (%MVC)')
        plt.xlabel('Time within window (s)')
        if signal_ptp == 0:
            text = '{}, {}%; Max - min: {:.4} - {:.4} = {}'.format(sub_info.sub, key, sig_max, sig_min, signal_ptp)
        else:
            text = '{}, {}%; Max - min: {:.4} - {:.4} = {:.4}'.format(sub_info.sub, key, sig_max, sig_min, signal_ptp)
        plt.annotate(text, xy=(0.01, 1.01), xycoords='axes fraction', fontsize=8)
        plt.tight_layout()
        plt.savefig('sit_' + key + '.png', dpi=300)
        shutil.move('sit_' + key + '.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'sit_' + key + '.png'))
        plt.close()
        return signal_ptp, index_min, index_max

    def _plot_signals(sub_data, sub_info, key, rest_idx1, rest_idx2, sit_idx1, sit_idx2):
        torque = sub_data[key].sig['torque'].proc
        emgSO = sub_data[key].sig['emgSO'].proc
        emgMG = sub_data[key].sig['emgMG'].proc
        emgLG = sub_data[key].sig['emgLG'].proc
        fig = plt.figure(figsize=(11, 7))
        # torque
        plt.subplot(2, 1, 1)
        plt.grid()
        plt.plot(torque, 'k')
        rest_idxs, sit_idxs = np.arange(rest_idx1, rest_idx2), np.arange(sit_idx1, sit_idx2)
        rest_torque, sit_torque = torque[rest_idx1: rest_idx2], torque[sit_idx1: sit_idx2]
        plt.plot(rest_idxs, rest_torque, 'b', label='REST')
        plt.plot(sit_idxs, sit_torque, 'r', label='SIT')
        plt.legend()
        plt.autoscale(enable=True, axis='x', tight=True)
        text = sub_info.sub + ': ' + key + '%'
        plt.annotate(text, xy=(0, 1), xycoords='axes fraction', fontsize=8)
        plt.ylabel('Torque (%MVC)')
        # EMG
        plt.subplot(2, 1, 2)
        plt.grid()
        plt.plot(emgSO, 'k', label='SO')
        plt.plot(emgMG, 'g', label='MG')
        plt.plot(emgLG, 'b', label='LG')
        plt.legend()
        plt.ylabel('EMG (mV)')
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.tight_layout()
        plt.savefig(key + '.png', dpi=300)
        shutil.move(key + '.png', os.path.join('.', 'data', 'proc', sub_info.sub, key + '.png'))
        plt.close()

    print('\n' + sub_info.sub)
    activations = dict()
    for key in list(sub_data.keys())[4:]:
        torque = sub_data[key].sig['torque'].proc
        index_rest, index_sit = _determine_sit_rest_indexes(sub_info, sub_data, key) # , max_vals_and_indexes

        rest_ptp, rest_idx1, rest_idx2 = _calculate_peak_to_peak_amplitude(sub_info.sub, key, sub_info.freq, index_rest, type='rest')
        sit_ptp, sit_idx1, sit_idx2 = _calculate_peak_to_peak_amplitude(sub_info.sub, key, sub_info.freq, index_sit, type='sit')
        print(f'key: {key}%, SIT is after REST (1-75% MVC): {index_sit-index_rest > 0}')

        activation = (1 - (sit_ptp / rest_ptp)) * 100
        activation_ = {key: {'rest_ptp': rest_ptp, 'sit_ptp': sit_ptp, 'activation': activation}}
        activations.update(activation_)

        _plot_signals(sub_data, sub_info, key, rest_idx1, rest_idx2, sit_idx1, sit_idx2)

    return activations

def _calculate_torque_EMG_at_activations(sub_info, sub_data):

    nsamples_before_trig = int(sub_info.freq * 0.05) # mean EMG over 50 ms window
    torques_emgs = dict()

    for key in list(sub_data.keys())[4:]:
        index_rest, index_sit = _determine_sit_rest_indexes(sub_info, sub_data, key)
        # shift indexed EMG region away from filter artefact close to stimulus artefact
        filter_artefact_length = int(sub_info.freq * 0.05)
        index_start, index_stop = index_sit - (filter_artefact_length + nsamples_before_trig), index_sit - filter_artefact_length

        torque = np.mean(sub_data[key].sig['torque'].proc[index_start: index_stop])
        emgSO = np.mean(sub_data[key].sig['emgSO'].rect[index_start: index_stop])
        emgMG = np.mean(sub_data[key].sig['emgMG'].rect[index_start: index_stop])
        emgLG = np.mean(sub_data[key].sig['emgLG'].rect[index_start: index_stop])
        emgTA = np.mean(sub_data[key].sig['emgTA'].rect[index_start: index_stop])

        torques_emgs_ = {key: {'torque': torque, 'emgSO': emgSO, 'emgMG': emgMG, 'emgLG': emgLG, 'emgTA': emgTA}}
        torques_emgs.update(torques_emgs_)

        # plot and check indexing of EMG
        i = index_start - int(sub_info.freq * 0.01)
        j = index_sit + int(sub_info.freq * 0.01)
        plt.figure()
        # EMG SO
        plt.subplot(3, 1, 1)
        emg = sub_data[key].sig['emgSO'].rect[i: j]
        time = sub_data[key].sig['emgSO'].times[i: j]
        plt.plot(time, emg, 'k')
        plt.plot(index_start / sub_info.freq, 0, 'g|', linewidth=5, label='start')
        plt.plot(index_stop / sub_info.freq, 0, 'r|', linewidth=5, label='stop')
        plt.ylabel('EMG SO (mV)')
        plt.legend()
        # EMG MG
        plt.subplot(3, 1, 2)
        emg = sub_data[key].sig['emgMG'].rect[i: j]
        time = sub_data[key].sig['emgMG'].times[i: j]
        plt.plot(time, emg, 'k')
        plt.plot(index_start / sub_info.freq, 0, 'g|', linewidth=5, label='start')
        plt.plot(index_stop / sub_info.freq, 0, 'r|', linewidth=5, label='stop')
        plt.ylabel('EMG MG (mV)')
        # EMG LG
        plt.subplot(3, 1, 3)
        emg = sub_data[key].sig['emgLG'].rect[i: j]
        time = sub_data[key].sig['emgLG'].times[i: j]
        plt.plot(time, emg, 'k')
        plt.plot(index_start / sub_info.freq, 0, 'g|', linewidth=5, label='start')
        plt.plot(index_stop / sub_info.freq, 0, 'r|', linewidth=5, label='stop')
        plt.ylabel('EMG LG (mV)')
        plt.xlabel('Time within window (s)')
        plt.tight_layout()
        plt.savefig('emg_' + key + '.png', dpi=300)
        shutil.move('emg_' + key + '.png', os.path.join('.', 'data', 'proc', sub_info.sub, 'emg_' + key + '.png'))
        plt.close()

    return torques_emgs

