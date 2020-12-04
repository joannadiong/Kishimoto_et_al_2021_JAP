import os, shutil
import numpy as np
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import pandas as pd

import utilities as utils

LENOVO = '/home/joanna/Dropbox/Projects/activation'
EXTDRV = '/media/joanna/Elements/Projects/activation'
REPO = '/home/joanna/Dropbox/Projects/activation_repo'
os.chdir(REPO); print(os.getcwd())

df = pd.read_csv(os.path.join('.', 'data', 'proc', 'subjects_data.csv'))
subjects = utils.generate_subject_list()

slope_SO, intercept_SO = 0.027, -5.315
slope_MG, intercept_MG = 0.025, -5.512
slope_LG, intercept_LG = 0.028, -5.828

df['lnemgSO_pred'] = slope_SO * df['activations'] + intercept_SO
df['lnemgMG_pred'] = slope_MG * df['activations'] + intercept_MG
df['lnemgLG_pred'] = slope_LG * df['activations'] + intercept_LG

# perform linear regression on predicted vs observed EMG
# note, R2 values assume data are independent whereas they are correlated within participants
file = os.path.join('.', 'data', 'proc', 'results_pred-obs.txt')
open(file, 'w').close()
with open(file, 'a') as file:
    file.write('\n')
    for pred, obs, muscle in zip(['lnemgSO_pred', 'lnemgMG_pred', 'lnemgLG_pred'],
                                 ['lnemgSO', 'lnemgMG', 'lnemgLG'],
                                 ['SO', 'MG', 'LG']):
        emg_pred = df[pred]
        emg_obs = df[obs]
        md = smf.ols('emg_pred ~ emg_obs', df)
        md_fit = md.fit(reml=True)
        file.write('\nRegression for muscle: {}'.format(muscle))
        file.write('\n' + str(md_fit.summary()) + '\n')

# plot figure
fig = plt.figure(figsize=(12, 4.5))
ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)

lims = [-7.5, -2]

for subject in subjects:
    emg_obs = df['lnemgSO'][df['subject'] == subject]
    emg_pred = df['lnemgSO_pred'][df['subject'] == subject]
    ax1.scatter(emg_obs, emg_pred, s=4, c='0.5')
    # lims = [
    #     np.min([ax1.get_xlim(), ax1.get_ylim()]),  # min of both axes
    #     np.max([ax1.get_xlim(), ax1.get_ylim()]),  # max of both axes
    # ]

    emg_obs = df['lnemgMG'][df['subject'] == subject]
    emg_pred = df['lnemgMG_pred'][df['subject'] == subject]
    ax2.scatter(emg_obs, emg_pred, s=4, c='0.5')

    emg_obs = df['lnemgLG'][df['subject'] == subject]
    emg_pred = df['lnemgLG_pred'][df['subject'] == subject]
    ax3.scatter(emg_obs, emg_pred, s=4, c='0.5')

ax1.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax1.set_aspect('equal')
ax1.set_xlim(lims)
ax1.set_ylim(lims)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylabel('Predicted EMG SO [ln(mV)]')
ax1.set_xlabel('Observed EMG SO [ln(mV)]')
ax1.text(-8.5, -1, 'A', fontsize=14, weight='bold')

ax2.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax2.set_aspect('equal')
ax2.set_xlim(lims)
ax2.set_ylim(lims)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_ylabel('Predicted EMG MG [ln(mV)]')
ax2.set_xlabel('Observed EMG MG [ln(mV)]')
ax2.text(-8.5, -1, 'B', fontsize=14, weight='bold')

ax3.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax3.set_aspect('equal')
ax3.set_xlim(lims)
ax3.set_ylim(lims)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.set_ylabel('Predicted EMG LG [ln(mV)]')
ax3.set_xlabel('Observed EMG LG [ln(mV)]')
ax3.text(-8.5, -1, 'C', fontsize=14, weight='bold')

plt.tight_layout()
plt.savefig('fig-7' + '.png', dpi=300)
shutil.move('fig-7' + '.png', os.path.join('.', 'data', 'proc', 'fig-7' + '.png'))
plt.savefig('fig-7' + '.svg', dpi=300)
shutil.move('fig-7' + '.svg', os.path.join('.', 'data', 'proc', 'fig-7' + '.svg'))
