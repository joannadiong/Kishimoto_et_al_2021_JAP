import os, shutil
import matplotlib.pyplot as plt

import process

subject = 'sub01'

sub_info, sub_data, sub_info_short = process._import_signals(subject)
sub_data = process._calibrate_loadcell_signals(sub_info, sub_data)
sub_data = process._remove_loadcell_offset_start_each_trial(sub_info, sub_data)
sub_data, max_vals_and_indexes, signals_above_threshold = process._find_MVC_normalize_torque_signals(sub_info, sub_data)

keys = ['05', '100']

fig = plt.figure(figsize=(11, 7))
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)

# torque: 5% activation
torque = sub_data[keys[0]].sig['torque'].proc
ymax = max(torque) + 5
times = sub_data[keys[0]].sig['torque'].times
stim_times = sub_data[keys[0]].sig['trig'].times
ax1.plot(times, torque, 'k')
ax1.plot(stim_times[0], ymax, 'k|', linewidth=5)
ax1.plot(stim_times[1], ymax, 'k|', linewidth=5)
ax1.autoscale(enable=True, axis='x', tight=True)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylabel('Torque [%MVC]')
ax1.set_xlabel('Time [s]')

# torque: 100% activation
torque = sub_data[keys[1]].sig['torque'].proc
ymax = max(torque) + 5
times = sub_data[keys[1]].sig['torque'].times
stim_times = sub_data[keys[1]].sig['trig'].times
ax2.plot(times, torque, 'k')
ax2.plot(stim_times[0], ymax, 'k|', linewidth=5)
ax2.plot(stim_times[1], ymax, 'k|', linewidth=5)
ax2.autoscale(enable=True, axis='x', tight=True)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_ylabel('Torque [%MVC]')
ax2.set_xlabel('Time [s]')

plt.tight_layout()
plt.savefig('fig-2' + '.svg', dpi=300)
shutil.move('fig-2' + '.svg', os.path.join('.', 'data', 'proc', sub_info.sub, 'fig-2' + '.svg'))
plt.close()
