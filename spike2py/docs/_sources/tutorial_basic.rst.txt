Basic Example
-------------
We will demonstate the basic usage of :mod:`spike2py` using data from a simple
experiment. The data was collected in Spike2 and then exported as a Matlab data
file using `Batch_export_MATLAB.s2s`_. The trial lasted approximately
10s and included three signals: a trigger signal, an EMG signal, and a pressure
signal.

The trigger signal was recorded as a `trigger`. This means Spike2 saved the times
when the trigger pulse occurred. The EMG and pressure signal were both recorded
as `waveform` signals, which is what Spike2 calls continues signals sampled at
a given sampling rate.

Here is the basic structure to retrieve the data and obtain a :class:`~spike2py.Trial`
object. :class:`~spike2py.SigInfo` will be used in its most basic form to give
:mod:`spike2py` the required information to properly import our signals. We
will also use :class:`~spike2py.Trial` to contain all the data and information
about our trial, including the signals themselves. ::

    import spike2py

    exp_cond = 'maximum pressure trial 1'
    filename = 'max_push1.mat'

    trig_info = spike2py.SigInfo(name='trig',
                                 stype='trig',
                                 s2name='trig',
                                 )

    emg_info = spike2py.SigInfo(name='biceps',
                                stype='sEMG',
                                s2name='EMG',
                                )

    pressure_info = spike2py.SigInfo(name='pressure',
                                     stype='waveform',
                                     s2name='signal',
                                     )

    signals = [trig_info, emg_info, pressure_info]

    trial_info = spike2py.TrialInfo(cond=exp_cond,
                                    path='./tutorials/',
                                    filename=filename,
                                    signals=signals,
                                    )

    trial = spike2py.Trial(trial_info)

We now have a `trial` object that contains information about our trial, and the
various signals we recorded.

>>> trial.cond
'maximum pressure trial 1'
>>> trial.filename
'max_push1.mat'
>>> trial.sig
{'trig': Signal(sig=sig, name='trig_info', s2name='trig', stype='trig'),
 'biceps': Signal(sig=sig, name='biceps', s2name='EMG', stype='sEMG'),
 'pressure': Signal(sig=sig, name='pressure', s2name='signal', stype='waveform')}

We can access our data to process or plot. For example, we can inspect the
trigger times.

>>> trial.sig['trig'].times
array([[0.701],
       [3.501],
       [3.801],
       [4.701],
       [4.901],
       [5.101],
       [5.301],
       [5.501],
       [5.701],
       [5.901],
       [6.101],
       [6.301],
       [6.501]])

We can plot the biceps EMG data. Note that by specifying `stype=sEMG` for our
EMG data, :mod:`spike2py` automatically creates a rectified and envelop version
of our EMG data.::

    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,4), dpi=300)
    plt.subplot(121)
    plt.plot(trial.sig['biceps'].times, trial.sig['biceps'].raw, 'k', label='raw EMG')
    plt.ylabel('amplitude (a.u.)')
    plt.xlabel('time (s)')
    plt.legend()
    plt.subplot(122)
    plt.plot(trial.sig['biceps'].times, trial.sig['biceps'].rect, label='rectified EMG')
    plt.plot(trial.sig['biceps'].times, trial.sig['biceps'].envel, 'r', label='envelop EMG')
    plt.ylabel('amplitude (a.u.)')
    plt.xlabel('time (s)')
    plt.legend()

.. image:: ./biceps_emg.png

We can also plot the pressure.

>>> plt.figure(figsize=(6, 4), dpi=300)
>>> plt.plot(trial.sig['pressure'].times, trial.sig['pressure'].raw, 'k')
>>> plt.ylabel('amplitude (a.u.)')
>>> plt.xlabel('time (s)')

.. image:: ./pressure.png
   :height: 400
   :width: 600

.. _`Batch_export_MATLAB.s2s`: https://github.com/MartinHeroux/spike2py/blob/master/tutorials/Batch_export_MATLAB.s2s
