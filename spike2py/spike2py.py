import sys
import os
import numpy as np
import scipy.io as sio
from collections import namedtuple
from scipy.signal import butter, filtfilt


class Signal:
    """ Signal object

    Signal and accompanying information. Most often this is a signal recorded
    using Spike2, but can also be an external signal that needs to be
    processed alongside the signal recorded with Spike2.

    Parameters
    ----------
    sig : np.array or list or list of lists
        When used in the context of the `Trial` class and a signal recorded
        in Spike2, `sig` will be the signal obtained from the imported Matlab
        file. The signal is obtained from `data[s2name]`, where `s2name` is the
        name of the Spike2 channel. If using with a signal not recorded in
        Spike2, `sig` will be either the signal on its own, or the signal and an
        accompanying time axis. `sig` = [`signal`, `times`]  where `times` is
        optional.
    stype : str
        Signal type. `external`: signal not recorded by Spike2 but needed for analysis.
        `waveform` : signal sampled at regular intervals. `trig` : trigger signal
        `sEMG` : surface EMG signal
    s2name : str
        Spike2 channel name
    name : str
        Informative name for signal

    Other Parameters
    ----------------
    fs : int, optional
        Sampling frequency of signal. If `times` is not present, `fs` is used
        to create a time axis for the signal.
    filtInfo : FiltInfo, optional
        Contains filter information
    calibInfo : CalibInfo, optional
        Contains calibration information
    offsetInfo : OffsetInfo, optional
        Contains offset information
    normInfo : NormInfo, option
        Contains normalization information

    """

    def __init__(self, signal, stype, name, fs=None, s2name=None,
                 filtInfo=None, calibInfo=None, offsetInfo=None, normInfo=None):

        self.name = name
        self.stype = stype
        self.s2name = s2name
        self.fs = fs
        self.interp = dict()
        self.envel = None
        self.proc = None
        self.raw = None
        self.rect = None
        self.filtInfo = filtInfo
        self.calibInfo = calibInfo
        self.offsetInfo = offsetInfo
        self.normInfo = normInfo
        if not isinstance(signal, np.ndarray):
            signal = np.array(signal)
        self._add_raw_values(signal)
        self._add_times(signal)
        self._proc_raw_values()

    def __repr__(self):
        return(f"Signal(signal=signal, name={repr(self.name)}, s2name={repr(self.s2name)}, "
               f"stype={repr(self.stype)})")

    def _add_raw_values(self, signal):
        if self.stype is 'external':
            if signal.shape[0] == 2:
                self.raw = signal[0]
                self.proc = signal[0]
            else:
                self.raw = signal
                self.proc = signal
        elif self.stype is 'waveform' or self.stype is 'sEMG':
            self.raw = signal['values'][0][0].flatten()
            self.proc = signal['values'][0][0].flatten()
            if self.fs is None:
                try:
                    self.fs = int(1 / signal[0][0][2].flatten())
                except:
                    self.fs = None
        else:
            self.raw = None

    def _add_times(self, signal):
        if self.stype is 'trig':
            self.times = signal[0][0][4]

        elif self.stype is 'external':
            if signal.shape[0] == 2:
                self.times = signal[1]
                self.fs = 1 / (signal[1][1] - signal[1][0])
            else:
                self.times = self._new_times(self.fs, len(self.raw))
        elif self.stype is 'waveform' or 'sEMG':
            self.times = self._new_times(fs=self.fs, len_signal=len(self.raw))

    def _proc_raw_values(self):
        if self.stype is not 'trig':
            if self.calibInfo is not None:
                if self.calibInfo.slope is not None:
                    self.calibrate()
            if self.offsetInfo is not None:
                if self.offsetInfo.type is not None:
                    self.remove_offset()
            if self.filtInfo is not None:
                if self.filtInfo.cutoff is not None:
                    self.filter()
            if self.stype is 'sEMG':
                if self.offsetInfo is None:
                    self.proc = self.proc - np.mean(self.proc)
                self.filter()
                self.rect = abs(self.proc)
                self.envelop()
            if self.normInfo is not None:
                self.normalize()

    def _new_times(self, fs, len_signal):
        """ Generate new time axis.

        Parameters
        ----------
        fs : int or float
            Sampling frequency
        len_signal : int
            Number of samples for new time axis

        Returns
        -------
        times : np.array
            Newly created time axis

        """
        times = list()
        time = 0
        step = round(1 / fs, 4)
        for i in range(len_signal):
            times.append(time)
            time = round(time + step, 4)
        return np.array(times)

    def remove_offset(self, type_=None, val=None, sig='proc'):
        """ Remove offset from signal.

        Done when signal first created if `offsetInfo` parameter provided.
        However, it can be performed after the signal has been created if the
        `type` and `val` parameters are provided; `offsetInfo` will be added to
        signal.

        By default offset is removed from `proc` signal. However,
        if sig=`raw`, then the offset will be removed from the raw signal and
        the results stored in `proc`. if sig=`rect`, then the offset will be
        removed from the `rect` signals and the results stored in `rect`.

        See `OffsetInfo` for details on `type`, `val` and `sig`.

        Notes
        -----
        If an offset is removed from an EMG signal when first imported, then the
        mean of the `proc` version of the EMG signal will not be removed. It is
        assumed that the user has provided offset details that replicate this step
        or remove the offsets from only a portion of the trial (e.g. only the
        first 5s to avoid the influence of a stimulation artifact.

        """
        if type_ is not None:
            self.offsetInfo = OffsetInfo(type=type_, val=val, sig=sig)

        if self.offsetInfo.type is 'val':
            if self.offsetInfo.sig is 'proc':
                self.proc = self.proc - self.offsetInfo.val
            elif self.offsetInfo.sig is 'raw':
                self.proc = self.raw - self.offsetInfo.val
            elif self.offsetInfo.sig is 'rect':
                self.rect = self.rect - self.offsetInfo.val
            else:
                print(f'Offset signal {self.sig} not recognized.')
                sys.exit(1)

        elif self.offsetInfo.type is 'mean':
            if self.offsetInfo.sig is 'proc':
                self.proc = self.proc - np.mean(self.proc)
            elif self.offsetInfo.sig is 'raw':
                self.proc = self.raw - np.mean(self.raw)
            elif self.offsetInfo.sig is 'rect':
                self.rect = self.rect - np.mean(self.rect)
            else:
                print(f'Offset signal {self.sig} not recognized.')
                sys.exit(1)

        elif self.offsetInfo.type is 'start':
            if self.offsetInfo.sig is 'proc':
                self.proc = self.proc - np.mean(self.proc[0:self.offsetInfo.val])
            elif self.offsetInfo.sig is 'raw':
                self.proc = self.raw - np.mean(self.raw[0:self.offsetInfo.val])
            elif self.offsetInfo.sig is 'rect':
                self.rect = self.rect - np.mean(self.rect[0:self.offsetInfo.val])
            else:
                print(f'Offset signal {self.sig} not recognized.')
                sys.exit(1)
        else:
            print(f'Offset type {type_} not recognized.')
            sys.exit(1)

    def filter(self, cutoff=None, order=None, type_=None):
        """Filter signal.

        Done when signal first created if `filtInfo` parameter provided.
        However, it can also be performed after the signal has been created if
        the `cutoff`, `order` and `type` parameters are provided. `filtInfo`
        will be added to sig.

        See `FiltInfo` for details on `cutoff`, `order` and `type`.

        """
        if order is not None:
            N = order
            self.filtInfo = FiltInfo(cutoff, order, type_)
        else:
            N = self.filtInfo.order
        if cutoff is not None:
            Wn = np.array(cutoff)/(self.fs/2)
        else:
            Wn = np.array(self.filtInfo.cutoff) / (self.fs / 2)
        if type_ is not None:
            b, a = butter(N=N,
                          Wn=Wn,
                          btype=type_,
                          )
        else:
            b, a = butter(N=N,
                          Wn=Wn,
                          btype=self.filtInfo.type,
                          )
        self.proc = filtfilt(b, a, self.proc)

    def calibrate(self, slope=None, offset=None):
        """Calibrate signal.

        Done when signal first created if `calibInfo` parameter provided.
        However, it can also be performed after the signal has been created if
        the `slope` and `offset` parameters are provided. `calibInfo` will be
        added to sig.

        See `CalibInfo` for details on `slope` and `offset`.

        """
        if slope is None:
            slope = self.calibInfo.slope
            offset = self.calibInfo.offset
        else:
            self.calibInfo = CalibInfo(slope, offset)
        self.proc = (self.proc * slope) - offset

    def normalize(self, type_=None, value=None, signal_version=None):
        """Normalize signal.

        Done when signal first created if `normInfo` parameter provided.
        However, it can also be performed after the signal has been created.
        `normInfo` will be added to sig.

        See `normInfo` for details on `type`, `value` and `sig_version`.

        """
        if type_ is not None:
            if type(signal_version) is str:
                signal_version = [signal_version]
            self.normInfo = NormInfo(type_=type_,
                                     value=value,
                                     sig_version=signal_version)

        if self.normInfo.type == 'percentage':
            if self.normInfo.sig_version is None or \
                    'proc' in self.normInfo.sig_version:
                max_val = np.max(self.proc)
                self.proc = (self.proc / max_val)*100
            if 'raw' in self.normInfo.sig_version:
                max_val = np.max(self.raw)
                self.raw = (self.raw / max_val)*100
            if 'rect' in self.normInfo.sig_version:
                max_val = np.max(self.rect)
                self.rect = (self.rect / max_val)*100
            if 'envel' in self.normInfo.sig_version:
                max_val = np.max(self.envel)
                self.envel = (self.envel / max_val)*100
            if 'interp' in self.normInfo.sig_version:
                interp_sig = namedtuple('interp_sig', 'times fs vals')
                keys = list(self.interp.keys())
                for key in keys:
                    max_val = np.max(self.interp[key].vals)
                    vals = (self.interp[key].vals / max_val)*100
                    times = self.interp[key].times
                    fs = self.interp[key].fs
                    self.interp[key + 'norm_100'] = interp_sig(times=times,
                                                               fs=fs,
                                                               vals=vals)
        elif self.normInfo.type == 'proportion':
            if self.normInfo.sig_version is None or 'proc' in self.normInfo.sig_version:
                max_val = np.max(self.proc)
                self.proc = (self.proc / max_val)
            if 'raw' in self.normInfo.sig_version:
                max_val = np.max(self.raw)
                self.raw = (self.raw / max_val)
            if 'rect' in self.normInfo.sig_version:
                max_val = np.max(self.rect)
                self.rect = (self.rect / max_val)
            if 'envel' in self.normInfo.sig_version:
                max_val = np.max(self.envel)
                self.envel = (self.envel / max_val)
            if 'interp' in self.normInfo.sig_version:
                Interp_sig = namedtuple('Interp_sig', 'times fs vals')
                keys = list(self.interp.keys())
                for key in keys:
                    max_val = np.max(self.interp[key].vals)
                    vals = (self.interp[key].vals / max_val)
                    times = self.interp[key].times
                    fs = self.interp[key].fs
                    self.interp[key + 'norm_100'] = Interp_sig(times=times,
                                                               fs=fs,
                                                               vals=vals)
        elif self.normInfo.type == 'value':
            if self.normInfo.sig_version is None or 'proc' in self.normInfo.sig_version:
                max_val = self.normInfo.value
                self.proc = (self.proc / max_val)*100
            if 'raw' in self.normInfo.sig_version:
                max_val = self.normInfo.value
                self.raw = (self.raw / max_val)*100
            if 'rect' in self.normInfo.sig_version:
                max_val = self.normInfo.value
                self.rect = (self.rect / max_val)*100
            if 'envel' in self.normInfo.sig_version:
                max_val = self.normInfo.value
                self.envel = (self.envel / max_val)*100
            if 'interp' in self.normInfo.sig_version:
                interp_sig = namedtuple('interp_sig', 'times fs vals')
                keys = list(self.interp.keys())
                for key in keys:
                    max_val = self.normInfo.value
                    vals = (self.interp[key].vals / max_val)*100
                    times = self.interp[key].times
                    fs = self.interp[key].fs
                    self.interp[key + 'norm_100'] = interp_sig(times=times,
                                                               fs=fs,
                                                               vals=vals)
        else:
            print(f'Unrecognized Normalization type <{type_}>.')
            sys.exit(1)

    def envelop(self, cutoff=5):
        """Create a signal envelop for sEMG.

        Parameters
        ----------
        cutoff : int, default 5
            Cutoff frequency for lowpass filter
        """
        b, a = butter(N=4,
                      Wn=np.array(cutoff) / (self.fs / 2),
                      btype='lowpass',
                      )
        self.envel = filtfilt(b, a, self.rect)

    def interpolate(self, name, signal_id, subset=None, new_x=None, new_fs=None,
                    old_x=None):
        """ Linear interpolation of signal.


        Parameters
        ----------
        name : str
            name to give to new interpolated signal
        sig_id: str
            Select which values to interpolation (e.g. 'raw' or 'proc').
        subset : list of two int and one bool, optional
            Index values to subset signal and True/False of whether or not to
            zero new x-axis
        new_x: np.array (optional)
            If provided, directly used as the new_x to interp values
        new_fs: int (optional)
            If provided, new_x generated based on new sampling fq
        old_x: np.array (optional, default is self.times)
            If provided, used as the old_x

        Examples
        -------

        >>> import spike2py

        Interpolate to a new sampling rate

        >>> trial.sig['torque'].interpolate(name='fs15', signal_id='proc', new_fs=15)

        Interpolate to a new x-axis

        >>> new_x = list(range(len(trial.sig['torque'].raw_vals)))
        >>> trial.sig['torque'].interpolate(name='interpolate', signal_id='raw', new_x=new_x)

        Notes
        -----

        Can provide `new_x` to interpolate values to this new x-axis or times,
        or `fs` to interpolate to a new x-axis created with a sampling rate of
        `fs` Hz.

        The default x-axis for the original (pre-interpolation) signal is
        `self.times`; however, a different `old_x` can be specified.
        """
        if subset is None:
            subset = [0, -1, False]
        if signal_id is 'raw':
            old_y = self.raw[subset[0]: subset[1]]
        elif signal_id is 'proc':
            old_y = self.proc[subset[0]: subset[1]]
        else:
            print(f"Interpolate sign_id {signal_id} not recognized.")
            sys.exit(1)

        if old_x is not None:
            old_x = old_x
        else:
            old_x = self.times[subset[0]: subset[1]]

        if new_x is not None:
            new_x = new_x
            new_fs = 1/(new_x[1] - new_x[0])
        elif new_fs is not None:
            new_x = np.arange(self.times[(subset[0])], self.times[(subset[1])],
                              1/new_fs)

        if subset[2]:
            old_x = old_x - old_x[0]
            new_x = new_x - new_x[0]

        new_y = np.interp(x=new_x, xp=old_x, fp=old_y)
        Interp_sig = namedtuple('Interp_sig', 'times fs vals')
        self.interp[name] = Interp_sig(times=new_x, fs=new_fs, vals=new_y)


class Trial:
    """Class for trial.

    Parameters
    ----------
    trialInfo : TrialInfo
        Class containing details about trial and its signals

    Examples
    --------

    >>> import spike2py

    >>> trial_pc1 = spike2py.Trial(trial_info1)  # Assumes `trial_info1` available

    >>> trial_pc2 = spike2py.Trial(trial_info2)  # Assumes `trial_info2 available
    >>> trial_pc1.merge('pc2', trial_pc2)

    >>> us = spike2py.SigInfo(name='us', stype='external', fs=15)  # ultrasound
    >>> us_sig = spike2py.Signal(us_dat, stype=us.stype, name=us.name, fs=us.fs)
    >>> trial_pc1.add_sig('us', us_sig)

    """

    def __init__(self, trialInfo):
        self.cond = trialInfo.cond
        self.path = trialInfo.path
        self.filename = trialInfo.filename
        self.infoSignals = trialInfo.signals
        self.sig = dict()
        self._read_sig()

    def __repr__(self):

        txt = (f"Trial(path={repr(self.path)}, "
              f"filename={repr(self.filename)}, "
              f"cond={repr(self.cond)}, "
              f"sig=sig)"
              f"\n# signals in trial:")
        for sig in self.sig.keys():
            txt = txt + f"\n#\t {sig}"
        return txt

    def _read_sig(self):
        file = os.path.join(self.path, self.filename)
        data = sio.loadmat(file)
        for sig in self.infoSignals:
            self.sig[sig.name] = Signal(signal=data[sig.s2name],
                                        stype=sig.stype,
                                        name=sig.name,
                                        fs=sig.fs,
                                        s2name=sig.s2name,
                                        filtInfo=sig.filtInfo,
                                        calibInfo=sig.calibInfo,
                                        offsetInfo=sig.offsetInfo,
                                        )

    def add_sig(self, sig_name, sig):
        """Add signal.

        Parameters
        ----------
        sig_name : str
            Short informative name for signal; will be key to retrieve
            signal from `Trial.sig[key]`
        sig : `Signal`, list, np.array, int, float, bool
            Data associated with added signal

        Notes
        -----

        Designed to have `sig` be a `Signal` instance, but this implementation
        provides flexibility to the user to add any type of signal or value.

        """
        self.sig[sig_name] = sig

    def merge(self, trial, prefix=None):
        """Merge trials.

        Adds all signals from another `Trial` instance.

        Parameters
        ----------
        trial: Trial
            Trial containing signals to be merged
        prefix : str, optional
            Prefix to add to merged signal names. Default is `merged\_`

        """
        if prefix is None:
            prefix = 'merged_'
        else:
            prefix = prefix + '_'
        for key, vals in trial.sig.items():
            self.add_sig(prefix + key, vals)


class SigInfo:
    """ Information to process signal recorded with Spike2.

    Passed to spike2.TrialInfo.sig, which is a list of all sig,
    or a single signal.

    Parameters
    ----------
    name : str
        User defined name, because Spike2 channel names are not always clear.
    stype : str
        Signal type. `external`for a signal not recorded by Spike2 but needed
        for analysis. `waveform` for a signal sampled at regular intervals.
        `trig` for a trigger signal. `sEMG` for a surface EMG signal.
    s2name : str
        Spike2 channel name.

    Other Parameters
    ----------------
    fs : int, optional
        Sampling frequency of signal; used when `stype` is `external`
    filt_cutoff : int, optional
        Value of filter
        Either single value (e.g. 30) or list of 2 values (e.g. [20, 450])
    filt_order: int, optional
        Order of filter
    filt_type : str, optional
        Define type of filter: `lowpass`, `highpass`, `bandpass`, `bandstop`
    calib_slope : float or int, optional
        Slope of the calibration equation
    calib_offset : float or int, optional
        Offset of the calibration equation (i.e. intercept)
    offset_type : str, optional
        `val` to remove value provided. `mean` to remove mean of signal.
        `start` to remove mean of n points from start of trial.
    offset_val: int, optional
        Value required for `val` and `start` offset removal types.
    offset_sig: str, default 'proc'
        Specify which signal from which to remove offset
        offset removal types.
    norm_type: str, optional
        Specify type of normalization to apply. Possibilities are
        `percentage`: Signal will range from 0-100, based on maximum in signal
        `proportional`: Signal will range from 0-1, based on maximum in signal
        `value`: Signal will range from 0-100, based on maximum in signal
    norm_val: float, optional
        Value required for `value` normalization type.
    norm_sig_version: str or list of str, optional
        Possibilities include `raw`, `_extract_data`, `rect`, `envel`, `interp`.
        Default is to use `_extract_data`. `interp` will apply normilization to all `interp` signals

    Examples
    --------

    >>> import spike2py

    >>> gas = spike2py.SigInfo(name='gas', stype='sEMG', s2name='GAS',
						filt_cutoff=5, filt_order=4, filt_type='lowpass',
						calib_slope=0.5, calib_offset=55,
						offset_type='val', offset_val=22, offset_sig='proc'
					   norm_type='percentage', norm_sig_version=['raw', '_extract_data'])
    >>> gas
    SigInfo(name='gas',
		stype='sEMG',
		s2name='GAS',
		filt_cutoff=[20, 450],
		filt_order=4,
		filt_type='bandpass',
		calib_slope=0.5,
		calib_offset=55,
		offset_type='val',
		offset_val=22,
		offset_sig='proc',
		norm_type='percentage',
		norm_val=None
		 norm_sig_version=['raw', '_extract_data'])
    """

    def __init__(self, name, stype, s2name=None, fs=None,
                 filt_cutoff=None, filt_order=None, filt_type=None,
                 calib_slope=None, calib_offset=None,
                 offset_type=None, offset_val=None, offset_sig='proc',
                 norm_type=None, norm_val=None, norm_sig_version=None):

        self.name = name
        self.stype = stype
        self.s2name = s2name
        self.fs = fs


        if filt_cutoff is not None:
            self.filtInfo = FiltInfo(stype, filt_cutoff, filt_order, filt_type)
        elif stype is 'sEMG':
            self.filtInfo = FiltInfo(stype)
        else:
            self.filtInfo = FiltInfo()

        if calib_slope is not None:
            self.calibInfo = CalibInfo(calib_slope, calib_offset)
        else:
            self.calibInfo = CalibInfo()

        if offset_type is not None:
            self.offsetInfo = OffsetInfo(offset_type, offset_val, offset_sig)
        else:
            self.offsetInfo = OffsetInfo()
        if norm_type is not None:
            self.normInfo = NormInfo(type_=norm_type,
                                     value=norm_val,
                                     sig_version=norm_sig_version)
        else:
            self.normInfo= NormInfo()

    def __repr__(self):
        return (f"SigInfo(name={repr(self.name)},\n\t\tstype={repr(self.stype)},\n\t\t"
                f"s2name={repr(self.s2name)},\n\t\tfilt_cutoff={repr(self.filtInfo.cutoff)},\n\t\t"
                f"filt_order={repr(self.filtInfo.order)},\n\t\tfilt_type={repr(self.filtInfo.type)},\n\t\t"
                f"calib_slope={repr(self.calibInfo.slope)},\n\t\tcalib_offset={repr(self.calibInfo.offset)},\n\t\t"
                f"offset_type={repr(self.offsetInfo.type)},\n\t\toffset_val={repr(self.offsetInfo.val)},\n\t\toffset_sig={repr(self.offsetInfo.sig)}\n\t\t"
                f"norm_type={repr(self.normInfo.type)},\n\t\tnorm_val={repr(self.normInfo.value)}\n\t\t norm_sig_version={repr(self.normInfo.sig_version)})")


class TrialInfo:
    """ Information to process trial recorded with Spike2.

    Passed to spike2.Trial(). Assumes the data has been exported using
    the Spike2 script `Batch export MATLAB.2s2` available from CED.

    Parameters
    ----------
    cond : str
        Experimental condition (e.g. 'cs-100mvc' or 'grasp_all_digits')
    path : str
        path to .mat file exported from Spike2
    filename : str
        Full name of .mat file (e.g. 'cs-100mvc.mat')
    signals : list of SignalInfo
        Each item in list is a `SigInfo` instance.
        These are the sig to be analysed.

    Examples
    --------

    >>> import spike2py

    >>> cond = 'grasp 10s'
    >>> path = '/home/martin/Documents/grasp_study/data/'
    >>> filename = 'grasp_10s.mat'
    >>> signal = spike2py.SigInfo(name='FDI EMG', stype='sEMG', s2name='FDI_emg')
    >>> trial_info = spike2py.TrialInfo(cond, path, filename, signal)
    >>> trial_info
    spike2.TrialInfo(cond='grasp 10s',
		path='/home/martin/Documents/grasp_study/data/',
		filename='grasp_10s.mat',
		signal=signal)
		# signal is an instance of spike2.SigInfo or a list of instances.

    """

    def __init__(self, cond, path, filename, signals):
        self.cond = cond
        self.path = path
        self.filename = filename
        self.signals = signals

    def __repr__(self):
        return f"spike2.TrialInfo(cond={repr(self.cond)}, \n\t\tpath={repr(self.path)}," \
               f"\n\t\tfilename={repr(self.filename)}, \n\t\tsignal=signal)" \
               f"\n\n\t\t# signal is an instance of spike2.SigInfo or a list of instances."


class FiltInfo:
    """ Filter information.

    Can be used independently, but most often passed to `SigInfo`.

    Parameters
    ----------
    stype : str, optional
        Signal type. `external`, `waveform`, `sEMG`
    cutoff : int, float, or list, optional
        Cutoff frequency of filter. Can be a single value for `lowpass` and
        `highpass` filters (e.g. 5), or a pair of `int` or `float` in a list
        (e.g. [20, 450]).
    order : int, optional
        Filter order; usually 2, 4 or 8.
    type : str, optional
        `lowpass`, `highpass`, `bandpass`, `bandstop`

    Examples
    --------

    >>> import spike2py

    >>> filter = spike2py.FiltInfo(stype='sEMG')
    >>> filter
    FiltInfo(cutoff=[20, 450], order=4, type='bandpass')


    >>> filter = spike2py.FiltInfo(cutoff=50, order=4, type='highpass')
    >>> filter
    FiltInfo(cutoff=50, order=4, type='highpass')


    >>> filter = spike2py.FiltInfo(cutoff=[2, 10], order=8, type='bandstop')
    >>> filter
    FiltInfo('cutoff=[2, 10], order=8, type='bandstop')

    Notes
    -----
    The default for `stype` sEMG assumes the filter will be used on a raw
    signal. When a sEMG signal is created, a rectified and 5Hz lowpass
    envelop version is also created. If an envelop with a different cutoff
    frequency is desired, run `Signal.envelop(cutoff=val)` with the
    desired cutoff value specified.

    """

    def __init__(self, stype=None, cutoff=None, order=None, type=None):
        if stype is None:
            self.cutoff = cutoff
            self.order = order
            self.type = type
        elif stype is 'sEMG':
            self.cutoff = [20, 450]
            self.order = 4
            self.type = 'bandpass'
        else:
            self.cutoff = cutoff
            self.order = order
            self.type = type

    def __repr__(self):
        return f"FiltInfo(cutoff={repr(self.cutoff)}, order={repr(self.order)}, " \
               f"type={repr(self.type)})"


class CalibInfo:
    """ Calibration information.

    Can be used independently, but most often passed to `SigInfo`.

    Parameters
    ----------
    slope : int or float, optional
        Slope of calibration equation
    offset : int or float, optional
        Offset of calibration equation (i.e. intercept)

    Examples
    --------

    >>> import spike2py

    >>> calibrate = spike2py.CalibInfo(slope=54.8, offset=2.6)
    >>> calibrate
    CalibInfo(slope=54.8, offset=2.6)

    """

    def __init__(self, slope=None, offset=None):
        self.slope = slope
        self.offset = offset

    def __repr__(self):
        return f'CalibInfo(slope={repr(self.slope)}, offset={repr(self.offset)})'


class OffsetInfo:
    """ Offset removal information.

    Can be used independently, but most often passed to `SigInfo`.

    Parameters
    ----------
    type : str, optional
        `val` to remove value provided. `mean` to remove mean of signal.
        `start` to remove mean of n points from start of trial.
    val: int, optional
        Value required for `val` and `start` offset removal types.

    sig: str, default (`proc`)
        Identifies which signal to remove the offset. `proc` and `rect` remove
        the offset in-place (i.e. result will overwrite original signal.
        However, `raw` will cause the offset to be removed from the `raw` signal
        and the result stored in `proc`.


    Examples
    --------

    >>> import spike2py

    >>> offset = spike2py.OffsetInfo(type='val', val=44.5)
    >>> offset
    OffsetInfo(type='val', val=44.5, sig='proc')


    >>> offset = spike2py.OffsetInfo(type='mean', val='')
    >>> offset
    OffsetInfo(type='val', val='', sig='proc')


    >>> offset = spike2py.OffsetInfo(type='start', val=100, sig='raw')
    >>> offset
    OffsetInfo(type='start', val=100, sig='raw')

    >>> offset = spike2py.OffsetInfo(type='start', val=100, sig='rect')
    >>> offset
    OffsetInfo(type='start', val=100, sig='rect')

    Notes
    -----
    If an offset is removed from an EMG signal when first imported, then the
    mean of the `proc` version of the EMG signal will not be removed. It is
    assumed that the user has provided offset details that replicate this step
    or remove the offsets from only a portion of the trial (e.g. only the
    first 5s to avoid the influence of a stimulation artifact.

    """

    def __init__(self, type=None, val=None, sig='proc'):
        self.type = type
        self.val = val
        self.sig = sig

    def __repr__(self):
        return f"OffsetInfo(type={repr(self.type)}, val={repr(self.val)}, " \
            f"sig={repr(self.sig)})"


class NormInfo:
    """ Normalization information.

    Can be used independently, but most often passed to `SigInfo`.

    Parameters
    ----------
    type: str
        Specify type of normalization to apply. Possibilities are
        `percentage`: Signal will range from 0-100, based on maximum in signal
        `proportional`: Signal will range from 0-1, based on maximum in signal
        `value`: Signal will range from 0-100, based on maximum in signal
    val: float, optional
        Value required for `value` normalization type.
    sig_version : list of str items, optional
        Possibilities include `raw`, `_extract_data`, `rect`, `envel`, `interp`.
        Default is to use `proc`. `interp` will apply
        normilization to all `interp` signals


    Examples
    --------

    >>> import spike2py

    >>> norm_info = spike2py.NormInfo(type_='proportion')
    >>> norm_info
    NormtInfo(type='proportion', value=None, sig_version=['_extract_data'])


    >>> norm_info =spike2py.NormInfo(type_='percentatge', sig_version=['raw', 'proc', 'rect', 'envel', 'interp'])
    >>> norm_info
    NormtInfo(type='percentatge', value=None, sig_version=['raw', 'proc', 'rect', 'envel', 'interp'])

    >>> norm_info =spike2py.NormInfo(type_='value', value=4.53)
    >>> norm_info
    NormtInfo(type='value', value=4.53, sig_version= ['proc'])

    """

    def __init__(self, type_=None, value=None, sig_version=None):
        self.type = type_
        self.value = value
        if type(sig_version) is str:
            self.sig_version = [sig_version]
        elif type(sig_version) is list:
            self.sig_version = sig_version
        elif sig_version is None:
            self.sig_version = ['proc']
        else:
            print(f'The parameter sig_version for the NormInfo class'
                  f'was of an unexpected type <{type(sig_version)}>')
            print('Expect a str or a list of str.')

    def __repr__(self):
        return f"NormtInfo(type={repr(self.type)}, value={repr(self.value)}, " \
            f"sig_version={repr(self.sig_version)})"
