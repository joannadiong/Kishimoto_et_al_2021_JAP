Overview
========

`Spike2`_ is the software used by many scientists to collect data using a
`Cambridge Electronics Design (CED)`_ data aquisition boards. While some
scientists use Spike2 to analyse their data, other scientists prefer to export
their data to other programs such as Matlab or Python.

There are ways to directly open Spike2 data files in Python. For example,
`Neo`_ is a package for analysing electrophysiology data that supports reading
a wide range of neurophysiology file formats, including Spike2 files.
However, package updates have been known to cause `issues`_ with reading Spike2
files. Moreover, the Neo framework was specifically designed to handle data from
electrophysiological experiments (cellular and animal recordings; it is less well
suited, or simply overkill, for data collected for life sciences, biomechanics,
and human neurophysiology experiments.

As an alternative, Spike2 data can be exported to Matlab data files (.mat). While
this requires the user to export data prior to opening it in Python, this
approach has some advantages. First, CED makes available a Spike2 script that
can batch export all Spike2 data files contained within a folder. A copy of this
script file is available `here`_. Second, opening Matlab data files has been
supported for a long time in Python through the `scipy.io`_ module.

The one downside to this approach is that when these Matlab data files are
opened in Python, it is not intuitive how the data are organized, nor where
the various details of your signal might be stored (e.g. sampling rate, sample
times, etc). Rarely would someone ever work directly with the default structure
of nested numpy arrays in a dictionary.

The present module provides a simple interface to signals and triggers recorded
in Spike2. The main building blocks are a :class:`~spike2py.Signal` class and
a :class:`~spike2py.Trial` class. The :class:`~spike2py.Trial` class can store
all the signals recorded during a trial, as well
as other details related to that trial. Basic signal processessing, such as
calibration, filtering, offset removal, can be applied to the signals. There is
also a default routine for surface electromyography (EMG) recordings (filtering,
rectifying, envelop filtering). The module include several helper classes to
streamline how :class:`~spike2py.Trial` and :class:`~spike2py.Signal` are accessed.

.. _Spike2: http://ced.co.uk/products/spkovin
.. _Cambridge Electronics Design (CED): http://ced.co.uk/
.. _Neo: https://github.com/NeuralEnsemble/python-neo
.. _`issues`: https://scientificallysound.org/2018/04/05/import-spike2-into-python/
.. _`here`: https://github.com/MartinHeroux/spike2py/blob/master/tutorials/Batch_export_MATLAB.s2s
.. _`scipy.io`: https://docs.scipy.org/doc/scipy/reference/io.html