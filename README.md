# Estimation of maximal muscle electromyographic activity from the relationship between muscle activity and voluntary activation

Kenzo C Kishimoto<sup>1</sup>
Martin E Héroux<sup>2,3</sup>
Simon C Gandevia<sup>2,3</sup>
Jane E Butler<sup>2,3</sup>
Joanna Diong<sup>2,4</sup>

<sup>1</sup>Discipline of Physiotherapy, Faculty of Health Sciences, The University of Sydney, New South Wales, Australia.  
<sup>2</sup>Neuroscience Research Australia (NeuRA), Sydney, New South Wales, Australia  
<sup>3</sup>School of Medical Sciences, University of New South Wales, New South Wales, Australia  
<sup>4</sup>Discipline of Anatomy and Histology, School of Medical Sciences, Faculty of Medicine and Health, The University of Sydney, New South Wales, Australia

## Suggested citation

Kishimoto KC, Héroux ME, Gandevia SC, Butler JE, Diong J (2021) Estimation of maximal muscle electromyographic activity from the relationship between muscle activity and voluntary activation. Journal of Applied Physiology (in press).

## Data

Raw data to generate Fig 2 and 3 are stored in on the Open Science Framework (OSF) repository [https://osf.io/wt7z8/][osf] and in **data -> raw** in these formats (in zipped folders):
    Spike2 .smr 
    Matlab .mat 
    Text .txt

Unzip the folders and extract the folders **sub01** and **sub13** into **data -> raw** in the project folder.

Processed data for statistical analysis and Stata code files are stored in the **stats** folder in these formats: 
    Text .csv
    Text .txt
    Stata .do

**stats -> data -> subjects_data.csv**
The Stata .do file analyses processed data in this dataset.

## Code

Python code files were written by Joanna Diong with input from Martin Héroux (Python v3.7). 

Stata code files were written by Joanna Diong (Stata 16, compatible with v13 and up). 

### Python

Python code is run in the **activation** virtual environment that contains the dependency files. The environment can be installed from the Terminal (Mac or Linux, or using the PyCharm Terminal in Windows) using `activation.yml`. To install, type in Terminal:

  `conda env create -f activation.yml`

To activate the environment, type: 

  `conda activate activation`

To deactivate, type:

  `conda deactivate`

However, the best way to reproduce the Python analysis would be to run the code in an integrated development environment (e.g. [PyCharm][pycharm]). Set the interpreter to operate in the virtual environment within the project folder.

**process.py** calls the deprecated Python package **spike2py** written by Martin Héroux. Download the package to a location outside of the project folder, and point (i.e. add root) the Python interpreter towards that location.

Note, for the new and revised packaged, see the [spike2py GitHub page][spike2py].

**fig-2.py**: Generate SVG file of single participant trial data (Fig 2) in **data -> proc -> sub01**

**fig-3.py**: Generate SVG file of single participant trial data of all signals (Fig 3) in **data -> proc -> sub13**

**fig-7.py**: Generate SVG file of predicted on observed EMG data (Fig 7) in **data -> proc**

**process.py, utilities.py, trials_key.py**: Modules containing functions used to clean data and plot figures. 

### Running Python code

Download all code files, data.zip, and stats.zip into a single folder. Unzip the data file into the same location.
Download the **spike2py** package. Point the Python interpreter to the location of **spike2py**.
Set file path for `REPO` to local directory:

1. **process.py** line 18
2. **fig-7.py** line 11

Run each file separately: **fig-2.py**, **fig-3.py**, **fig-7.py**

### Stata

Running **stats -> script.do** imports **subjects_data.csv** and models activation and EMG. Output and locations: 

1. **stats -> log_files**: log files of results
2. **stats -> graphs**: figures

Fig 4, 5 and 6 are generated from the Stata do file.

### Running Stata code

Open **stats -> script.do**.

File paths to generate results and figures are currently configured for Linux or Mac using forward single slashes `/`. To configure for Windows, use double backslashes `\\` (**script.do**, lines 28, 35-37). 

In line 28 of **script.do**, set file path to local directory.

Uncomment line 22 of `ssc install mlt` to install `mlt` package, unless already installed.

Uncomment lines 47-48, 417 to save a log file in **stats -> log_files**.

Run **script.do**.

[osf]: https://osf.io/wt7z8/ 
[spike2py]: https://github.com/MartinHeroux/spike2py
[pycharm]: https://www.jetbrains.com/pycharm/promo/?gclid=Cj0KCQiAtqL-BRC0ARIsAF4K3WFahh-pzcvf6kmWnmuONEZxi544-Ty-UUqKa4EelnOxa5pAC9C4_d4aAisxEALw_wcB 
