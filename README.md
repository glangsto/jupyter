# jupyter
## Jupyter Notebooks for analyzing Horn Radio Astronomy Observations
### Glen Langston - National Scienc foundation - 2021 September 17
This project, a work in progress, is intended to provide a convenient method of checking Horn Radio Telescope observations.
The horn radio telescopes are student, teacher and hobbiest built radio telescopes that can easily observe our place in the Milky Way Galaxy.   The horns are funnels a few feet in diameter, which channel radio wavelength signals to amplifiers and software defined radios. 
The calibration is based on Jupyter Notebooks.   Start with notebook <b>_plotting.ipynb_</b>.
<p>
The horn construction is described at www.WVURAIL.org
<p>
These notebooks assume the observations are taken using Gnuradio software systems available on GitHub.
Please see http://www.github.com/WVURAIL/gr-radio_astro for the software to make the observations.
These notebooks only work with NSF designs found there.
<p>
The Jupyter notebooks found here are only for analyzing already completed observations.
<p>
### To Do:
1. Provide utilities to convert from and to Comma separated value files.
2. Provide a utility program to convert spectra taken with _spectrometer_w_cal.grc_ observations (in WVURAIL/gr-radio_astro/examples) to the .ast and .hot format files.
3. Divide ras.py into smaller units.
4. Expand the mapping code to spectral line cubes

