# jupyter
## Jupyter Notebooks for analyzing Horn Radio Astronomy Observations
### Glen Langston - National Scienc foundation - 2021 October 14

Here we provide a sequence of Aficionado/Student/Teacher/Hobbist projects, to provide a convenient method for learning how to use your Horn Radio Telescope observations.
Each of the Horn Telescopes use GnuRadio software designs to observe our Milky Way Galaxy.  Different versions of the software record different types measurements.

The horn radio telescopes are student, teacher and hobbiest built radio telescopes that can easily observe our place in the Milky Way Galaxy.   The horns are funnels a few feet in diameter, which channel radio wavelength signals to amplifiers and software defined radios. 
The calibration is based on Jupyter Notebooks.   Start with notebook <b>https://github.com/glangsto/jupyter/blob/master/MilkyWayPlotting.ipynb</b>.
<p>
The horn construction is described at www.WVURAIL.org
<p>
These notebooks assume the observations are taken using Gnuradio software systems available on GitHub.
Please see http://www.github.com/WVURAIL/gr-radio_astro for the software to make the observations.
These notebooks only work with NSF designs found there.
<p>
The Jupyter notebooks here work with the provided data.  You can edit these notebooks to work with your own observations.
<p>
  
  
### To Do:
1. Provide utilities to convert from and to Comma separated value files.
2. Provide a utility program to convert spectra taken with _spectrometer_w_cal.grc_ observations (in WVURAIL/gr-radio_astro/examples) to the .ast and .hot format files.
3. Divide ras.py into smaller units.
4. Expand the mapping code to spectral line cubes

