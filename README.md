# jupyter
## Jupyter Notebooks for analyzing Horn Radio Astronomy Observations
### Glen Langston - National Scienc foundation - 2022 January 19

Here we provide a sequence of Aficionado/Student/Teacher/Hobbist projects.
These plotting and computing tools guide you in looking at your own ombservations.

You can update these tools to make discoveries with your own Horn Radio Telescopes.
Each of the Horn Telescopes use GnuRadio software designs to observe our Milky Way Galaxy.  Different versions of the software record different types measurements.
The horn radio telescopes are student, teacher and hobbiest built radio telescopes that can easily observe our place in the Milky Way Galaxy.   The horns are funnels a few feet in diameter, which channel radio wavelength signals to amplifiers and software defined radios. 
<p>
The horn construction is described at www.WVURAIL.org
<p>
These notebooks assume the observations are taken using Gnuradio software systems available on GitHub.
Please see http://www.github.com/WVURAIL/gr-radio_astro for the software to make the observations.
These notebooks only work with NSF designs found there.
<p>
<h3> A notebook for each topic</h3>
The Jupyter notebooks here work with the provided data.  You can edit these notebooks to work with your own observations.
We intend to have a notebook for each of our "Lessons".   Here are the first lessons
<p><ol>
<li>Chart the Shape of the Milky Way:  <b>https://github.com/glangsto/jupyter/blob/master/MilkyWay1Intensity.ipynb</b>.
<li>Hydrogen, its everywhere!:  A days observations: <b>https://github.com/glangsto/jupyter/blob/master/MilkyWayPlotting.ipynb</b>.
</ol>
<p>
<p>
  
  
### To Do:
1. Provide utilities to convert from and to Comma separated value files.
2. Provide a utility program to convert spectra taken with _spectrometer_w_cal.grc_ observations (in WVURAIL/gr-radio_astro/examples) to the .ast and .hot format files.
3. Divide ras.py into smaller units.
4. Expand the mapping code to spectral line cubes

