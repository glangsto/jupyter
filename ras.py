#Python class to plot raw NSF spectra.
#HISTORY
#21Sep16 GIL implement averaging in raw() start on tsys()
#21Sep15 GIL updates for python 3 and add arguements for tsys plotting
#21Jun07 GIL rr -> ras, complete class
#21Jun06 GIL initial version 
#
import matplotlib as mpl
import numpy as np
import copy
import radioastronomy
import interpolate
import gainfactor as gf
import rasnames

class Plot(object):
    """
    Define the different types of plots for Radio Astronomy Service (RAS) Horns
    """
    
    def setdefaults( self):
        """
        Set (and reset) default plot parameters
        """
        self.plotFrequency = True,
        self.doBaseline = False
        self.writeTsys = False
        self.writeKelvin = False
        self.doSave = False
        self.flagRfi = False
        self.firstdate = ""
        self.lastdate = ""
        self.nave = 0
        self.minGlat = 0.
        self.maxGlat = 0.
        self.minel = 0.
        self.maxel = 0.
        self.lastel = 999.
        self.lastaz = 999.
        self.lastgain = 0.
        self.beginutc = 0.
        self.endutc = 0.
        # end of setdefaults()
        return

#   some SDRs put spike in center of spectrum; indicate spike flagging here
    def __init__( self, verbose=False,  
                  flagCenter = False, 
                  doSave = False, 
                  flagRfi = True, 
                  plotFrequency = True, 
                  myTitle = "",
                  doPlotFile = False, 
                  outFileDir = "../out",
                  plotFileDir = "../plot",
                  keepFileDir = "../keep",
                  aveTimeSec=3600.,
                  telIndex = 0,
                  doBaseline = False,
                  writeTsys = False,
                  writeKelvin = False,
                  ):
        self.verbose = verbose          # print extra messages
        self.flagCenter = flagCenter    # if flag spike in center of spectrum,
        self.doSave = doSave            # if saving intermediate/average files
        self.flagRfi = flagRfi          # flag flagging RFI
        self.plotFrequency = plotFrequency   # plot frequency not velocity
        self.myTitle = myTitle          # provide plot title
        self.doPlotFile = doPlotFile    # if creating a plot file
        self.plotFileDir = plotFileDir  # location of plot Files
        self.outFileDir = outFileDir    # location of intermediate files
        self.keepFileDir = keepFileDir  # location to keep calibration (hot,cold) files
        self.aveTimeSec = aveTimeSec    # averaging time for spectra
        self.telIndex = telIndex        # add telescope number to plots
        self.doBaseline = doBaseline    # optionally subtract a baseline
        self.doZero = False             # optionally add a zero line
        self.writeTsys = writeTsys      # write log if integrated intenities
        self.writeKelvin = writeKelvin  # write calibrated intensities
        self.doKeep = False             # optionally write hot + cold raw files
        self.lowel = 10.                # minimum elevation for cold obs.
        
        # constants
        self.EPSILON = 1.E-10           # define a small number
        # put your list of known RFI features here.
        #Must have at least two, if flagRfi is true.
        self.linelist = [1400.00, 1420.0]  # RFI lines in MHz
        self.linewidth = [ 5, 5] # number of channels to flag

        # define reference frequency for velocities (MHz)
        self.nuh1 = 1420.40575 # neutral hydrogen frequency (MHz)
        self.nuoh1= 1612.231   # OH line
        self.nuoh2= 1665.402   # OH line
        self.nuoh3= 1667.359   # OH Line
        self.nuoh4= 1720.530   # OH Line

        # select the default frequency for plotting velocities
        self.nuRefFreq = self.nuh1
        self.c = 299792.458  # (Speed of light  km/sec)

        self.myTitle = myTitle    # input plot Title

        # set default names for files for hot and cold loads
        # null means compute 
        self.hotFileName = ""
        self.coldFileName = ""

        # currently used velocities for plotting range
        self.maxvel = 220.
        self.minvel = -220.
        self.maxSvel = 150.   # max velocity for numerical integration 
        self.minSvel = - self.maxSvel # min velocity for numerical integration 
        self.fitOrder = 3
        self.maxPlot = int(25)
        self.fileTag = ""
        self.firstdate = ""
        self.lastdate = ""
        self.doScaleAve = False
        self.xa = -1     # initialize indices corresponding to bandwidth to plot
        self.xb = -1
        self.scalefactor = 1.0
        self.nplot = 0   # count of plots so far
        # initialize of range for search
        self.setdefaults()
        self.lowel = 30.
        self.lowGlat = 30.

        # create places to store hot and cold files
        self.ave_hot  = radioastronomy.Spectrum()
        self.ave_cold = radioastronomy.Spectrum()
        self.nhot = 0
        self.ncold = 0
        self.thot = 295. # Kelvins
        self.tcold = 10. # Kelvins
        self.nData = 128
        self.hv = np.zeros(self.nData)
        self.cv = np.zeros(self.nData)
        self.gain = np.zeros(self.nData)
        self.vel = np.zeros(self.nData)
        self.tRxMiddle = 100. # kelvins, estimated
        self.gainAve = 1.
        self.tint = 30.
        
        # parameters for determining a change of obs
        self.nave = 0
        self.lastel = 999.
        self.lastaz = 999.
        self.lastgain = 0.
        self.beginutc = 0.
        self.endutc = 0.
        # store new spectra here if a change is found
        self.next = radioastronomy.Spectrum()
        self.next.tint = 0.   # indicate no data yet
        
        # end of init
        return

    def help(self, argstring=""):
        ###
        # provide user help info also parses and updates configuration
        ###
        args = argstring.split()
        nargs = len(args)

        if (nargs < 1): 
            print("ras.help(flags): Plotting Inputs for telescope obs.")
            print("Usage: .help('<flags> <files>')")
            print("Where <flags> are:")
            print("-A optionally scale intensities by count of spectra averaged")
            print("-B <sample> Set first sample to plot (default is 1/4 of samples)")
            print("-BASE  Fit and remove a spectral baseline")
            print("-C optionally flag the center of the band")
            print("-E <sample> Set last sample to plot (default is end of samples)")
            print("-H optionally set the high velocity region for baseline fit")
            print("-I <integration time> Time (seconds) to average observations before plotting")
            print("-K <dir> optionally keep average hot and cold load calibration observations")
            print("-L optionally set the low velocity region for baseline fit")
            print("-N <number> optionally set the number of spectra to plot")
            print("-O <dir> optionally output intermediate, average files")
            print("-P <dir> write PNG and PDF files instead of showing plot")
            print("-Q optionally plot intensity versus freQuency, instead of velocity")
            print("-R optionally flag known RFI Lines")
            print("-S <filename> optionally set summary file name")
            print("-U optionally update reference frequency for a different line")
            print("   ie -U 1612.231, 1665.402, 1667.349, 1720.530 or 1420.40575")
            print("-V optionally plot velocity")
            print("-Z <file tag> optionally add tag to PDF and PNG file names")
            print("-0 optionally plot zero intensity line(s)")
            print("-MINEL optionally set lowest elevation for calibration obs (default %7.1f)" % (self.lowel))            
            print("")
            print("Glen Langston - NSF - September 16, 2021")
            return ["",""]

        iarg = 0
        # must be some arguments, parse them
        while iarg < nargs:

            if self.verbose:
                print("Arg[%d] = %s" % (iarg, args[iarg]))
            
            # if folding data
            if args[iarg].upper() == '-F':
                print('Folding specectra')
                self.doFold = True
            elif args[iarg].upper() == '-A':
                self.toScalAve = True
            elif args[iarg].upper() == '-AVE':
                iarg = iarg + 1
                self.aveTimeSec = int( args[iarg])
                print('Averaging for %6.1f seconds' % (self.aveTimeSec))
            elif args[iarg].upper() == '-B':   # if setting beginning sample
                iarg = iarg + 1
                self.xa = int( args[iarg])
                print('Plotting starting at channel: %4d' % (self.xa))
            elif args[iarg].upper() == '-BASE':
                self.doBaseline = True
                print("Fitting and subtracking a baseline")
            elif args[iarg].upper() == '-C':
                self.flagCenter = True           
            elif args[iarg].upper() == '-F':
                iarg = iarg+1
                self.fitOrder = int( sys.argv[iarg])
                if self.fitOrder < 0:
                    self.fitOrder = 0
                elif fitOrder > 10:
                    self.fitOrder = 10
                doPoly = True
                if self.fitOrder == 0:
                    print("Fitting a constant baseline")
                elif self.fitOrder == 1:
                    print("Fitting a linear baseline")
                else:
                    print(("Fitting a %d-nd order polynomical baseline" % \
                          (self.fitOrder)))
            elif args[iarg].upper() == '-E':   # if setting ending sample
                iarg = iarg + 1
                self.xb = int( args[iarg])
            elif args[iarg].upper() == '-H':
                iarg = iarg+1
                self.maxvel = np.float( args[iarg])
                print('Maximum (high) velocity for sum: %7.2f km/sec'\
                      % (self.maxvel))
            elif args[iarg].upper() == "-I":
                iarg=iarg+1
                self.tint = np.float(args[iarg])
                print('Spectral integration time for averaging: %f' % (self.tint))
            elif args[iarg].upper() == '-K':
                self.doKeep = True
                iarg = iarg+1
                self.keepFileDir = str(args[iarg])
            elif args[iarg].upper() == '-L':
                iarg = iarg+1
                self.minvel = np.float( args[iarg])
                print('Minium (low)  velocity for sum: %7.2f km/sec' \
                      % (self.minvel))
            elif args[iarg].upper() == '-N':   # if number of spectra to plot
                iarg = iarg+1
                self.maxPlot = int(args[iarg])
                if self.maxPlot < 1:
                    print("Not Plotting")
                else:
                    print("Plot will have a maximum of %d spectra" \
                          % (self.maxPlot))
            elif args[iarg].upper() == '-O':
                self.doOutFile = True
                iarg = iarg+1
                self.outFileDir = str(args[iarg])
            elif args[iarg].upper() == '-P':
                self.doPlotFile = True
                iarg = iarg+1
                self.plotFileDir = str(args[iarg])
            elif args[iarg].upper() == '-Q':
                self.plotFrequency = True
            elif args[iarg].upper() == '-R':
                self.flagRfi = True
                print("Flagging RFI")
            elif args[iarg].upper() == '-T':   # if plot title provided
                iarg = iarg+1
                self.myTitle = str(args[iarg])
                print('Plot Title : ', self.myTitle)
            elif args[iarg].upper() == '-V':   # default is plotting Frequency
                self.plotFrequency = False              # plot velocity
            elif args[iarg].upper() == '-VA':   # now look for flags with arguments
                iarg = iarg+1
                self.minvel = float(args[iarg])
                print('Minimum velocity for baseline fit: %7.2f km/sec ' \
                      % (self.minvel))
            elif args[iarg].upper() == '-VB':   # now look for flags with arguments
                iarg = iarg+1
                self.maxvel = float(args[iarg])
                print('Maximum velocity for baseline fit: %7.2f km/sec '  % (self.maxvel))
            elif args[iarg].upper() == '-CEN':   # if nU ref is provided (in MHz)
                iarg = iarg+1
                self.nuRefFreq = float(args[iarg])
                print( 'Reference Frequency : %9.3f MHz' % (self.nuRefFreq))
            elif args[iarg].upper() == '-MINEL':  # if min elevation
                iarg = iarg + 1
                self.lowel = float( args[iarg])
                print(( "Using elevations > %7.2f (d) for Cold load obs." \
                        % (self.lowel)))
            elif args[iarg].upper() == '-W':
                self.writeKelvin = True
            elif args[iarg].upper() == '-Z':     # label written files
                iarg = iarg+1
                self.fileTag = str(args[iarg])
                print( 'File tag: %s' % (self.fileTag))
            elif args[iarg].upper() == '-0':
                doZero = True
                print('Plotting zero intensity lines')
            else:
                break
            iarg = iarg + 1
        # returns the file names/directories
        if iarg >= nargs:
            return ["",""]
        names = args[iarg:]
        if self.verbose:
            print("Observation Files and Directories:")
            print(names)
        return names
    # end of help/parsing arguments

    def Help(self, argstring=""):
        ###
        # Alias for help()
        ###
        names = self.help(argstring)
        return names

    def check_obs( self, ave_spec, in_spec):
        """
        check_obs() looks for changes in the observing parameters and/or
        integration time exceeded.
        This program keeps track of the last spectrum, that was not
        averaged because of the gain/angle/time change.  The last spectrum
        becomes "next"
        Inputs:
        ave_spec - average spectrum if self.nave > 0
        in_spec  - new input spectrum, by comparision a new observations is detected
        self.nave count > 1       indicates an average had been started.
        self.next.durationSec > 0 indicates the next spectrum must be processed.
        """
 
        newobs = False    # start assuming this is a continueed observation.
        if self.nave <= 0.:    # if no previous observatitions in avevarage
            # check if an observation change was detected in last round and
            # kept because it was different than the average.
            if self.next.durationSec > 0.:
                # if the previous file duration was longer than the integration time
                if self.next.durationSec > self.tint:
                    # just keep the input observation and return previous
                    ave_spec = self.next
                    self.next = in_spec
                    self.beginutc = self.next.utc
                    self.endutc =   self.next.utc
                    newobs = True
                    self.nave = 0
                    return newobs, ave_spec
                else: # observations are shorter than the averaging time
                    self.beginutc = self.next.utc
                    self.endutc = self.next.utc
                    ave_spec, self.nave, self.beginutc, self.endutc = \
                        self.average_spec( ave_spec, self.next, self.nave, \
                                           self.beginutc, self.endutc)
                    self.next.durationSec = 0. # don't sum the next spectrum again
            else:  # else check if input spectrum exceeds average time
                self.beginutc = in_spec.utc
                self.endutc = in_spec.utc
                #self.nave is updated by self.average_spec()
                if in_spec.durationSec > self.tint:
                    newobs = True
                    return newobs, in_spec
                # else no previous spectra, just save the weighted input spectrum
                ave_spec, self.nave, self.beginutc, self.endutc = \
                    self.average_spec( ave_spec, in_spec, self.nave, self.beginutc, self.endutc)
        # now check if input spectrum is the start of a new spectrum
        
        # check if an elevation change, which requires completing last and keeping new
        if in_spec.telel != ave_spec.telel or in_spec.telaz != ave_spec.telaz:
            newobs = True
            if self.verbose:
                if in_spec.telel != ave_spec.telel:
                    print("Elevation change from: %7.1f to %7.1f (deg)" % \
                      (ave_spec.telel, in_spec.telel))
                elif in_spec.telaz != ave_spec.telaz:
                    print("Azimuth change from: %7.1f to %7.1f (deg)" % \
                          (ave_spec.telaz, in_spec.telaz))
            self.next = in_spec
            ave_spec = self.normalize_spec( ave_spec, self.beginutc, self.endutc)
            self.nave = 0
            return newobs, ave_spec
        
        # if here, then not yet a new observation, integrate input into average                  
        ave_spec, self.nave, self.beginutc, self.endutc = \
            self.average_spec( ave_spec, in_spec, self.nave, self.beginutc, self.endutc) 
        
        # compute the duration from begin and end utcs.
        aveutc, duration = radioastronomy.aveutcs( self.beginutc, self.endutc)
        # if duration exceeded, then flag this is end of last observation
        if duration > self.tint:
            ave_spec = self.normalize_spec( ave_spec, self.beginutc, self.endutc)
            self.nave = 0
            newobs = True
                          
        # end of check_obs()
        return newobs, ave_spec
        
    def average_spec( self, ave_spec, in_spec, nave, firstutc, lastutc):
        """
        Averages spectra/weight functions by observing duration
        input/output
        ave_spec   spectrum structure containing weighted average
        input:
        in_spec    spectrum to be added to the average
        in/out:
        nave       Count of number of spectra averaged so far.
                   nave = 0 implies initializing the sum.
        firstutc   date of first observation averaged
        lastutc    date of current last date of spectrum to be averaged
        """

        if self.verbose:
            nData = len(in_spec.ydataA)
            n6 = int(nData/6)
            n56 = int(5*n6)            
            medianData = np.median( in_spec.ydataA[n6:n56])
            print(( "Input: %8.3f, count: %d" % (medianData, in_spec.count)))

        # if restarting the sum
        if nave == 0:
            ave_spec = copy.deepcopy(in_spec)  # initial spectrum is 1st
            firstutc = in_spec.utc
            lastutc = in_spec.utc
            self.nData = len(in_spec.ydataA)
            nave = 1
            # replace with duration scaled intensities
            ave_spec.ydataA = (in_spec.ydataA * in_spec.durationSec)/float(in_spec.nave)
            # keep track of observing time for weighted sum
            ave_spec.durationSec = in_spec.durationSec
        else: # else not enough time yet, average observations
            if in_spec.utc < firstutc:
                firstutc = in_spec.utc
            elif in_spec.utc > lastutc:
                lastutc = in_spec.utc
            nave = nave + 1
            ave_spec.ydataA = ave_spec.ydataA + \
                (in_spec.ydataA * in_spec.durationSec/float(in_spec.nave))
            # keep track of observing time for weighted sum
            ave_spec.durationSec = ave_spec.durationSec + in_spec.durationSec
            ave_spec.count = ave_spec.count + in_spec.count

        # end of average_spec()
        return ave_spec, nave, firstutc, lastutc

    def normalize_spec( self, ave_spec, firstutc, lastutc):
        """
        Normaize the average after completing sum
        input/output
        ave_spec   input raw sum of observation 
               output normalized sum of observations
        input      first and last utcs in observation
        Note:  The count of observations is not used, rather the integration
        is weighted by durations.  This corrects for observations with 
        different integration times.
        """
        # now renormalize for total integration time
        if ave_spec.durationSec > 0.:
            ave_spec.ydataA = ave_spec.ydataA/float(ave_spec.durationSec)
        else:
            print( "Normalizing spectrum with no integration time!")
        # compute average time from first and last utcs
        aveutc, duration = radioastronomy.aveutcs( firstutc, lastutc)
        ave_spec.utc = aveutc
        #need to re-calculate representative RA,Dec for average time
        ave_spec.azel2radec()

        # end of normalize_spec()
        return ave_spec

    def read_hot( self, names):
        """
        read_hot() reads in all files in the names list and averages hot load
        observations.   The hot load file is returned.
        While reading the files, the minmum elevation and 
        galactic latitudes are recorded
        """

        rs = radioastronomy.Spectrum()
        nhot = 0       # init count of hot files
        ncold = 0


        # if reading a hot file
        if (self.hotFileName != ""):
            self.ave_hot.read_spec_ast(self.hotFileName)
            self.ave_hot.ydataA = self.ave_hot.ydataA/self.ave_hot.nave
            self.nhot = 1
            # prepare transfer for gain calculation
            self.hv = self.ave_hot.ydataA
            return self.nhot
        # only process the hot files
        hotnames, nhot = rasnames.splitNames( names, ".hot", "")
        if nhot < 1:
            print("No Hot files: can not calibrate")
            return
        
        nhot = 0;  # now restart the average, in case some files are not OK
        self.nave = 0;  # make sure not files remain from last average
        self.next.durationSec = 0.
        
        # now run through and find hot and cold loads obs        
        for filename in hotnames:

            rs.read_spec_ast(filename)
          
            nChan = len( rs.ydataA)
            if nChan != 32 and nChan != 64 and nChan != 128 and nChan != 256 \
               and nChan != 512 and nChan != 1024 and nChan != 2048 and \
                   nChan != 4096:
                print("Unusual data length (%d), file %s" % (nChan, filename))
                continue
            
            # if a hot load observation
            if rs.telel < 0:
                if nhot == 0:

                    self.firstutc = rs.utc
                    self.lastutc = rs.utc
                    # nhot will be updated by average_spec()
                # accumulate spectra, nhot is updated also
                self.ave_hot, nhot, self.firstutc, self.lastutc = \
                    self.average_spec( self.ave_hot, rs, nhot, self.firstutc, self.lastutc)
        if nhot > 0:
            print(( "Found %3d Hot load observations" % (nhot)))
            self.ave_hot = self.normalize_spec( self.ave_hot, self.firstutc, self.lastutc)
        else:
            print( "No Hot load data, can not calibrate")
            return 0
            # exit()

        self.nhot = nhot
        # do more cleanup on spectra for RFI
        yv = copy.deepcopy(self.ave_hot.ydataA)
        if self.flagRfi:
            xv = self.ave_hot.xdata * 1.E-6
            # interpolate rfi
            hv = interpolate.lines( self.linelist, self.linewidth, xv, yv)
        else:
            hv = yv

        if self.flagCenter:             # if flagging spike in center of plot
    # remove spike in center of the plot
            icenter = int(self.nData/2)
            hv[icenter] = (hv[icenter-2] + hv[icenter+2])*.5
            hv[icenter-1] = (3.*hv[icenter-2] + hv[icenter+2])*.25
            hv[icenter+1] = (hv[icenter-2] + 3.*hv[icenter+2])*.25

        # all outputs are part of object
        self.ave_hot.ydataA = hv
        self.hv = hv
        
        # if keeping hot and cold files
        if self.doKeep and self.hotFileName == "":
            outname = radioastronomy.utcToName( self.ave_hot.utc)
            outname = outname + ".hot"  # output in counts
            # add telescope index
            outname = ("T%d-" % cpuIndex)+outname
            n = self.ave_hot.nChan
            print("Ave Hot: %d: %.6f" % \
                  (n, np.median( self.ave_hot.ydataA[int(n/3):int(2*n/3)])))
            # multiply by number of spectra averaged to yield a plot
            # with sufficient numerical resolution.
            # the nave factor is removed on read
            self.ave_hot.ydataA = self.ave_hot.ydataA*self.ave_hot.nave
            self.ave_hot.write_ascii_file(self.outputDir, outname)
            print( "Wrote Average Hot  Load File: %s%s" % ("../", outname))

        # end of read hot
        return nhot

    def read_angles( self, names):
        """
        read_angles() reads all files and counts number of files with 
        high elevation and high galactic latitude
        Inputs:
        lowel   minimum elevation to accept for cold load 
        """

        names, count = rasnames.splitNames(names, ".ast", ".hot")
        
        # count types
        ncold = 0
        rs = radioastronomy.Spectrum()
        nName = len(names)

        count = 0
        # now average coldest data for calibration
        for filename in names:

            rs.read_spec_ast(filename)
            rs.azel2radec()    # compute ra,dec from az,el
            if count == 0:
                minel = rs.telel
                maxel = rs.telel
                minGlat = rs.gallat
                maxGlat = rs.gallat
            count = count + 1    

            if rs.telel < self.lowel:  #if elevation too low for a cold load obs
                continue
            if rs.gallat > maxGlat:
                maxGlat = rs.gallat
            if rs.gallat < minGlat:
                minGlat = rs.gallat
            if minel > rs.telel:
                minel = rs.telel
            if maxel < rs.telel:
                maxel = rs.telel
        
            # count number of high elevation files
            if rs.telel > self.lowel:
                ncold = ncold + 1

        print(( "Found %d high elevation obs in %d files" % (ncold, nName)))
        self.ncold = ncold
        self.minel = minel
        self.maxel = maxel
        self.minGlat = minGlat
        self.maxGlat = maxGlat

        # end of read_angles()
        return ncold

    def read_cold( self, names, lowel, lowGlat):
        """
        read_cold() checks files and averages selected files with 
        high elevation and galactic Latitude.
        Inputs:
        lowel   minimum elevation to accept for cold load 
        lowGlat minimum galactic latitude
        """
        # starting new sum of cold (high elevation and galactic latitude) obs
        ncold = 0
        rs = radioastronomy.Spectrum()
        nName = len(names)

        nread = 0
        # now average coldest data for calibration
        for filename in names:

            rs.read_spec_ast(filename)
            if rs.telel < self.lowel:  #if elevation too low for a cold load obs
                continue
                
            # note this test excludes low galactic latitude ranges
            Glat = rs.gallat   # find min and max absolute value galactic Latitudes.
            if Glat < 0.:
                Glat = - Glat
            # if the Galactic latitude is above the minimum for calibration
            if nread == 0:
                self.minGlat = Glat
                self.maxGlat = Glat
            nread = nread + 1
            if Glat < self.minGlat:
                self.minGlat = Glat
            if Glat > self.maxGlat:
                self.maxGlat = Glat
            if Glat > lowGlat:
                # if first acceptable cold file, then init values
                if ncold == 0:
                    self.nave = 0
                    self.firstutc = rs.utc
                    self.lastutc = rs.utc
                    # initial values for min/max searches
                    self.minel = rs.telel
                    self.maxel = rs.telel
                    self.next.durationSec = 0.
                    # ncold is updated by aveerage_spec()
                self.ave_cold, ncold, self.firstutc, lastutc = \
                    self.average_spec( self.ave_cold, rs, ncold, self.firstutc, self.lastutc)
                if rs.telel < self.minel:
                    self.minel = rs.telel
                elif rs.telel > self.maxel:
                    self.maxel = rs.telel

            # end of all files to average
        if ncold < 1:
            print( "No high Galactic Latitude data")
            return 0
        else:
            self.ave_cold = self.normalize_spec( self.ave_cold, self.firstutc, self.lastutc)

        if self.flagRfi:
            cv = self.ave_cold.ydataA
            xv = self.ave_cold.xdata*1.E-6
            # interpolate rfi
            yv = interpolate.lines( self.linelist, self.linewidth, xv, cv)
        else:
            yv = self.ave_cold.ydataA
                
        if self.flagCenter:             # if flagging spike in center of plot
            nData = len( self.ave_cold.ydataA)
            # remove spike in center of the plot
            icenter = int(nData/2)
            yv[icenter] = (yv[icenter-2] + yv[icenter+2])*.5
            yv[icenter-1] = (3.*yv[icenter-2] + yv[icenter+2])*.25
            yv[icenter+1] = (yv[icenter-2] + 3.*yv[icenter+2])*.25
            
        self.ave_cold.ydataA = yv
        self.cv = yv
        self.ncold = ncold
            
        if self.verbose:
            print( "Found %3d High Galactic Latitude spectra" % (ncold))
            print("Min, Max Galactic Latitude: %7.1f,%7.1f" \
                   % (self.minGlat, self.maxGlat))
            print("Min, Max Elevation:         %7.1f,%7.1f" \
                   % (self.minel, self.maxel))
            
        # if keeping hot and cold files and did not read this file already
        if self.doKeep and self.coldFileName == "":
            outname = radioastronomy.utcToName( self.ave_cold.utc)
            outname = outname + ".cld"  # output in counts            
            # add telescope index
            outname = ("T%d-" % cpuIndex)+outname
            n = self.ave_cold.nChan
            print("Ave Cold: %d: %.6f" % \
                  (n, np.median( self.ave_hot.ydataA[int(n/3):int(2*n/3)])))
            # multiply by number of spectra averaged to yield a plot
            # with sufficient numerical resolution.
            # the nave factor is removed on read
            self.ave_cold.ydataA = self.ave_cold.ydataA*self.ave_cold.nave
            self.ave_cold.write_ascii_file(self.outputDir, outname)
            print( "Wrote Average Cold Load File: %s%s" % ("../", outname))

        # end of read_cold()
        return ncold

    def computeGain( self):
        """
        Compute the telescope gain (Kelvins per Count) based on hot and cold
        The hot and cold files must already be computed and intensities in self.hv and self.cv
        """

        nData = len( self.ave_hot.ydataA)
        n9 = int(nData/9)
        n39 = int(3*n9)
        n49 = int(4*n9)
        n59 = int(5*n9)
        n69 = int(6*n9)


        self.vel = np.zeros(nData)
        xv = self.ave_cold.xdata * 1.E-6 # convert to MHz
        # create index array
        for jjj in range (0, nData):
            self.vel[jjj] = self.c * (self.nuRefFreq - xv[jjj])/self.nuRefFreq

        self.xa, self.xb = gf.velocity_to_indicies( self.vel, \
                                                    self.minvel, self.maxvel)
        # compute gain on a channel by channel basis for tRx calculation
        gainHC = np.zeros(nData)
        for iii in range(nData):
            gainHC[iii] = (self.hv[iii] - self.cv[iii])/(self.thot - self.tcold)
            if gainHC[iii] < self.EPSILON:
                gainHC[iii] = self.EPSILON

        # the hot/cold gain ratios are only used to compute tRxMiddle
        trx = np.zeros(nData)
        for iii in range(nData):
            trx[iii] = (self.cv[iii]/gainHC[iii]) - self.tcold

        #Prepare to compute tRx, which is based only on cold load observations
        tRxA = np.median(trx[n39:n49])
        tRxB = np.median(trx[n59:n69])
        self.tRxMiddle = (tRxA + tRxB)*.5

        tStdA = np.std(trx[n39:n49])
        tStdB = np.std(trx[n59:n69])
        tRms  = (tStdA + tStdB) * .5

        #at this point, only the hot load observations are used to compute T  sys
        #No cold observations used for calibration, except for tRxMiddle

        print("Median Receiver Temp: %7.2f +/- %5.2f (%5.2f %5.2f) (K)" % \
              (self.tRxMiddle, tRms, tStdA, tStdB))

        # for remainder of calculations only use hot counts for calibration
        # Using hot load only reduces interference effects
        self.gain = np.zeros(nData)
        for iii in range(nData):
            self.gain[iii] = (self.thot + self.tRxMiddle)/self.hv[iii]
            if self.gain[iii] < self.EPSILON:
                self.gain[iii] = self.EPSILON

        gainA = np.median(self.gain[n39:n49])
        gainB = np.median(self.gain[n59:n69])
        self.gainAve = 2.0/(gainA + gainB)  # Report gain in K per Count 

        if self.verbose:
            print('Min Vel  %7.1f, Max Vel  %7.1f' % \
                  (self.minvel, self.maxvel))
            print('Min Chan %7d, Max Chan %7d' % (self.xa, self.xb))
            print('Min,Max Galactic Latitudes %7.1f,%7.1f (d)' % \
                  (self.minGlat, self.maxGlat))

        # end of computeGain
        return
 
    def writeTsys( self, in_spec):
        """
        write the fully calibrated Spectrum
        """
        tsky = in_spec.ydataA
        nData = in_spec.nChan
        # prepare indicies for statitics
        n6 = int(nData/6)
        n26 = 2*n6
        n46 = 4*n6
        n56 = 5*n6
            
        # get tsys from averages of ends of spectra
        tSys = np.median(tsky[n6:n56])
        tStdA = np.std(tsky[n6:n26])
        tStdB = np.std(tsky[n46:n56])
        cA = np.median(cv[n6:n26])
        cB = np.median(cv[n46:n56])
        counts = (cA+cB)/2.
        tStd = (tStdA+tStdB)/2.
        ave_spec.tSys = tSys
        ave_spec.tRx = self.tRxMiddle
        ave_spec.tRms = tStd
        ave_spec.tint = ave_spec.durationSec
        ave_spec.bunit = 'Kelvins'
        ave_spec.KperC = self.gainAve

    def raw(self, names=[""]):
        """
        Plot all files in the names list
        """
        # to create plots in cronjobs, must use a different backend
        if self.doPlotFile:
            mpl.use('Agg')
        import matplotlib.pyplot as plt

        if self.plotFrequency:
            print( "Ploting Intensity versus Frequency")
        else:
            print( "Ploting Intensity versus Velocity")

        linestyles = ['-','-','-', '--', '--','--','-.', '-.', '-.']
        colors = ['b','g','r','m','c']
        nstyles = len(linestyles)
        ncolors = len(colors)
        xallmax = -9.e9
        xallmin =  9.e9
        yallmax = -9.e9
        yallmin =  9.e9

        # initialize spectrum for reading and plotting
        rs = radioastronomy.Spectrum()

        # find all relevant files in directories and lists
        files, count = rasnames.splitNames(names, ".ast", ".hot", doDebug=self.verbose)
        if count < 1:
            print("No RAS files, can not calibrate")
            return
        
        count = 0 # restart the counting
        # look at all the files, but will stop after max plots
        for filename in files:
            # a file name must have a 3 letter extension (ie .hot, .ast)
            if len(filename) < 5:
                continue
            if self.verbose:
                print('%s' % (filename))

            rs.read_spec_ast( filename)
            if count == 0:
                ave_spec = rs
                self.nave = 0
            count = count + 1
            
            # now compute strings for plotting
            utcstr = rs.utc.isoformat()
            time, date, self.firstdate, self.lastdate = \
                rasnames.parsetime( utcstr, self.firstdate, self.lastdate)

            newobs, ave_spec = self.check_obs( ave_spec, rs)
            # wait until a new obs to plot
            if not newobs:
                continue
                
            xv = ave_spec.xdata  * 1.E-6 # convert to MHz
            # deal with spectra with different numbers of channels.
            nData = len( xv)
            n6 = int(nData/6)
            n56 = 5*n6

            if self.xa < 0:
                self.xa = 0
            if self.xb < 0:
                self.xb = nData
            if self.xb > nData:
                self.xb = nData

            yv = ave_spec.ydataA
            # ignore the 1/6 of the data at each end of spectrum
            ymedian = np.median(yv[n6:n56])

            if self.flagRfi:   # RFI is always flagged in topocentric frequencies
                # interpolate rfi
                yv = interpolate.lines( self.linelist, self.linewidth, xv, yv) 

            if not self.plotFrequency:   # if plotting velocity
                self.vel = np.zeros(nData)
                for jjj in range (0, nData):
                    self.vel[jjj] = self.c * (self.nuRefFreq - xv[jjj])/self.nuRefFreq
                self.xa, self.xb = gf.velocity_to_indicies( self.vel, \
                                                            self.minvel, self.maxvel)
                xv = self.vel
            # keep track of x axis for plot sizing
            xmin = min(xv[self.xa:self.xb])
            xmax = max(xv[self.xa:self.xb])
            xallmin = min(xmin,xallmin)
            xallmax = max(xmax,xallmax)

            if self.flagCenter:             # if flagging spike in center of plot
                # remove spike in center of the plot
                icenter = int(nData/2)
                yv[icenter] = (yv[icenter-2] + yv[icenter+2])*.5
                yv[icenter-1] = (3.*yv[icenter-2] + yv[icenter+2])*.25
                yv[icenter+1] = (yv[icenter-2] + 3.*yv[icenter+2])*.25

            # prepare to report intensity summary
            ymin = min(yv)
            ymax = max(yv)
            ymed = np.median(yv)
            coordtitle = '  Time   AZ,EL (deg)  Lon,Lat (deg)' 
            if self.nplot == 0:
                print("%s    Max   Median    Count  " % (coordtitle))
            coordlabel = '%s %5s,%5s  %5.1f,%5.1f' % \
                (time, ave_spec.telaz, ave_spec.telel, ave_spec.gallon, ave_spec.gallat)        
            # summarize the minimum and median values
            print('%s %8.3f %8.3f  %8d' % (coordlabel, ymax, ymed, ave_spec.count))
            if self.nplot <= 0 and self.maxPlot > 0:
                fig,ax1 = plt.subplots(figsize=(10,6))
                # fig.canvas.set_window_title(date)
                for tick in ax1.xaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 
                for tick in ax1.yaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 

            # if maximum number of plots completed, exit to show the results
            if self.nplot > self.maxPlot:
                break
            note = ave_spec.site
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            plt.xlim(xallmin,xallmax)

            if self.plotFrequency:
                plt.plot(xv[self.xa:self.xb], yv[self.xa:self.xb], \
                         colors[(self.nplot % ncolors) ], \
                         linestyle=linestyles[(self.nplot % nstyles)], label=coordlabel, lw=2)
            else:
                plt.plot(xv[self.xa:self.xb], yv[self.xa:self.xb], \
                         colors[(self.nplot % ncolors) ], \
                         linestyle=linestyles[(self.nplot % nstyles)], label=coordlabel, lw=2)
            self.nplot = self.nplot + 1
        # end for all names
        # clean up averaging parameters
        self.nave = 0
        self.next.durationSec = 0.
        
        if (self.maxPlot < 1) or (self.nplot < 1):
            print("No Plots, exiting")
            return
            # exit()
    
        if self.myTitle == "":
            self.myTitle = note
  
        # scale min and max intensities for nice plotting
        plt.ylim(0.9*yallmin,1.25*yallmax)

        plt.title(self.myTitle, fontsize=16)
        plt.xlabel('Frequency (MHz)',fontsize=16)
        ylabel = 'Intensity (%s)' % rs.bunit
        plt.ylabel(ylabel, fontsize=16)
        plt.legend(loc='upper right')
        # zero the plot count for next execution
        self.nplot = 0 
 
        # if writing files
        if self.doPlotFile:
            self.firstdate = self.firstdate[2:]
            if self.fileTag == "":
                self.fileTag = "R-" + self.firstdate
            outpng = self.outFileDir + self.fileTag + ".png"
            plt.savefig(outpng,bbox_inches='tight')
            outpdf = self.outFileDir + self.fileTag + ".pdf"
            plt.savefig(outpdf,bbox_inches='tight')
            print( "Wrote files %s and %s" % (outpng, outpdf))
        else:
            # else show the plots
            plt.show()

        #end of ras.raw(names)
        return

    def Tintegrate( self, in_spec, iVmin, iVmax):
        """
        Tintegrate() computes two integrals of spectrum
        1) The kelvins-km/sec value for the velocity range
        2) The intensity weighted average velocity
        The program assumes the spectral velocity is pre-computed with corrections
        The indicies iVmin and iVmax are regions of the spectra containing the sources
        """
        # self.xa, self.xb are indices to the velocity range to plot.  pre-computed
        icenter = int(in_spec.nChan/2)
        dV = (self.vel[icenter+4] - self.vel[icenter - 4])*.125
        if dV < 0:
            dV = - dV

        iVmin = int(iVmin)
        iVmax = int(iVmax)
        if iVmin < 0:
            iVmin = 0
        if iVmax > in_spec.nChan:
            iVmax = in_spec.nChan
        nv = iVmax - iVmin
        
        # create sub-arrays of intensity and velocity
        tSs = in_spec.ydataA[iVmin:iVmax]
        # the spectra only have a frequency axis.  the Velocity axis is separate.
        vSs = self.vel[iVmin:iVmax]
        tSourcemin = min(tSs)
        tSourcemax = max(tSs)
        # get index to maximum value; then get velocity
        iSourcemax = np.argmax(tSs)
        velSource = vSs[iSourcemax]
        # integrate over spectrum for required velocity range
        tSum = np.sum(tSs)
        tSumRms = np.std(tSs)

        # computed the integrated velocity
        tvs = tSs*vSs
        tVSum = np.sum(tvs)
        tVSumRms = np.std(tvs)

        # Integration is reported in Kelvin*Km/Sec;
        # Multiply by source velocity range
        tSumKmSec = tSum * ( self.maxSvel - self.minSvel)/float(nv)
        dTSumKmSec = tSumRms * dV * np.sqrt(float(nv))
        if self.verbose:
            print("T-V Sum :  %10.3f +/- %7.3f" % (tSumKmSec, dTSumKmSec))
            print("Velocity:  %10.3f +/- %7.3f" % (tSVsum, tVSumRms))
            
        # end of Tintegrate()
        return tSumKmSec, dTSumKmSec, tVSum, tVSumRms
               
    def tsys(self, names):
        ###
        # Plot Tsys calibrated spectra
        ###
        astnames, nast = rasnames.splitNames(names,".ast", "", doDebug=self.verbose)
        hotnames, nhot = rasnames.splitNames(names,".hot", "", doDebug=self.verbose)

        print("Reading Hot Load observations")
        nhot = self.read_hot( hotnames)
        if nhot < 1:
            print("No Hot load Observations, can not calibrate")
 
        print("Reading Cold Sky observations")
        ncold = self.read_cold( astnames, self.lowel, self.lowGlat)
        if ncold < 1:
            self.lowel = 1.
            self.lowGlat = self.maxGlat - 5.
            ncold = self.read_cold( astnames, self.lowel, self.lowGlat)
            print("Min, Max El: %6.1f,%6.1f deg; Min, Max Glat: %6.2f, %6.2f deg" % \
                  (self.minel, self.maxel, self.minGlat, self.maxGlat))

            if ncold < 1:
                print("No Cold load files above minimum eleation %8.1f deg" % (self.lowel))
                return
            else:
                print("Found %d Cold Sky Observations" % (ncold))
            
        # now, based on hot and cold loads, compute gains for remaining spectra
        self.computeGain()

        # to create plots in cronjobs, must use a different backend
        if self.doPlotFile:
            mpl.use('Agg')
        import matplotlib.pyplot as plt

        if self.plotFrequency:
            print( "Ploting Intensity versus Frequency")
        else:
            print( "Ploting Intensity versus Velocity")

        linestyles = ['-','-','-', '--', '--','--','-.', '-.', '-.']
        colors = ['b','g','r','m','c']
        nstyles = len(linestyles)
        ncolors = len(colors)
        xallmax = -9.e9
        xallmin =  9.e9
        yallmax = -9.e9
        yallmin =  9.e9

        # initialize spectrum for reading and plotting
        rs = radioastronomy.Spectrum()
        
        count = 0
        # plot no more than N spectra
        for filename in astnames:
            # a file name must have a 3 letter extension (ie .hot, .ast)
            if len(filename) < 5:
                continue
            if self.verbose:
                print('%s' % (filename))

            rs.read_spec_ast( filename)
            # need to initialize the ave spectrum similar to read spectrum
            if count == 0:
                ave_spec = rs
                self.nave = 0
            count = count + 1
            
            # now compute strings for plotting
            utcstr = rs.utc.isoformat()
            # provide an ascii string of UTC time and update dates
            time, date, self.firstdate, self.lastdate = \
                rasnames.parsetime( utcstr, self.firstdate, self.lastdate)

            newobs, ave_spec = self.check_obs( ave_spec, rs)
            # wait until a new obs to plot
            if not newobs:
                continue
                
            xv = ave_spec.xdata  * 1.E-6 # convert to MHz
            # deal with spectra with different numbers of channels.
            nData = len( xv)
            n6 = int(nData/6)
            n56 = 5*n6

            if self.xa < 0:
                self.xa = 0
            if self.xb < 0:
                self.xb = nData
            if self.xb > nData:
                self.xb = nData

            yv = ave_spec.ydataA
            # ignore the 1/6 of the data at each end of spectrum
            ymedian = np.median(yv[n6:n56])

            if self.flagRfi:   # RFI is always flagged in topocentric frequencies
                # interpolate rfi
                yv = interpolate.lines( self.linelist, self.linewidth, xv, yv) 

            if self.flagCenter:             # if flagging spike in center of plot
                # remove spike in center of the plot
                icenter = int(nData/2)
                yv[icenter] = (yv[icenter-2] + yv[icenter+2])*.5
                yv[icenter-1] = (3.*yv[icenter-2] + yv[icenter+2])*.25
                yv[icenter+1] = (yv[icenter-2] + 3.*yv[icenter+2])*.25                
            # finally correct for the system gain
            tsky = yv * self.gain
            
            if not self.plotFrequency:   # if plotting velocity
                self.vel = np.zeros(nData)
                for jjj in range (0, nData):
                    self.vel[jjj] = self.c * (self.nuRefFreq - xv[jjj])/self.nuRefFreq

                # compute velocity correction for this direction and date
                corr = gf.compute_vbarycenter( ave_spec)
                if self.verbose:
                    print("Barycentric Correction: %9.3f km/sec" % (corr))
                # make the topcentric velocity correction
                self.vel = self.vel + corr
                          
                self.xa, self.xb = gf.velocity_to_indicies( self.vel, \
                                                                self.minvel, self.maxvel)
             
                # Computes and subtracts baseline for source intensities
                # the indicies self.xa, xb are centers for two regions for fit.
                # The fit excludes the central velocities.
                nchanfit = 20
                baseline = gf.fit_baseline( self.vel[0:nData], tsky[0:nData], \
                                           self.xa, self.xb, nchanfit, self.fitOrder)
            
                # indicies to range of channels containing source.
                iVmin, iVmax = gf.velocity_to_indicies( self.vel, self.minSvel, self.maxSvel)

                # remove baseline to get the source spectrum
                tSource = tsky[0:nData] - baseline[0:nData]

                # if plotting/keeping the baseline subtracted spectra, transfer to Sky
                if self.doBaseline:
                    tsky = tSource
               
                # now compute integrated intensity and noise estimates
                nv = iVmax - iVmin
                if nv < 0:           
                    print(("Velocity Index Error: %d > %d" % (iVmin, iVMax)))
                    nv = -nv
                    
                # transfer for saving/calculating
                ave_spec.ydataA = tSource
                # compute/log integrals of this spectrum
                tSumKmSec, dTSumKMSec, tVsum, tVSumRm = \
                    self.Tintegrate( ave_spec, iVmin, iVmax)
                
                xv = self.vel[self.xa:self.xb]
                yv = tsky[self.xa:self.xb]
            else:
                yv = tsky
                

            # keep track of x axis for plot sizing
            xmin = min(xv)
            xmax = max(xv)                          
            xallmin = min(xmin,xallmin)
            xallmax = max(xmax,xallmax)

            # prepare to report intensity summary
            ymin = min(yv)
            ymax = max(yv)
            ymed = np.median(yv)
            coordtitle = '  Time   AZ,EL (deg)  Lon,Lat (deg)' 
            if self.nplot == 0:
                print("%s  Max-(K)-Median     Count  " % (coordtitle))
            coordlabel = '%s %5s,%5s  %5.1f,%5.1f' % \
                (time, ave_spec.telaz, ave_spec.telel, ave_spec.gallon, ave_spec.gallat)        
            # summarize the minimum and median values
            print('%s %8.2f %8.2f  %8d' % (coordlabel, ymax, ymed, ave_spec.count))
            if self.nplot <= 0 and self.maxPlot > 0:
                fig,ax1 = plt.subplots(figsize=(10,6))
                # fig.canvas.set_window_title(date)
                for tick in ax1.xaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 
                for tick in ax1.yaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 

            if self.nplot > self.maxPlot:
                continue
            note = ave_spec.site
            yallmin = min(ymin,yallmin)
            yallmax = max(ymax,yallmax)
            plt.xlim(xallmin,xallmax)
            ave_spec.bunit = "Kelvins"
            
            plt.plot(xv, yv, \
                         colors[(self.nplot % ncolors) ], \
                         linestyle=linestyles[(self.nplot % nstyles)], label=coordlabel, lw=2)
            self.nplot = self.nplot + 1
            
        # end for all names
        # clean up averaging parameters
        self.nave = 0
        self.next.durationSec = 0.
        
        if (self.maxPlot < 1) or (self.nplot < 1):
            print("No Plots, exiting")
            return
            # exit()
    
        if self.myTitle == "":
            self.myTitle = note
  
        # scale min and max intensities for nice plotting
        dy = yallmax - yallmin
        plt.ylim(yallmin-(dy*.1),yallmax+(dy*.1))

        plt.title(self.myTitle, fontsize=16)
        if self.plotFrequency:
            plt.xlabel('Frequency (MHz)',fontsize=16)
        else:
            plt.xlabel('Velocity (km/sec,  Corr=%9.2f)' % (corr),fontsize=16)
                
        ylabel = 'Intensity (%s)' % ave_spec.bunit
        plt.ylabel(ylabel, fontsize=16)
        plt.legend(loc='upper right')
        # zero the plot count for next execution
        self.nplot = 0 

        # if writing files
        if self.doPlotFile:
            self.firstdate = self.firstdate[2:]
            if self.fileTag == "":
                self.fileTag = "T-" + self.firstdate
            outpng = self.outFileDir + self.fileTag + ".png"
            plt.savefig(outpng,bbox_inches='tight')
            outpdf = self.outFileDir + self.fileTag + ".pdf"
            plt.savefig(outpdf,bbox_inches='tight')
            print( "Wrote files %s and %s" % (outpng, outpdf))
        else:
            # else show the plots
            plt.show()
      
        self.plotFrequency = True  
        # End of ras.tsys()
        return
    
