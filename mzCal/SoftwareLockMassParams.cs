using MassSpectrometry;
using Spectra;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace mzCal
{
    public class SoftwareLockMassParams
    {
        #region isotopologue parameters
        // THIS PARAMETER IS FRAGILE!!!
        // TUNED TO CORRESPOND TO SPECTROMETER OUTPUT
        // BETTER SPECTROMETERS WOULD HAVE BETTER (LOWER) RESOLUIONS
        // Parameter for isotopolouge distribution searching
        public double fineResolution = 0.1;

        #endregion

        public event EventHandler<OutputHandlerEventArgs> outputHandler;
        public event EventHandler<ProgressHandlerEventArgs> progressHandler;
        public event EventHandler<OutputHandlerEventArgs> watchHandler;

        public HashSet<int> MS2spectraToWatch;
        public HashSet<int> MS1spectraToWatch;
        public DoubleRange mzRange;

        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
        public Identifications identifications;

        public delegate void PostProcessing(SoftwareLockMassParams p);
        public PostProcessing postProcessing;

        public delegate string GetFormulaFromDictionary(string dictionary, string acession);
        public GetFormulaFromDictionary getFormulaFromDictionary;
        public bool calibrateSpectra = true;
        internal int randomSeed;
        internal bool deconvolute;
        public string paramString = "";
        public int minMS2 = 2;
        public int minMS1 = 3;

        #region Constructors

        public SoftwareLockMassParams(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int randomSeed, bool deconvolute)
        {
            this.myMsDataFile = myMsDataFile;
            MS1spectraToWatch = new HashSet<int>();
            MS2spectraToWatch = new HashSet<int>();
            this.randomSeed = randomSeed;
            this.deconvolute = deconvolute;

        }

        #endregion

        public virtual void OnOutput(OutputHandlerEventArgs e)
        {
            var handler = this.outputHandler;
            if (handler != null)
                handler(this, e);
        }

        public virtual void OnProgress(ProgressHandlerEventArgs e)
        {
            var handler = this.progressHandler;
            if (handler != null)
                handler(this, e);
        }

        public virtual void OnWatch(OutputHandlerEventArgs e)
        {
            var handler = this.watchHandler;
            if (handler != null)
                handler(this, e);
        }
    }

    public class OutputHandlerEventArgs : EventArgs
    {
        public string output { get; private set; }
        public OutputHandlerEventArgs(string output)
        {
            this.output = output;
        }
        public OutputHandlerEventArgs()
        {
            output = "";
        }
    }

    public class ProgressHandlerEventArgs : EventArgs
    {
        public int progress { get; private set; }
        public ProgressHandlerEventArgs(int progress)
        {
            Debug.Assert(progress >= 0 && progress <= 100);
            this.progress = progress;
        }
    }
}