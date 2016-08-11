using mzCal;
using System;
using System.Reflection;

namespace mzCalCommandLine
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Version " + Assembly.GetExecutingAssembly().GetName().Version);

            string origDataFile = args[0];
            string mzidFile = args[1];
            bool deconvolute = false;
            if (args.Length > 2)
                deconvolute = args[2].Equals("deconvolute");

            mzCalIO.mzCalIO.Load();

            SoftwareLockMassParams a = mzCalIO.mzCalIO.GetReady(origDataFile, P_outputHandler, P_progressHandler, P_watchHandler, mzidFile, deconvolute);

            SoftwareLockMassRunner.Run(a);
            

            //Console.Read();
        }

        private static void P_progressHandler(object sender, ProgressHandlerEventArgs e)
        {
            Console.Write(e.progress + "% ");
        }

        private static void P_outputHandler(object sender, OutputHandlerEventArgs e)
        {
            Console.WriteLine(e.output);
        }

        private static void P_watchHandler(object sender, OutputHandlerEventArgs e)
        {
            Console.WriteLine(e.output);
        }
    }
}