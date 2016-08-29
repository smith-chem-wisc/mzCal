using mzCal;
using System;
using System.Reflection;

namespace mzCalCommandLine
{
    class Program
    {
        static void Main(string[] args)
        {
            var version = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            if (version.Equals("1.0.0.0"))
                Console.WriteLine("Not a release version");
            else
                Console.WriteLine(version);

            string origDataFile = args[0];
            string mzidFile = args[1];

            mzCalIO.mzCalIO.Load();

            SoftwareLockMassParams a = mzCalIO.mzCalIO.GetReady(origDataFile, P_outputHandler, P_progressHandler, P_watchHandler, mzidFile);

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