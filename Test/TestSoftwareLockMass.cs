using mzCal;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public sealed class TestSoftwareLockMass
    {

        [OneTimeSetUp]
        public void setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
            EventHandler<OutputHandlerEventArgs> p_outputHandler = P_outputHandler;
            outputHandler += p_outputHandler;
            EventHandler<ProgressHandlerEventArgs> p_progressHandler = P_progressHandler;

        }

        [Test]
        public void TestSoftwareLockMassRunner()
        {
            mzCalIO.mzCalIO.Load();

            SoftwareLockMassParams a = mzCalIO.mzCalIO.GetReady(@"myFakeFile.mzML", P_outputHandler, P_progressHandler, P_outputHandler, @"myIdentifications.mzid", true);

            SoftwareLockMassRunner.Run(a);
        }

        event EventHandler<OutputHandlerEventArgs> outputHandler;

        public void OnOutput(OutputHandlerEventArgs e)
        {
            var handler = this.outputHandler;
            if (handler != null)
                handler(this, e);
        }

        private static void P_outputHandler(object sender, OutputHandlerEventArgs e)
        {
            Console.WriteLine(e.output);
        }
        private static void P_progressHandler(object sender, ProgressHandlerEventArgs e)
        {
            Console.Write(e.progress + "% ");
        }



        [Test]
        public void TestConstantCalibration()
        {
            List<LabeledDataPoint> trainingList = new List<LabeledDataPoint>();
            trainingList.Add(new LabeledDataPoint(new double[3] { 1, 1, 1 }, 0.5));
            trainingList.Add(new LabeledDataPoint(new double[3] { 1, 2, 1 }, 0.5));
            trainingList.Add(new LabeledDataPoint(new double[3] { 2, 1, 1 }, 0.5));
            trainingList.Add(new LabeledDataPoint(new double[3] { 2, 2, 1 }, 0.5));
            IdentityCalibrationFunction cf = new IdentityCalibrationFunction(OnOutput);
            Assert.AreEqual(4 * Math.Pow(0.5, 2) / 4, cf.getMSE(trainingList));
            ConstantCalibrationFunction cfconst = new ConstantCalibrationFunction(OnOutput, trainingList);
            Assert.AreEqual(0, cfconst.getMSE(trainingList));
        }
    }
}