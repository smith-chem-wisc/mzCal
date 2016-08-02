using Chemistry;
using System;
using System.IO;

namespace mzCal
{
    static class Calibrators
    {
        public static void CalibrateTSV(CalibrationFunction cf, SoftwareLockMassParams p)
        {
            using (StreamReader reader = new StreamReader(p.tsvFile))
            using (StreamWriter file = new StreamWriter(p.tsvFile + ".calibrated"))
            {
                string line;
                int lines = 0;
                while ((line = reader.ReadLine()) != null)
                {
                    if (lines > 0)
                    {
                        // Do something with the line.
                        string[] parts = line.Split('\t');

                        // Needed for getting RT for precursor
                        int MS2spectrumNumber = Convert.ToInt32(parts[1]);
                        string precursorID;
                        p.myMsDataFile.GetScan(MS2spectrumNumber).TryGetPrecursorID(out precursorID);
                        int MS1spectrumNumber = MS2spectrumNumber - 1;
                        while (!precursorID.Equals(p.myMsDataFile.GetScan(MS1spectrumNumber).id))
                            MS1spectrumNumber--;

                        double precursorRT = p.myMsDataFile.GetScan(MS1spectrumNumber).RetentionTime;

                        // To Calibrate
                        double precursorMZ = Convert.ToDouble(parts[5]);

                        // To Recompute
                        double precursorMass = Convert.ToDouble(parts[8]);
                        double precursorMassErrorDA = Convert.ToDouble(parts[18]);
                        double precursorMassErrorppm = Convert.ToDouble(parts[19]);

                        // Other needed info
                        int precursorCharge = Convert.ToInt32(parts[7]);
                        double theoreticalMass = Convert.ToDouble(parts[17]);
                        double precursorIntensity = Convert.ToDouble(parts[6]);
                        double TotalIonCurrent = p.myMsDataFile.GetScan(MS1spectrumNumber).TotalIonCurrent;
                        double InjectionTime = p.myMsDataFile.GetScan(MS1spectrumNumber).InjectionTime;

                        precursorMZ -= cf.Predict(new double[6] { 1, precursorMZ, precursorRT, precursorIntensity, TotalIonCurrent, InjectionTime });
                        precursorMass = precursorMZ.ToMass(precursorCharge);
                        precursorMassErrorDA = precursorMass - theoreticalMass;
                        precursorMassErrorppm = precursorMassErrorDA / theoreticalMass * 1e6;

                        parts[5] = precursorMZ.ToString();
                        parts[8] = precursorMass.ToString();
                        parts[18] = precursorMassErrorDA.ToString();
                        parts[19] = precursorMassErrorppm.ToString();

                        file.WriteLine(string.Join("\t", parts));
                    }
                    else
                        file.WriteLine(line);
                    lines++;
                }
            }
        }

    }
}
