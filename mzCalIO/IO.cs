using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using mzCal;
using MzIdentML;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace mzCalIO
{
    public static class mzCalIO
    {
        public static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized;
        public static UsefulProteomicsDatabases.Generated.obo psimodDeserialized;
        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

        public static string unimodLocation = @"unimod_tables.xml";
        public static string psimodLocation = @"PSI-MOD.obo.xml";
        public static string elementsLocation = @"elements.dat";
        public static string uniprotLocation = @"ptmlist.txt";

        private static int GetLastNumberFromString(string s)
        {
            return Convert.ToInt32(Regex.Match(s, @"\d+$").Value);
        }

        public static SoftwareLockMassParams GetReady(string origDataFile, EventHandler<OutputHandlerEventArgs> p_outputHandler, EventHandler<ProgressHandlerEventArgs> p_progressHandler, EventHandler<OutputHandlerEventArgs> p_watchHandler, string mzidFile, bool deconvolute)
        {
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
            if (Path.GetExtension(origDataFile).Equals(".mzML"))
                myMsDataFile = new Mzml(origDataFile);
            else
                myMsDataFile = new ThermoRawFile(origDataFile);
            int randomSeed = 1;
            var a = new SoftwareLockMassParams(myMsDataFile, randomSeed, deconvolute);
            a.outputHandler += p_outputHandler;
            a.progressHandler += p_progressHandler;
            a.watchHandler += p_watchHandler;
            a.postProcessing = MzmlOutput;
            a.getFormulaFromDictionary = getFormulaFromDictionary;
            a.identifications = new MzidIdentifications(mzidFile);

            //a.MS2spectraToWatch.Add(13208);
            //a.MS2spectraToWatch.Add(15143
            //a.MS2spectraToWatch.Add(14296);
            //a.MS2spectraToWatch.Add(7228);
            //a.MS2spectraToWatch.Add(15221);
            //a.MS2spectraToWatch.Add(8192);

            //a.MS1spectraToWatch.Add(17639);
            //a.MS2spectraToWatch.Add(17706);



            //a.MS2spectraToWatch.Add(5998);
            //a.MS2spectraToWatch.Add(5987);
            //a.MS2spectraToWatch.Add(15727);
            //a.MS2spectraToWatch.Add(15511);
            //a.MS2spectraToWatch.Add(4869);
            //a.MS2spectraToWatch.Add(8491);
            //a.MS2spectraToWatch.Add(11502);
            //a.MS2spectraToWatch.Add(14795);
            //a.MS2spectraToWatch.Add(8928);
            //a.MS2spectraToWatch.Add(15737);
            //a.MS2spectraToWatch.Add(15187);
            //a.MS2spectraToWatch.Add(18092);
            //a.MS2spectraToWatch.Add(5326);
            //a.MS2spectraToWatch.Add(4648);
            //a.MS2spectraToWatch.Add(4636);
            //a.MS2spectraToWatch.Add(13465);
            //a.MS2spectraToWatch.Add(13449);
            //a.MS2spectraToWatch.Add(13353);
            //a.MS2spectraToWatch.Add(13343);
            //a.MS2spectraToWatch.Add(10457);
            //a.MS2spectraToWatch.Add(11823);
            //a.MS2spectraToWatch.Add(16893);
            //a.MS2spectraToWatch.Add(17454);
            //a.MS2spectraToWatch.Add(17486);
            //a.MS2spectraToWatch.Add(17258);
            //a.MS2spectraToWatch.Add(11791);
            //a.MS2spectraToWatch.Add(15694);
            //a.MS2spectraToWatch.Add(11018);
            //a.MS2spectraToWatch.Add(6682);
            //a.MS2spectraToWatch.Add(14091);
            //a.MS2spectraToWatch.Add(7676);
            //a.MS2spectraToWatch.Add(14736);
            //a.MS2spectraToWatch.Add(17784);
            //a.MS2spectraToWatch.Add(11010);



            a.mzRange = new DoubleRange(0, 0);
            return a;
        }

        public static void Load()
        {

            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
            unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
            psimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(psimodLocation);
            uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);
        }

        public static string getFormulaFromDictionary(string dictionary, string acession)
        {
            if (dictionary == "UNIMOD")
            {
                string unimodAcession = acession;
                var indexToLookFor = GetLastNumberFromString(unimodAcession) - 1;
                while (unimodDeserialized.modifications[indexToLookFor].record_id != GetLastNumberFromString(unimodAcession))
                    indexToLookFor--;
                return Regex.Replace(unimodDeserialized.modifications[indexToLookFor].composition, @"[\s()]", ""); ;
            }
            else if (dictionary == "PSI-MOD")
            {
                string psimodAcession = acession;
                UsefulProteomicsDatabases.Generated.oboTerm ksadklfj = (UsefulProteomicsDatabases.Generated.oboTerm)psimodDeserialized.Items[GetLastNumberFromString(psimodAcession) + 2];

                if (GetLastNumberFromString(psimodAcession) != GetLastNumberFromString(ksadklfj.id))
                    throw new Exception("Error in reading psi-mod file, acession mismatch!");
                else
                {
                    foreach (var a in ksadklfj.xref_analog)
                    {
                        if (a.dbname == "DiffFormula")
                        {
                            return Regex.Replace(a.name, @"[\s()]", "");
                        }
                    }
                    return uniprotDeseralized[GetLastNumberFromString(psimodAcession)].ThisChemicalFormula.Formula;
                }
            }
            else
                throw new Exception("Not familiar with modification dictionary " + dictionary);
        }

        public static void MzmlOutput(SoftwareLockMassParams p)
        {
            p.OnOutput(new OutputHandlerEventArgs("Creating _indexedmzMLConnection, and putting data in it"));
            var path = Path.Combine(Path.GetDirectoryName(p.myMsDataFile.FilePath), Path.GetFileNameWithoutExtension(p.myMsDataFile.FilePath) + "-Calibrated.mzML");
            MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(p.myMsDataFile, path);
        }

    }
}
