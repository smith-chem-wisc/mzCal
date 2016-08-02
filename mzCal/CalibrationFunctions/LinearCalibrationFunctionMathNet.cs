using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    public class LinearCalibrationFunctionMathNet : CalibrationFunction
    {
        Func<double[], double> f;
        private Action<OutputHandlerEventArgs> onOutput;
        private bool[] useFeature;
        private bool[] logVars;
        private int numFeatures;
        
        public LinearCalibrationFunctionMathNet(Action<OutputHandlerEventArgs> onOutput, List<LabeledDataPoint> trainList2, bool[] varsBool, bool[] logVars)
        {
            this.onOutput = onOutput;
            this.logVars = logVars;
            useFeature = varsBool;
            numFeatures = useFeature.Where(b => b == true).Count();
            Train(trainList2);
        }

        public override double Predict(double[] t)
        {
            return f(IndexMap(t));
        }


        private double[] IndexMap(double[] input)
        {
            double[] output = new double[numFeatures];
            int featInd = 0;
            for (int k = 0; k < useFeature.Length; k++)
                if (useFeature[k])
                {
                    if(logVars[k])
                        output[featInd] = Math.Log(input[k]);
                    else
                        output[featInd] = input[k];
                    featInd++;
                }
            return output;
        }

        public void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            double[][] ok = new double[trainingList.Count()][];
            int k = 0;
            foreach (LabeledDataPoint p in trainingList)
            {
                ok[k] = IndexMap(p.inputs);
                k++;
            }
            var ok2 = trainingList.Select(b => b.output).ToArray();
            
            var ye = new Func<double[], double>[numFeatures+1];
            ye[0] = a => 1;
            for (int i = 0; i< numFeatures; i++) {
                int j = i;
                ye[j + 1] = a => a[j];
            }
            f = Fit.LinearMultiDimFunc(ok, ok2, ye);
            onOutput(new OutputHandlerEventArgs("Finished fitting a linear"));
        }
    }
}