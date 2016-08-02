using System;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    public class ConstantCalibrationFunction : CalibrationFunction
    {
        private double a;
        private Action<OutputHandlerEventArgs> onOutput;

        public ConstantCalibrationFunction(Action<OutputHandlerEventArgs> onOutput, IEnumerable<LabeledDataPoint> trainingList)
        {
            this.onOutput = onOutput;
            Train(trainingList);
        }

        public override double Predict(double[] inputs)
        {
            return a;
        }

        public void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            a = trainingList.Select(b => b.output).Average();
            onOutput(new OutputHandlerEventArgs("a = " + a));
        }
    }
}