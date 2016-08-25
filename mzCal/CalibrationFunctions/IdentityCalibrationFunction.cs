using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    public class IdentityCalibrationFunction : CalibrationFunction
    {
        private Action<OutputHandlerEventArgs> onOutput;

        public IdentityCalibrationFunction(Action<OutputHandlerEventArgs> onOutput)
        {
            this.onOutput = onOutput;
        }

        public override double Predict(double[] inputs)
        {
            return 0;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            onOutput(new OutputHandlerEventArgs("Sucessfully trained IdentityCalibrationFunction"));
        }
    }
}