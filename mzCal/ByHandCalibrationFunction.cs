using System;
using System.Collections.Generic;

namespace mzCal
{
    internal class ByHandCalibrationFunction : CalibrationFunction
    {
        private Action<OutputHandlerEventArgs> onOutput;

        public ByHandCalibrationFunction(Action<OutputHandlerEventArgs> onOutput, List<LabeledDataPoint> trainList1)
        {
            this.onOutput = onOutput;
        }

        public override double Predict(double[] t)
        {
            return -t[1] / 200000;
        }
    }
}