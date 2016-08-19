using System;
using System.Collections.Generic;

namespace mzCal
{
    internal class ByHandCalibrationFunction : CalibrationFunction
    {
        public ByHandCalibrationFunction(Action<OutputHandlerEventArgs> onOutput, List<LabeledDataPoint> trainList1)
        {
        }

        public override double Predict(double[] t)
        {
            return -t[1] / 200000;
        }
    }
}