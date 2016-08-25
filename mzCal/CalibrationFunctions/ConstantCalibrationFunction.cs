﻿using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    public class ConstantCalibrationFunction : CalibrationFunction
    {
        public double a;
        private Action<OutputHandlerEventArgs> onOutput;

        public ConstantCalibrationFunction(Action<OutputHandlerEventArgs> onOutput)
        {
            this.onOutput = onOutput;
        }

        public override double Predict(double[] inputs)
        {
            return a;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            a = trainingList.Select(b => b.output).Median();
            onOutput(new OutputHandlerEventArgs("a = " + a));
        }
    }
}