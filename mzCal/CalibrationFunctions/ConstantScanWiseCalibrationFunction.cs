using System;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    public class ConstantScanWiseCalibrationFunction : CalibrationFunction
    {
        private Action<OutputHandlerEventArgs> onOutput;
        private MyBinaryTree myBinaryTree;

        public ConstantScanWiseCalibrationFunction(Action<OutputHandlerEventArgs> onOutput, IEnumerable<LabeledDataPoint> trainingList)
        {
            this.onOutput = onOutput;
            myBinaryTree  = new MyBinaryTree(trainingList.ToList());

        }

        public override double Predict(double[] inputs)
        {
            return myBinaryTree.GetValue(inputs[2]);
        }
    }
}