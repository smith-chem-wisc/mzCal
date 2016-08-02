namespace mzCal
{
    internal class OppositeCalibrationFunction : CalibrationFunction
    {
        private CalibrationFunction combinedCalibration;

        public OppositeCalibrationFunction(CalibrationFunction combinedCalibration)
        {
            this.combinedCalibration = combinedCalibration;
        }

        public override double Predict(double[] t)
        {
            return -combinedCalibration.Predict(t);
        }
    }
}