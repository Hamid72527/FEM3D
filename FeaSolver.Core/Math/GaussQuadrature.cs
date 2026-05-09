namespace FeaSolver.Core.Math
{
    public static class GaussQuadrature
    {
        public static (double[,] points, double[] weights) GetTetPoints(int ngl)
        {
            if (ngl == 1)
            {
                return (new double[,] { { 0.25, 0.25, 0.25, 0.25 } }, new double[] { 1.0 });
            }
            if (ngl == 4)
            {
                double l1 = 0.585410196624968;
                double l2 = 0.138196601125010;
                double w = 0.25;
                return (new double[,] {
                    { l1, l2, l2, l2 },
                    { l2, l1, l2, l2 },
                    { l2, l2, l1, l2 },
                    { l2, l2, l2, l1 }
                }, new double[] { w, w, w, w });
            }
            throw new System.NotImplementedException("Only 1 and 4 point rules are implemented.");
        }
    }
}
