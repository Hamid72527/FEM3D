using MathNet.Numerics.LinearAlgebra;

namespace FeaSolver.Core.Math
{
    public static class Kinematics
    {
        public static Matrix<double> ComputeBMatrix(int nnel, double[] dhdx, double[] dhdy, double[] dhdz)
        {
            var B = Matrix<double>.Build.Dense(6, nnel * 3);
            for (int i = 0; i < nnel; i++)
            {
                int i1 = i * 3;
                int i2 = i1 + 1;
                int i3 = i2 + 1;

                B[0, i1] = dhdx[i];
                B[1, i2] = dhdy[i];
                B[2, i3] = dhdz[i];

                B[3, i1] = dhdy[i];
                B[3, i2] = dhdx[i];

                B[4, i2] = dhdz[i];
                B[4, i3] = dhdy[i];

                B[5, i1] = dhdz[i];
                B[5, i3] = dhdx[i];
            }
            return B;
        }

        public static Matrix<double> GetElasticDMatrix(double E, double nu)
        {
            double factor = E / ((1 + nu) * (1 - 2 * nu));
            var D = Matrix<double>.Build.Dense(6, 6);
            D[0, 0] = D[1, 1] = D[2, 2] = 1 - nu;
            D[0, 1] = D[0, 2] = D[1, 0] = D[1, 2] = D[2, 0] = D[2, 1] = nu;
            D[3, 3] = D[4, 4] = D[5, 5] = (1 - 2 * nu) / 2.0;
            return D * factor;
        }
    }
}
