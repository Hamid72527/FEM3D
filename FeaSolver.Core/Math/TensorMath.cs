using MathNet.Numerics.LinearAlgebra;
using System;

namespace FeaSolver.Core.Math
{
    public static class TensorMath
    {
        public static Matrix<double> ToMatrix(Vector<double> v)
        {
            var m = Matrix<double>.Build.Dense(3, 3);
            m[0, 0] = v[0];
            m[1, 1] = v[1];
            m[2, 2] = v[2];
            m[0, 1] = m[1, 0] = v[3];
            m[1, 2] = m[2, 1] = v[4];
            m[0, 2] = m[2, 0] = v[5];
            return m;
        }

        public static (double J2, double J3) GetInvariants(Vector<double> stress)
        {
            double mean = (stress[0] + stress[1] + stress[2]) / 3.0;
            double s11 = stress[0] - mean;
            double s22 = stress[1] - mean;
            double s33 = stress[2] - mean;
            double s12 = stress[3];
            double s23 = stress[4];
            double s13 = stress[5];

            double J2 = 0.5 * (s11 * s11 + s22 * s22 + s33 * s33) + s12 * s12 + s23 * s23 + s13 * s13;

            var sMat = Matrix<double>.Build.Dense(3, 3);
            sMat[0, 0] = s11; sMat[1, 1] = s22; sMat[2, 2] = s33;
            sMat[0, 1] = sMat[1, 0] = s12;
            sMat[1, 2] = sMat[2, 1] = s23;
            sMat[0, 2] = sMat[2, 0] = s13;

            double J3 = sMat.Determinant();

            return (J2, J3);
        }
    }
}
