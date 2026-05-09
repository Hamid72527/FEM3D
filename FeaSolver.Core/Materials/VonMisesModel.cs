using MathNet.Numerics.LinearAlgebra;
using FeaSolver.Core.Models;
using System;

namespace FeaSolver.Core.Materials
{
    public class VonMisesModel
    {
        public static (Vector<double> stress, Vector<double> backStress, double peeq, Matrix<double> ddsdde) UpdateStress(
            Material props,
            Vector<double> strainInc,
            Vector<double> prevStress,
            Vector<double> prevBackStress,
            double prevPeeq)
        {
            double E = props.E;
            double nu = props.Nu;
            double yield0 = props.YieldStress;
            double H = props.HardeningModule;

            var D = Kinematics.GetElasticDMatrix(E, nu);
            var stressTrial = prevStress + D * strainInc;

            // Deviatoric operator
            var DEV = Matrix<double>.Build.Dense(6, 6);
            DEV[0, 0] = DEV[1, 1] = DEV[2, 2] = 2.0/3.0;
            DEV[0, 1] = DEV[0, 2] = DEV[1, 0] = DEV[1, 2] = DEV[2, 0] = DEV[2, 1] = -1.0/3.0;
            DEV[3, 3] = DEV[4, 4] = DEV[5, 5] = 1.0;

            var stressDev = DEV * (stressTrial - prevBackStress);
            double qTrial = Math.Sqrt(1.5 * stressDev.DotProduct(DEV_Scale(stressDev)));

            double f = qTrial - (yield0 + H * prevPeeq);

            if (f <= 0)
            {
                return (stressTrial, prevBackStress, prevPeeq, D);
            }
            else
            {
                // Simple Radial Return
                double deltaP = f / (3 * props.G + H);
                var flowDir = (1.5 / qTrial) * stressDev;
                var stress = stressTrial - 2 * props.G * deltaP * flowDir;
                double peeq = prevPeeq + deltaP;

                // Tangent (Simplified for small strain)
                var tangent = D; // Should be elasto-plastic tangent for full NR
                return (stress, prevBackStress, peeq, tangent);
            }
        }

        private static Vector<double> DEV_Scale(Vector<double> v)
        {
            var res = v.Clone();
            res[3] /= 1.0; // In my B matrix and D matrix, engineering shear is used?
            // MATLAB uses [s11 s22 s33 s12 s23 s13]
            return res;
        }
    }
}
