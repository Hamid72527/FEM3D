using MathNet.Numerics.LinearAlgebra;
using FeaSolver.Core.Models;
using FeaSolver.Core.Math;
using System;
using System.Linq;

namespace FeaSolver.Core.Materials
{
    public class TrescaModel
    {
        public static (Vector<double> stress, Vector<double> backStress, double peeq, Matrix<double> ddsdde) UpdateStress(
            Material props,
            Vector<double> strainInc,
            Vector<double> prevStress,
            Vector<double> prevBackStress,
            double prevPeeq)
        {
            var D = Kinematics.GetElasticDMatrix(props.E, props.Nu);
            var stressTrial = prevStress + D * strainInc;

            var sMat = TensorMath.ToMatrix(stressTrial - prevBackStress);
            var evals = sMat.Evd().EigenValues.Select(e => e.Real).OrderBy(v => v).ToArray();

            double trescaTrial = evals[2] - evals[0];
            double yield = props.YieldStress + props.HardeningModule * prevPeeq;

            if (trescaTrial <= yield)
            {
                return (stressTrial, prevBackStress, prevPeeq, D);
            }
            else
            {
                // Simplified return map for Tresca
                double deltaP = (trescaTrial - yield) / (4 * props.G + props.HardeningModule);
                // This is an approximation. Full Tresca return map is complex.

                // For the sake of this port, we will use a radial-like return
                double scale = yield / trescaTrial;
                var stress = prevBackStress + scale * (stressTrial - prevBackStress);

                return (stress, prevBackStress, prevPeeq + deltaP, D);
            }
        }
    }
}
