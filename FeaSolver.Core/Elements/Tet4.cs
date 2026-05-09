using MathNet.Numerics.LinearAlgebra;
using FeaSolver.Core.Models;
using FeaSolver.Core.Math;

namespace FeaSolver.Core.Elements
{
    public class Tet4 : IElement
    {
        public int Id { get; }
        public int[] NodeIndices { get; }

        public Tet4(int id, int[] nodeIndices)
        {
            Id = id;
            NodeIndices = nodeIndices;
        }

        public (double[] shape, double[] dhdr, double[] dhds, double[] dhdt) GetShapeFunctions(double r, double s, double t)
        {
            double[] shape = { 1 - r - s - t, r, s, t };
            double[] dhdr = { -1, 1, 0, 0 };
            double[] dhds = { -1, 0, 1, 0 };
            double[] dhdt = { -1, 0, 0, 1 };
            return (shape, dhdr, dhds, dhdt);
        }

        public Matrix<double> ComputeStiffness(Node[] allNodes, Material material)
        {
            int nnel = 4;
            int edof = nnel * 3;
            var K = Matrix<double>.Build.Dense(edof, edof);

            var (points, weights) = GaussQuadrature.GetTetPoints(1);
            var D = Kinematics.GetElasticDMatrix(material.E, material.Nu);

            for (int i = 0; i < weights.Length; i++)
            {
                double r = points[i, 0];
                double s = points[i, 1];
                double t = points[i, 2];
                var (_, dhdr, dhds, dhdt) = GetShapeFunctions(r, s, t);

                var jacob = Matrix<double>.Build.Dense(3, 3);
                for (int j = 0; j < nnel; j++)
                {
                    var node = allNodes[NodeIndices[j]];
                    jacob[0, 0] += dhdr[j] * node.X;
                    jacob[0, 1] += dhdr[j] * node.Y;
                    jacob[0, 2] += dhdr[j] * node.Z;
                    jacob[1, 0] += dhds[j] * node.X;
                    jacob[1, 1] += dhds[j] * node.Y;
                    jacob[1, 2] += dhds[j] * node.Z;
                    jacob[2, 0] += dhdt[j] * node.X;
                    jacob[2, 1] += dhdt[j] * node.Y;
                    jacob[2, 2] += dhdt[j] * node.Z;
                }

                double detJ = jacob.Determinant();
                var invJ = jacob.Inverse();

                double[] dhdx = new double[nnel];
                double[] dhdy = new double[nnel];
                double[] dhdz = new double[nnel];

                for (int j = 0; j < nnel; j++)
                {
                    dhdx[j] = invJ[0, 0] * dhdr[j] + invJ[0, 1] * dhds[j] + invJ[0, 2] * dhdt[j];
                    dhdy[j] = invJ[1, 0] * dhdr[j] + invJ[1, 1] * dhds[j] + invJ[1, 2] * dhdt[j];
                    dhdz[j] = invJ[2, 0] * dhdr[j] + invJ[2, 1] * dhds[j] + invJ[2, 2] * dhdt[j];
                }

                var B = Kinematics.ComputeBMatrix(nnel, dhdx, dhdy, dhdz);
                K += B.Transpose() * D * B * (weights[i] * detJ / 6.0);
            }
            return K;
        }
    }
}
