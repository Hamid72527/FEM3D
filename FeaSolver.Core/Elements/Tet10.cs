using MathNet.Numerics.LinearAlgebra;
using FeaSolver.Core.Models;
using FeaSolver.Core.Math;

namespace FeaSolver.Core.Elements
{
    public class Tet10 : IElement
    {
        public int Id { get; }
        public int[] NodeIndices { get; }

        public Tet10(int id, int[] nodeIndices)
        {
            Id = id;
            NodeIndices = nodeIndices;
        }

        public (double[] shape, double[] dhdr, double[] dhds, double[] dhdt) GetShapeFunctions(double r, double s, double t)
        {
            double l1 = 1 - r - s - t;
            double l2 = r;
            double l3 = s;
            double l4 = t;

            double[] shape = new double[10];
            shape[0] = l1 * (2 * l1 - 1.0);
            shape[1] = l2 * (2 * l2 - 1.0);
            shape[2] = l3 * (2 * l3 - 1.0);
            shape[3] = l4 * (2 * l4 - 1.0);
            shape[4] = 4.0 * l2 * l1;
            shape[5] = 4.0 * l2 * l3;
            shape[6] = 4.0 * l3 * l1;
            shape[7] = 4.0 * l4 * l1;
            shape[8] = 4.0 * l2 * l4;
            shape[9] = 4.0 * l3 * l4;

            double[] dhdr = new double[10];
            dhdr[0] = -4 * l1 + 1.0;
            dhdr[1] = 4 * l2 - 1.0;
            dhdr[4] = 4 * l1 - 4 * l2;
            dhdr[5] = 4 * l3;
            dhdr[6] = -4 * l3;
            dhdr[7] = -4 * l4;
            dhdr[8] = 4 * l4;

            double[] dhds = new double[10];
            dhds[0] = -4 * l1 + 1.0;
            dhds[2] = 4 * l3 - 1.0;
            dhds[4] = -4 * l2;
            dhds[5] = 4 * l2;
            dhds[6] = 4 * l1 - 4 * l3;
            dhds[7] = -4 * l4;
            dhds[9] = 4 * l4;

            double[] dhdt = new double[10];
            dhdt[0] = -4 * l1 + 1.0;
            dhdt[3] = 4 * l4 - 1.0;
            dhdt[4] = -4 * l2;
            dhdt[6] = -4 * l3;
            dhdt[7] = 4 * l1 - 4 * l4;
            dhdt[8] = 4 * l2;
            dhdt[9] = 4 * l3;

            return (shape, dhdr, dhds, dhdt);
        }

        public Matrix<double> ComputeStiffness(Node[] allNodes, Material material)
        {
            int nnel = 10;
            int edof = nnel * 3;
            var K = Matrix<double>.Build.Dense(edof, edof);

            var (points, weights) = GaussQuadrature.GetTetPoints(4);
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
