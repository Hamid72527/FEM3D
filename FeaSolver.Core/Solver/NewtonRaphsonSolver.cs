using MathNet.Numerics.LinearAlgebra;
using FeaSolver.Core.Models;
using FeaSolver.Core.Elements;
using System.Collections.Generic;
using System.Linq;

namespace FeaSolver.Core.Solver
{
    public class NewtonRaphsonSolver
    {
        public static Vector<double> Solve(Mesh mesh, Material material, Dictionary<int, double> dirichletBCs, Vector<double> externalForce)
        {
            int sdof = mesh.Nodes.Length * 3;
            var displ = Vector<double>.Build.Dense(sdof);
            var K_global = Matrix<double>.Build.Dense(sdof, sdof);

            // For simple linear case (first iteration of NR)
            foreach (var elem in mesh.Elements)
            {
                var ke = elem.ComputeStiffness(mesh.Nodes, material);
                Assemble(K_global, ke, elem.NodeIndices);
            }

            var f_global = externalForce.Clone();
            ApplyBCs(K_global, f_global, dirichletBCs);

            displ = K_global.Solve(f_global);
            return displ;
        }

        private static void Assemble(Matrix<double> K, Matrix<double> ke, int[] nodeIndices)
        {
            int nnel = nodeIndices.Length;
            for (int i = 0; i < nnel; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    int row = nodeIndices[i] * 3 + j;
                    for (int k = 0; k < nnel; k++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            int col = nodeIndices[k] * 3 + l;
                            K[row, col] += ke[i * 3 + j, k * 3 + l];
                        }
                    }
                }
            }
        }

        private static void ApplyBCs(Matrix<double> K, Vector<double> f, Dictionary<int, double> dirichletBCs)
        {
            foreach (var bc in dirichletBCs)
            {
                int dof = bc.Key;
                double val = bc.Value;

                for (int j = 0; j < K.ColumnCount; j++) K[dof, j] = 0;
                for (int i = 0; i < K.RowCount; i++) K[i, dof] = 0;

                K[dof, dof] = 1.0;
                f[dof] = val;
            }
        }
    }
}

namespace FeaSolver.Core.Models
{
    public class Mesh
    {
        public Node[] Nodes { get; set; }
        public IElement[] Elements { get; set; }

        public Mesh(Node[] nodes, IElement[] elements)
        {
            Nodes = nodes;
            Elements = elements;
        }
    }
}
