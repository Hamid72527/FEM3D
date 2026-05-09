using Xunit;
using FeaSolver.Core.Elements;
using FeaSolver.Core.Models;
using MathNet.Numerics.LinearAlgebra;

namespace FeaSolver.Tests
{
    public class ElementTests
    {
        [Fact]
        public void Tet4_StiffnessMatrix_ShouldHaveCorrectSize()
        {
            var nodes = new[]
            {
                new Node(1, 0, 0, 0),
                new Node(2, 1, 0, 0),
                new Node(3, 0, 1, 0),
                new Node(4, 0, 0, 1)
            };
            var material = new Material(2e11, 0.3);
            var tet = new Tet4(1, new[] { 0, 1, 2, 3 });

            var K = tet.ComputeStiffness(nodes, material);

            Assert.Equal(12, K.RowCount);
            Assert.Equal(12, K.ColumnCount);
            // Volume of this tet is 1/6. DetJ should be 1.
            // Sum of stiffness matrix should be 0 (equilibrium)
            for(int i=0; i<12; i++)
            {
                Assert.True(System.Math.Abs(K.RowSums()[i]) < 1e-6);
            }
        }
    }
}
