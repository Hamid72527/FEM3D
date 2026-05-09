using MathNet.Numerics.LinearAlgebra;
using FeaSolver.Core.Models;

namespace FeaSolver.Core.Elements
{
    public interface IElement
    {
        int Id { get; }
        int[] NodeIndices { get; }
        Matrix<double> ComputeStiffness(Node[] allNodes, Material material);
    }
}
