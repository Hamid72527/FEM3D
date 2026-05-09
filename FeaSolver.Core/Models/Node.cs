namespace FeaSolver.Core.Models
{
    public class Node
    {
        public int Id { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public int[] DofIndices { get; set; } = new int[3];

        public Node(int id, double x, double y, double z)
        {
            Id = id;
            X = x;
            Y = y;
            Z = z;
        }
    }
}
