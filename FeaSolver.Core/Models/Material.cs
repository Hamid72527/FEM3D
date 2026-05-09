namespace FeaSolver.Core.Models
{
    public class Material
    {
        public string Name { get; set; } = "Steel";
        public double E { get; set; } // Elastic Modulus
        public double Nu { get; set; } // Poisson's Ratio
        public double YieldStress { get; set; }
        public double HardeningModule { get; set; }

        public double G => E / (2 * (1 + Nu));
        public double K => E / (3 * (1 - 2 * Nu));

        public Material(double e, double nu, double yieldStress = 0, double h = 0)
        {
            E = e;
            Nu = nu;
            YieldStress = yieldStress;
            HardeningModule = h;
        }
    }
}
