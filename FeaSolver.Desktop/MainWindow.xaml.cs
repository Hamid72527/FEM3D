using System;
using System.Windows;
using HelixToolkit.Wpf;
using System.Windows.Media.Media3D;
using Microsoft.Win32;

namespace FeaSolver.Desktop
{
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void ImportSTL_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new OpenFileDialog { Filter = "STL files (*.stl)|*.stl" };
            if (dialog.ShowDialog() == true)
            {
                var reader = new StlReader();
                var modelGroup = reader.Read(dialog.FileName);
                GeometryModel.Content = modelGroup;
                Viewport.ZoomExtents();
            }
        }

        private void GenerateMesh_Click(object sender, RoutedEventArgs e)
        {
            MessageBox.Show("Meshing logic will be triggered here.");
        }

        private void RunSim_Click(object sender, RoutedEventArgs e)
        {
            MessageBox.Show("Solver logic will be triggered here.");
        }
    }
}
