# Running FeaSolver in Visual Studio

1. **Prerequisites**:
   - Install Visual Studio 2022 (latest version).
   - Ensure the ".NET desktop development" workload is installed.
   - Ensure you have .NET 8.0 or 10.0 SDK installed.

2. **Opening the Project**:
   - Locate the `FeaSolver.sln` file in the root directory.
   - Double-click it to open in Visual Studio.

3. **Restoring Packages**:
   - Visual Studio should automatically restore NuGet packages (`MathNet.Numerics` and `HelixToolkit.Wpf`).
   - If not, right-click the solution in Solution Explorer and select "Restore NuGet Packages".

4. **Running the App**:
   - Set `FeaSolver.Desktop` as the Startup Project (right-click it > Set as Startup Project).
   - Press **F5** or click the "Start" button (Green arrow) to run.

5. **Using the App**:
   - Click "Import STL" and select the provided `Part.stl`.
   - Use the mouse to rotate/zoom in the 3D viewport.
   - Click "Generate Mesh" to create the tetrahedral mesh.
   - Select material and click "Run Simulation" to solve.
