# AutoVent
Generate automated vent holes in STL geometry for vacuum forming workflows.

## 🚀 Features
- Detects surface depressions in STL geometry
- Automatically places **conical holes (frustums)** at local minima
- Configurable hole size, angle, spacing, and depth
- Uses **OpenSCAD** for robust frustum generation
- Outputs STL files as negative volumes for slicer

## 🛠 Requirements
- Python 3.9+
- [OpenSCAD](https://openscad.org/) (must be installed and accessible by path)
- Python libraries:
  - `trimesh`
  - `numpy`
  - `scipy`
  - `shapely`
  - `scikit-learn`

Install dependencies with:
```bash
sudo apt update
sudo apt install openscad python3-pip
pip install trimesh numpy scipy shapely scikit-learn
