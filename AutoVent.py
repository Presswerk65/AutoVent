import trimesh
import numpy as np
import os
import sys
import shutil
import subprocess
from glob import glob
from sklearn.neighbors import KDTree
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point


# === CONFIGURATION ===

drill_tip_diameter_mm = 1.2       # Diameter at the tip of the frustum (mm)
drill_angle_deg = 15.0            # Opening angle of the frustum (degrees)
frustum_height = 10               # Lentgh of the frustum (mm)

local_min_radius_mm = 1.2         # Radius (mm) within which a point is considered a local minimum
min_hole_spacing_mm = 5.0         # Minimum distance (mm) between two holes at the same height

z_tolerance_mm = 0.4              # Max. height difference (mm) for points to be considered equal
elevation_threshold_mm = 1.0      # Height difference (mm) used to check elevations between points
intermediate_points = 5           # Number of intermediate points for elevation check between two points

grid_size_mm = 0.15               # Grid size (mm) for XY data reduction

edge_clearance_mm = 2             # Minimum clearance (mm) from the boundary for hole placement
surface_depth_mm = 10.0           # Depth (mm) from the top surface for hole detection

# Automatically set working directory to current directory
working_dir = os.getcwd()

# Determine OpenSCAD path depending on platform
if sys.platform.startswith("win"):
    openscad_path = r"C:\Program Files\OpenSCAD\openscad.exe"
else:
    openscad_path = "openscad"  # Linux / Codespace: OpenSCAD must be in PATH

# Subfolders for input, export, and temporary files
input_dir = working_dir
export_dir = os.path.join(working_dir, "export")
tmp_dir = os.path.join(working_dir, "temp")

# Ensure all required directories exist
os.makedirs(input_dir, exist_ok=True)
os.makedirs(export_dir, exist_ok=True)
os.makedirs(tmp_dir, exist_ok=True)

def generate_frustum_openscad(radius_top, radius_base, height, path):
    """Generate a frustum in OpenSCAD and export as STL."""
    scad_code = f"""
    cylinder(h={height}, r1={radius_base}, r2={radius_top}, $fn=64);
    """
    scad_file = path.replace(".stl", ".scad")
    with open(scad_file, "w") as f:
        f.write(scad_code)

    result = subprocess.run([openscad_path, "-o", path, scad_file], capture_output=True, text=True)
    if result.returncode != 0:
        print("OpenSCAD error:", result.stderr)
        return False
    return True


def process_stl(path):
    """Process a single STL file, detect depressions, and generate drill frustums."""
    name = os.path.basename(path)
    basename, ext = os.path.splitext(name)
    output_path = os.path.join(export_dir, f"{basename}_Holes{ext}")
    print(f"\n🔧 Processing: {name}")

    mesh = trimesh.load(path)
    vertices = mesh.vertices
    z_min, z_max = mesh.bounds[:, 2]

    # Select only top surface points
    top_points = vertices[vertices[:, 2] >= z_max - surface_depth_mm]

    # Create convex hull for boundary detection
    hull = ConvexHull(vertices[:, :2])
    hull_polygon = Polygon(vertices[hull.vertices, :2])

    # Reduce points to grid
    xy_grid = np.floor(top_points[:, :2] / grid_size_mm) * grid_size_mm
    _, unique_idx = np.unique(xy_grid, axis=0, return_index=True)
    sampled_points = top_points[unique_idx]

    # Filter points near the edge
    mask = []
    for pt in sampled_points:
        p = Point(pt[:2])
        mask.append(hull_polygon.exterior.distance(p) >= edge_clearance_mm)
    sampled_points = sampled_points[mask]

    # Detect local minima
    tree = KDTree(sampled_points[:, :2])
    local_minima = []
    for i, pt in enumerate(sampled_points):
        idx = tree.query_radius([pt[:2]], r=local_min_radius_mm)[0]
        if all(pt[2] <= sampled_points[j][2] for j in idx if j != i):
            local_minima.append(pt)

    print(f"→ {len(local_minima)} depressions found – filtering by spacing & elevation…")

    # Filter based on spacing and terrain elevation
    filtered_points = []
    point_tree = KDTree(sampled_points[:, :2])
    for pt in local_minima:
        allowed = True
        for f in filtered_points:
            dist = np.linalg.norm(f[:2] - pt[:2])
            if dist < min_hole_spacing_mm and abs(f[2] - pt[2]) < z_tolerance_mm:
                line = np.linspace(f[:2], pt[:2], num=intermediate_points)
                _, idxs = point_tree.query(line, k=1)
                max_z = max(f[2], pt[2])
                if all(sampled_points[i][2] <= max_z + elevation_threshold_mm for i in idxs.flatten()):
                    allowed = False
                    break
        if allowed:
            filtered_points.append(pt)

    print(f"✅ {len(filtered_points)} holes will be placed")

    # Generate frustums
    holes = []
    tip_radius = drill_tip_diameter_mm / 2.0
    opening_angle_rad = np.deg2rad(drill_angle_deg)

    for i, pt in enumerate(filtered_points):
        base_radius = tip_radius + np.tan(opening_angle_rad / 2.0) * frustum_height
        frustum_path = os.path.join(tmp_dir, f"frustum_{i}.stl")
        ok = generate_frustum_openscad(tip_radius, base_radius, frustum_height, frustum_path)
        if not ok:
            print(f"⚠️ Error generating frustum {i} with OpenSCAD")
            continue
        frustum = trimesh.load(frustum_path)
        frustum.apply_translation([pt[0], pt[1], pt[2] - frustum_height + 0.1])
        holes.append(frustum)

    if not holes:
        print("⚠️ No holes generated. Copying original file.")
        shutil.copy(path, output_path)
        return

    # Clean temporary directory
    for f in os.listdir(tmp_dir):
        os.remove(os.path.join(tmp_dir, f))

    all_holes = trimesh.util.concatenate(holes)
    all_holes.export(output_path)

    print(f"💾 Frustum saved in: {output_path}")

    # Clean temporary directory again
    # for f in os.listdir(tmp_dir):
    #    os.remove(os.path.join(tmp_dir, f))

    # Clean up temp directory after processing
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir, ignore_errors=True)        


# Process all STL files in the current directory
all_stl_files = glob("*.stl")
if not all_stl_files:
    print("📂 No STL files found in the current directory.")
else:
    for path in all_stl_files:
        process_stl(path)
    print("\n🏁 All STL files processed.")
