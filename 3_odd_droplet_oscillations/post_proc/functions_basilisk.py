import numpy as np
import matplotlib as mpl
import os
from scipy.interpolate import RegularGridInterpolator
from view_functions import *

def find_all_vtk_steps(results_folder):
    if( os.path.isdir(results_folder) ):
        list_files = os.listdir(results_folder)
        list_vtk_files = np.sort([int(file.split("Interface_")[1][:-4]) for file in list_files if file.startswith("Interface")])
    else:
        print("Folder not found: ", results_folder)
        exit()
    return list_vtk_files


def plot_mesh_basilisk_solution(file_name, property_name="pressure",
                                chosen_colormap=mpl.cm.coolwarm, crange=[None, None], ax=None,
                                plot_quiver = False):


    
    fields = ReadFields_Uniform(file_name, return_vertex_coordinates=True)

    basilisk_mesh_x = fields["x"]
    basilisk_mesh_y = fields["y"]
    basilisk_mesh_x_centers = fields["x-centers"]
    basilisk_mesh_y_centers = fields["y-centers"]
    basilisk_mesh_field = np.copy( fields[property_name] ) # making a copy

    
    # Taking only cells inside the droplet
    plot_only = basilisk_mesh_field > 1e+07
    basilisk_mesh_field[plot_only] = np.nan

    cb = None
    if( ax is not None ):
        # plot_only = np.sqrt(basilisk_mesh_x_centers**2 + basilisk_mesh_y_centers**2) >= 1.0
        cb = ax.pcolormesh(basilisk_mesh_x, basilisk_mesh_y, basilisk_mesh_field, cmap=chosen_colormap, vmin=crange[0], vmax=crange[1])

        if( plot_quiver ):
            # ax.quiver(basilisk_mesh_x_centers, basilisk_mesh_y_centers, basilisk_mesh_ux, basilisk_mesh_uy)
            basilisk_mesh_Ex = 0.5*basilisk_mesh_y_centers
            basilisk_mesh_Ey = 0.5*basilisk_mesh_x_centers
            basilisk_mesh_Rx = 0.5*basilisk_mesh_y_centers
            basilisk_mesh_Ry = - 0.5*basilisk_mesh_x_centers

            basilisk_mesh_ux -= basilisk_mesh_Rx
            basilisk_mesh_uy -= basilisk_mesh_Ry
            ax.streamplot(basilisk_mesh_x_centers, basilisk_mesh_y_centers, basilisk_mesh_ux, basilisk_mesh_uy)

    for field in fields:
        if( field!="x" and field!="y" and field!="x-centers" and field!="y-centers" ):
            fields[field][plot_only] = np.nan if field!="fractions" else 0.0

    return cb, fields



def get_deformation_and_angle(simulation_folder, color="black"):
    # Loading basilisk data
    list_vtk_files = None
    if( os.path.isdir(simulation_folder) ):
        list_files = os.listdir(simulation_folder)
        list_vtk_files = np.sort([int(file.split("MeshDump_")[1][:-4]) for file in list_files if file.startswith("MeshDump")])
    else:
        print("Folder not found")


    vtk_step = list_vtk_files[-1]
    file_name = "%s/Interface_%04d.vtk" % (simulation_folder, vtk_step)
    points, segments, _, line_coll, _, _, _ = \
        read_polydata(file_name, rotate=0.0, flip_y=False, color=color, linewidth=1.5, projected_rotation=False)
    centroid = np.average(points, axis=0)
    biggest_distance = 0.0
    farthest_point = None
    droplet_angle = None
    for p in points:
        p = p - centroid
        distance_to_center = np.linalg.norm(p)
        point_angle = np.rad2deg( np.arctan2(p[1], p[0]) )
        point_angle += 360.0 if (point_angle<0.0) else 0.0
        if( distance_to_center>biggest_distance and point_angle<=90 ): # Restricting to first quadrant
            biggest_distance = distance_to_center
            farthest_point = p
            droplet_angle = point_angle


    filedata = np.loadtxt("%s/log_file.txt" % (simulation_folder))
    array_deformation = filedata[:, 4]
    final_deformation = array_deformation[-1]

    return final_deformation, droplet_angle