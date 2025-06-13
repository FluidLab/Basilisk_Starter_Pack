import numpy as np
import matplotlib.pyplot as plt
import os
from functions_basilisk import *
plt.rcParams.update({'font.size': 15})

Oh_even = 0.01
Oh_odd = 0.1
mesh_level = 9
dt = 0.0001
base_folder = "D:/share/basilisk_codes/git_stuff/Basilisk_Starter_Pack/3_odd_droplet_oscillations/outputs/"
simulation_folder = "odd_oscillations_DT%g_mesh%d_OhE%g_OhO%g/" % (dt, mesh_level, Oh_even, Oh_odd)

list_vtk_steps = find_all_vtk_steps(base_folder + simulation_folder)

frames_folder = "frames"
if( not os.path.isdir(frames_folder) ):
    os.mkdir(frames_folder)
# os.system("rm frames/*.png")

# Loading the log_file which constains data over time like the average vorticity, droplet size, etc
filedata = np.loadtxt("%s/log_file.txt" % (base_folder + simulation_folder))
log_array_time = filedata[:, 1]
log_array_avg_vorticity = filedata[:, 8]

i_frame = 0
for vtk_step in list_vtk_steps:

    fig, (ax, ax_curve) = plt.subplots(1, 2, figsize=(12, 8))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim([-1.0, 1.0])
    ax.set_ylim([-1.0, 1.0])
    ax.tick_params(left = False, labelleft = False, labelbottom = False, bottom = False)

    # Plotting the interface of the droplet
    file_name = "%s/Interface_%04d.vtk" % (base_folder + simulation_folder, vtk_step)
    time, points, segments, normals, line_collection, _ = read_polydata(file_name, color="black", linewidth=3)
    ax.add_collection(line_collection)

    # Plotting a scalar field inside the droplet
    
    file_name = "%s/MeshDump_%04d.bin" % (base_folder + simulation_folder, vtk_step)
    if( os.path.isfile(file_name) ):
        cb, fields = plot_mesh_basilisk_solution(file_name, property_name="vorticity",
                                        chosen_colormap=mpl.cm.coolwarm, crange=[-2.0, 2.0], ax=ax,
                                        plot_quiver = False)
        fig.colorbar(cb, ax=ax, orientation="horizontal", shrink=0.8, label="Vorticity", pad=0.03)

    # Adding a text label showing the time
    ax.text(0.02, 0.98, "Time = %.2f" % (time), ha="left", va="top", transform=ax.transAxes, fontsize=17)

    # Plotting the curve with the data from the log_file
    ax_curve.plot(log_array_time, log_array_avg_vorticity, "-r", lw=3)
    ax_curve.axvline(x=time, ls="--", color="gray", lw=2.0)
    ax_curve.axhline(y=0.0, ls="--", color="gray", lw=1.5)
    ax_curve.set_xlabel("Time")
    ax_curve.set_ylabel("Average vorticity")

    plt.tight_layout()
    plt.savefig("%s/frame%04d.png" % (frames_folder, i_frame))
    plt.close()
    print("Made video frame %d out of %d ... " % (vtk_step, list_vtk_steps[-1]))
    i_frame += 1
    # exit()

ffmpeg_folder = ""
video_name = "video_OhEven%g_OhOdd%g.mp4" % (Oh_even, Oh_odd)
framerate = 10
os.system('%sffmpeg -y -framerate %s -i %s/frame%%04d.png -b 5000k %s' % (ffmpeg_folder, str(framerate), frames_folder, video_name))