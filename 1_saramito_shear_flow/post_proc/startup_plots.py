import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import solve_ivp
font = {'weight' : 'normal', 'size'   : 20}
matplotlib.rc('font', **font)

# Base folder where Basilisk is outputting the simulation data in your machine
base_folder = "/mnt/lustre/home/hugof/spreading_paper_code/shear_flow/outputs/"

# Parameters to plot. Doing a sweep in Wi for fixed Bi, beta, Re
# Note: make sure you ran the simulation for these parameters before using this script
# The script will look for the simulated data
array_Wi = [0.05, 0.1, 0.5, 1]
Bi = 3
beta = 0.111111
Re = 0.1
array_colors = ["red", "green", "blue", "orange", "purple", "black"]

# right hand side (RHS) of Saramito's ODE for the stress (equation 13 in his paper)
# Note: I leave on the time derivative on the left side. Everything else is on the RHS
# Note2: I assume a = 1 in Saramito's equation 13 in the paper
def ode_rhs(time, tau):
	txx = tau[0]
	txy = tau[1]

	norm_tau = np.sqrt( 0.5*txx*txx + 2.0*txy*txy ) + 1e-6
	k = np.max([0.0, 1.0 - Bi/norm_tau])
	rhs_tau11 = 2.0*txy - (k/Wi)*txx
	rhs_tau12 = - (k/Wi)*txy + (1.0-beta)/Wi

	return np.array( [rhs_tau11, rhs_tau12] )



fig, ax = plt.subplots(2, 1, figsize=(12, 9))
ax_txy = ax[0]
ax_txx = ax[1]
for i_Wi, Wi in enumerate(array_Wi):
	print("Processing data for Wi = %g" % (Wi))

	# === Exact solution
	tMin = 0.0001
	tMax = 100.0
	t_eval = np.logspace(np.log10(tMin), np.log10(tMax), 100)
	solution = solve_ivp(ode_rhs, [0.0, tMax], [0.0, 0.0], t_eval=t_eval)
	ax_txx.loglog(solution.t, solution.y[0], color=array_colors[i_Wi])
	ax_txy.loglog(solution.t, solution.y[1], color=array_colors[i_Wi])

	# Selecting some equally spaced points in the log scale
	def find_nearest(x):
		difference_array = np.absolute(array_t - x)
		return difference_array.argmin() #returns the index


	# Numerical solution
	simulation_folder = "startup_Wi%g_Bi%g_beta%g_Re%g/" % (Wi, Bi, beta, Re)
	file_name = base_folder + simulation_folder + "log_file.txt"
	file_data = np.loadtxt(file_name)
	array_t = file_data[:, 1]
	array_txx = file_data[:, 3]
	array_txy = file_data[:, 4]
	array_tyy = file_data[:, 5]
	selected_indices = np.array( list(map( find_nearest, t_eval )) )
	array_t = array_t[selected_indices]
	array_txx = array_txx[selected_indices]
	array_txy = array_txy[selected_indices]
	array_tyy = array_tyy[selected_indices]
	ax_txx.loglog(array_t, array_txx, 'o', color=array_colors[i_Wi], fillstyle="none")
	ax_txy.loglog(array_t, array_txy, 'o', label = "Wi = %g" % (Wi), color=array_colors[i_Wi], fillstyle="none")
	

ax_txy.legend()
ax_txy.set_xlim([tMin, 5.0*2.0*np.pi])
ax_txy.set_xlabel(r"Time: $t$ [nondim]")
ax_txy.set_ylabel(r"Stress: $\tau_{xy}$ [nondim]")
	

ax_txx.set_xlim([tMin, 5.0*2.0*np.pi])
ax_txx.set_xlabel(r"Time: $t$ [nondim]")
ax_txx.set_ylabel(r"Stress: $\tau_{xx}$ [nondim]")

# print( ax_txy.get_ylim() )
# print( ax_txx.get_ylim() )

# ax_txy.set_ylim((0.0005814090331328178, 4.516181846928209))
# ax_txx.set_ylim((1.7801207106242702e-08, 4.581777118489361))


plt.suptitle(r"Fixed: $Bi = %g$       $\beta = %g$   %s" % (Bi, beta, "(viscoelastic)" if Bi==0 else ""))
plt.tight_layout()
plt.savefig("startup_Bi%g.png" % (Bi))

print("\n\nFinished. A figure was created in the current folder.")