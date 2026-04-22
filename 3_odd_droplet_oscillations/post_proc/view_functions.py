import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection, PolyCollection
import matplotlib as mpl
import cmasher as cmr
import struct
from scipy.interpolate import Rbf



def ReadLine_Binary(file):
	line = ""
	read_stuff = file.read(1).decode('UTF-8')
	index = 0
	while( read_stuff!='\n' ):
		index += 1
		line += read_stuff
		read_stuff = file.read(1).decode('UTF-8')
		if(index>500):
			return None
	return line

def ReadPoint_Double(file):
	x = struct.unpack('>d', file.read(8))[0]
	y = struct.unpack('>d', file.read(8))[0]
	z = struct.unpack('>d', file.read(8))[0]
	return np.array([x, y])

def ReadPoint_Float(file):
	x = struct.unpack('>f', file.read(4))[0]
	y = struct.unpack('>f', file.read(4))[0]
	z = struct.unpack('>f', file.read(4))[0]
	return np.array([x, y, z])

def ReadInt(file):
	return struct.unpack('>i', file.read(4))[0]

def ReadDouble(file):
	return struct.unpack('>d', file.read(8))[0]

def ReadFloat(file):
	return struct.unpack('>f', file.read(4))[0]

def ReadFieldsVTK(file_name, color_fields="pressure", 
				  cmap=cmr.iceburn, crange=[None, None], xlim=[None, None], ylim=[None, None]):

	ReadPoint = ReadPoint_Float
	ReadNumber = ReadFloat
	
	try:
		file = open(file_name, 'rb')
	except FileNotFoundError:
		return None, None, None, None, None, None
	
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)

	num_points = int(line.split(" ")[1])
	

	num_cells = np.rint(num_points/4).astype(int)
	num_triangles = 2*num_cells
	array_cell_center = np.zeros(shape=(num_cells, 2))
	# array_cell_sizes = np.zeros(shape=(num_cells,))
	
	array_triangles = np.zeros(shape=(num_triangles, 3))
	array_points = np.zeros(shape=(3*num_triangles, 2))

	

	array_facecolors = []
	for i in range(num_cells):
		p1 = ReadPoint(file)
		p2 = ReadPoint(file)
		p3 = ReadPoint(file)
		p4 = ReadPoint(file)
		cell_center = 0.25*(p1 + p2 + p3 + p4)
		# radius = np.linalg.norm(cell_center - p1)

		array_cell_center[i, :] = cell_center

		index = 6*i
		array_points[index, :] = p3
		array_points[index+1, :] = p2
		array_points[index+2, :] = p1

		array_points[index+3, :] = p4
		array_points[index+4, :] = p3
		array_points[index+5, :] = p1

		
		array_triangles[2*i, :] = np.array([index, index+1, index+2])
		array_triangles[2*i + 1, :] = np.array([index+3, index+4, index+5])

		# array_triangles[2*i, :] = np.array([index, index+1, index+2])
		# array_triangles[2*i + 1, :] = np.array([index, index+2, index+3])

		array_facecolors.append([1.0, 0.0, 0.0])
		array_facecolors.append([1.0, 0.0, 0.0])
		array_facecolors.append([1.0, 0.0, 0.0])
		array_facecolors.append([1.0, 0.0, 0.0])
		array_facecolors.append([1.0, 0.0, 0.0])
		array_facecolors.append([1.0, 0.0, 0.0])

	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	for i in range(num_cells):
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)

	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	for i in range(num_cells):
		read_stuff = ReadInt(file)

	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)

	# Starting to read printed field data
	values = {}
	line = ReadLine_Binary(file)
	far_away_pressure = 0.0
	num_cells_far_away_pressure = 0
	while( (line is not None) and line.split(" ")[0]=="SCALARS" ):
		field_name = line.split(" ")[1]
		line = ReadLine_Binary(file)
		values[field_name] = np.zeros(shape=(num_triangles,))
		index = num_triangles - 1
		for i in range(num_cells):
			read_stuff = ReadNumber(file)
			values[field_name][index] = read_stuff
			values[field_name][index - 1] = read_stuff
			index = index - 2

			cell_center = array_cell_center[i]
			radius = np.sqrt(cell_center[0]**2 + cell_center[1]**2)
			if( field_name=="pressure" and radius>9.95 ):
				far_away_pressure += read_stuff
				num_cells_far_away_pressure += 1

		line = ReadLine_Binary(file)
		line = ReadLine_Binary(file)

	far_away_pressure = far_away_pressure/num_cells_far_away_pressure if color_fields=="pressure" else 0.0


	Ca = 0.01
	laplace_pressure = 1.0/(Ca) if color_fields=="pressure" else 0.0
	# laplace_pressure = 0.0

	# Setting cell colours
	dict_values = values
	values = dict_values[color_fields] - far_away_pressure
	crange[0] = np.min(values) if crange[0] is None else crange[0]
	crange[1] = np.max(values - laplace_pressure) if crange[1] is None else crange[1]
	norm_colors = mpl.colors.Normalize(vmin=crange[0], vmax=crange[1]) 
	
	print("min max: ", crange, far_away_pressure, color_fields)
	for i in range(num_cells):
		index_vertex = 6*i
		

		index_triangle = 2*i
		current_volume_fraction = dict_values["fractions"][index_triangle]
		current_value = values[index_triangle] - laplace_pressure*current_volume_fraction
		color = cmap( norm_colors(current_value) )
		array_facecolors[index_vertex + 0] = color
		array_facecolors[index_vertex + 1] = color
		array_facecolors[index_vertex + 2] = color

		index_triangle = 2*i + 1
		current_volume_fraction = dict_values["fractions"][index_triangle]
		current_value = values[index_triangle] - laplace_pressure*current_volume_fraction
		color = cmap( norm_colors(current_value) )
		array_facecolors[index_vertex + 3] = color
		array_facecolors[index_vertex + 4] = color
		array_facecolors[index_vertex + 5] = color


	# Applying transformations
	# rotate = 0.5*np.pi
	# flip_x = True
	# A = np.array([[np.cos(rotate), -np.sin(rotate)], [np.sin(rotate), np.cos(rotate)]])
	# flip_x = -1.0 if True else 1.0
	# flip_y = 1.0
	# for i in range(3*num_triangles):
	# 	array_points[i, :] = np.array([flip_x, flip_y])*A.dot(array_points[i, :])

	return array_points, array_triangles, array_facecolors, norm_colors


def ReadFields_Uniform(file_name, return_vertex_coordinates=False):
	try:
		file = open(file_name, "rb")
	except FileNotFoundError:
		return None


	time = float(file.readline().decode("ascii").split()[1])
	nx, ny = np.array( file.readline().decode('ascii').split() ).astype(int)
	field_names = file.readline().decode("ascii").split(";")

	num_cells = nx*ny
	fields = {}
	for field_name in field_names:
		fields[field_name.strip()] = np.reshape( np.fromfile(file, np.float32, num_cells), (nx, ny) ).T

	fields["x-centers"] = np.copy(fields["x"])
	fields["y-centers"] = np.copy(fields["y"])

	if( return_vertex_coordinates ):
		mesh_x = fields["x"]
		dx = mesh_x[0, 1] - mesh_x[0, 0]
		array_x = np.linspace(mesh_x[0, 0] - 0.5*dx, mesh_x[0, -1] + 0.5*dx, nx+1)
		mesh_y = fields["y"]
		dy = mesh_y[1, 0] - mesh_y[0, 0]
		array_y = np.linspace(mesh_y[0, 0] - 0.5*dy, mesh_y[-1, 0] + 0.5*dy, ny+1)
		fields["x"], fields["y"] = np.meshgrid(array_x, array_y)

	return fields

def InterpolateFromVTK(file_name, chosen_field, xy_lims):
	try:
		file = open(file_name, 'rb')
	except FileNotFoundError:
		return None, None, None, None, None, None
	
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)

	num_points = int(line.split(" ")[1])
	

	num_cells = np.rint(num_points/4).astype(int)
	array_cell_center = np.zeros(shape=(num_cells, 2))

	min_dx = 1e+10
	
	for i in range(num_cells):
		p1 = ReadPoint(file)
		p2 = ReadPoint(file)
		p3 = ReadPoint(file)
		p4 = ReadPoint(file)
		cell_center = 0.25*(p1 + p2 + p3 + p4)
		# radius = np.linalg.norm(cell_center - p1)

		dx = 2.0*np.abs(cell_center[0] - p1[0])
		if( dx<min_dx ):
			min_dx = dx

		array_cell_center[i, :] = cell_center

	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	for i in range(num_cells):
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)
		read_stuff = ReadInt(file)

	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	for i in range(num_cells):
		read_stuff = ReadInt(file)

	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)

	# Starting to read printed field data
	line = ReadLine_Binary(file)
	while( (line is not None) and line.split(" ")[0]=="SCALARS" ):
		field_name = line.split(" ")[1]
		line = ReadLine_Binary(file)
		cell_values = np.zeros(shape=(num_cells,))
		for i in range(num_cells):
			read_stuff = ReadDouble(file)
			if( field_name==chosen_field ):
				cell_values[i] = read_stuff

		line = ReadLine_Binary(file)
		line = ReadLine_Binary(file)




	print("mindx: ", min_dx)
	# Constructing the RBF interpolator
	print("Constructing rbf....\n")
	rbf_interpolator = Rbf(array_cell_center[:, 0], array_cell_center[:, 1], np.zeros(shape=(num_cells,)), cell_values)


	print("Interpolating values....\n")
	list_grid = []
	x = xy_lims[0]
	while(x<xy_lims[1]):
		y = xy_lims[2]
		while(y<xy_lims[3]):
			list_grid.append( [x, y, 0.0] )
			y += min_dx
		x += min_dx

	list_grid = np.array(list_grid)
	value = rbf_interpolator(list_grid[:, 0], list_grid[:, 1], list_grid[:, 2])
	print(value)





	return array_cell_center, cell_values

def read_tracer_particles(filename):
	file = open(filename, "rt")
	line = file.readline()
	line = file.readline()
	line = file.readline()
	line = file.readline()
	line = file.readline()
	num_points = int( line.split(" ")[1] )
	points = np.zeros(shape=(num_points, 2))
	for i in range(num_points):
		line = file.readline()
		line = line.split(" ")
		points[i, 0] = float(line[0])
		points[i, 1] = float(line[1])
	file.close()

	return points

def read_polydata(filename : str, only_2D : bool = True, rotate_2D : float = 0.0, flip : str | None = None, 
				  color : str = "black", color_range = [None, None], colormap = None, number_copies = 1):
	"""
	This function reads a VTK file "Interface_XXXX.vtk" coming from one of our Basilisk simulations.

	Parameters
		filename: the name of the interface file
		2D_only: if this file comes from a 2D simulation, set this as true. If 3D set as false. The output of the function is different in each case.
		rotate_2D: only relevant if 2D_only=True. You can rotate the visualization by a certain angle. This parameter is that angle (in degrees).
		flip_y: generate a second data set which is the same as the original, but with the y-coordinates flipped. Cute for axissymetric simulations where only the top half is simulated
		color: color of the lines of the interface. Can be a matplotlib solid color or the name of a property contained in the VTK file. In this case we will color based on that property with a colormap
	"""

	file = open(filename, "rb")
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	time = line.split(" ")[-1]
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	
	# Rotation matrix, which will be used later if requested
	sin_angle = np.sin(np.pi*rotate_2D/(180.0))
	cos_angle = np.cos(np.pi*rotate_2D/(180.0))
	rotation_matrix = np.array([[cos_angle, -sin_angle], [sin_angle, cos_angle]])

	flip = [-1, 1, 1] if flip=="x" else [1, -1, 1] if flip=="y" else None

	# Loading all the points into a numpy [n x 3] array
	num_points = int( str(line).split(" ")[1] )
	points = np.zeros(shape=(num_points, 3))
	for i in range(num_points):
		points[i, :] = ReadPoint_Float(file)
	flipped_points = flip*points if flip else None

	# Reading how many lines (in 2D) or polygons we will have
	line = ReadLine_Binary(file)
	line = ReadLine_Binary(file)
	num_polys = int(str(line).split(" ")[1])

	# If 2D: Discarding the third coordinate and applying rotation
	if( only_2D ):
		points = np.dot(points[:, :2], rotation_matrix.T)
		flipped_points = np.dot(flipped_points[:, :2], rotation_matrix.T) if flip else None

	# Looping over each line (in 2D) or polygon (in 3D) and reading which vertices define it
	collection = []
	flipped_collection = []
	segments = []
	degenerate_indices = []
	for i in range(num_polys):
		num_vertices = ReadInt(file) # How many vertices in this line/polygon
		
		indices_points = [] # List of vertices in this line/polygon
		for i in range(num_vertices):
			indices_points.append( ReadInt(file) )

		# Discarding degenerate cases that have less than 2 vertices
		if( num_vertices<2 ):
			degenerate_indices.append( i )
			continue

		segments.append( indices_points )

		if( only_2D ):
			collection.append(np.array([points[indices_points[0], :2], points[indices_points[1], :2]]))
			if( flip ):
				flipped_collection.append(np.array([flipped_points[indices_points[0], :], flipped_points[indices_points[1], :]]))
            

	# Now we are going to read the scalar data associated to each polygon
	line = ReadLine_Binary(file) # New line
	line = ReadLine_Binary(file) # CELL_DATA X
	dict_values = {}
	line = ReadLine_Binary(file) # SCALARS
	while( (line is not None) and line.split(" ")[0]=="SCALARS" ):
		field_name = line.split(" ")[1] 
		line = ReadLine_Binary(file) # LOOKUP_TABLE
		dict_values[field_name] = np.zeros(shape=(len(segments), ))
		i_segment = 0
		for i in range(num_polys):
			if( not(i in degenerate_indices) ):
				dict_values[field_name][i_segment] = ReadFloat(file)
				i_segment += 1
		line = ReadLine_Binary(file)
		line = ReadLine_Binary(file)

	file.close()

	if( color in dict_values.keys() ):
		values = dict_values[color]
		min_values, max_values = np.min(values), np.max(values)
		norm_colors = mpl.colors.Normalize(vmin=color_range[0] if color_range[0] else min_values, vmax=color_range[1] if color_range[1] else max_values) 
		colormap = Colormap('bids:plasma') if colormap is None else colormap
		color = colormap( norm_colors( values ) )

	if( only_2D ):
		line_collection = [LineCollection(collection, colors=color) for i in range(number_copies)]
		flipped_line_collection = [LineCollection(flipped_collection, colors=color) for i in range(number_copies) ] if len(flipped_collection)>0 else None
		line_collection = line_collection if len(line_collection)>1 else line_collection[0]
		if( flipped_line_collection is not None ):
			flipped_line_collection = flipped_line_collection if len(flipped_line_collection)>1 else flipped_line_collection[0]

		return float(time), dict_values, points, segments, line_collection, flipped_line_collection
	
	print("ERROR. read_polydata: need to implement 3D version.")
	exit()