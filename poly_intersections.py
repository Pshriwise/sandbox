#!/usr/bin/python 

from itaps import iMesh, iBase
import argparse
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt


import time
mesh = iMesh.Mesh()


CCW = 1
CW = -1

########


def point_match(pnt1, pnt2):
    """
    Returns true or false for matching points after rounding each coordinate to 6 digits.
    """
    b = False
    x1, y1, z1 = round(pnt1[0],6),round(pnt1[1],6),round(pnt1[2],6)
    x2, y2, z2 = round(pnt2[0],6),round(pnt2[1],6),round(pnt2[2],6)
    if x1 == x2 and y1 == y2 and z1 == z2:
        b = True
    return b

def surface_intersections(tris, axis, coord):
    """
    Returns a list of ordered intersections for a set of triangles.
    """
    lines=[]
    #get all intersections for these surface triangles
    for tri in tris:
        intersect, line = intersection(axis, coord, tri)
        if intersect:
            lines.append(line[0])

    #time to start building intersections
    intersections=[]

    while len(lines) !=0:
        #arbitrarily start a new intersection
        current_intersection = lines[0]
        del lines[0]

        i=0
        while i < len(lines):
            line = lines[i]
            #insers points in the appropriate place if there is a match at the front or the back
            if point_match(line[0],current_intersection[0]):
                current_intersection = np.insert(current_intersection, 0, line[1], 0)
                del lines[i]
                i=0
            elif point_match(line[0],current_intersection[-1]):
                current_intersection = np.append(current_intersection, [line[1]], 0)
                del lines[i]
                i=0
            elif point_match(line[1], current_intersection[0]):
                current_intersection = np.insert(current_intersection, 0, line[0], 0)
                del lines[i]
                i=0
            elif point_match(line[1],current_intersection[-1]):
                current_intersection = np.append(current_intersection, [line[0]], axis=0)
                del lines[i]
                i=0
            else:
                i+=1
        intersections.append(current_intersection)
    
    return intersections


def intersection(axis, coord, triangle):
    """
    Returns the intersection (if one exists) between a *triangle*
    and the plane based on the *axis* and *coord*. 
    Also returns a boolean indicating whether or not an
    intersection was found.
    """
    #create an array for the coordinates
    triangle_verts=np.array([], ndmin=3)
    triangle_verts.shape =(0, 3, 3)

    verts = mesh.getEntAdj(triangle, iBase.Type.vertex)
    vert1 = mesh.getVtxCoords(verts[0])
    vert2 = mesh.getVtxCoords(verts[1])
    vert3 = mesh.getVtxCoords(verts[2])

    temp = np.vstack( (vert1, vert2, vert3) )
    #insert the coordinates into the array
    triangle_verts = np.append(triangle_verts, [temp], axis = 0)

    #check for an intersection                                                                                                                                                   
    line = triangle_plane_intersect(axis, coord, triangle_verts)

    #if a line is returned, indicate that we have an intersection
    intersect = True if line.size is not 0 else False

    return intersect, line

def insert_pnt(line, coll):
    """
    Adds a point on the appropriate side of the ordered line
    collection. Only the non-conincident point is added.
    """
    
    if point_match(line[0],coll[0]):
        coll = np.insert(coll, 0, line[1], 0)
    elif point_match(line[0],coll[-1]):
        coll = np.append(coll, [line[1]], 0)
    elif point_match(line[1], coll[0]):
        coll = np.insert(coll, 0, line[0], 0)
    elif point_match(line[1],coll[-1]):
        coll = np.append(coll, [line[0]], axis=0)
    else:
        print "Warning: no points coincide in this tri, moving on"

    return coll
                

def get_surfaces():
    """
    Gets all surfaces in the current mesh instance
    and returns them as a list.
    """
    
    #get all the entsets
    sets = mesh.getEntSets(1)
    
    surfs = []
    #filter out the surfaces using the CATEGORY tag
    for set in sets:
        tags = mesh.getAllTags(set)
        
        for tag in tags:
            if tag.name == "CATEGORY":
                tag_handle = mesh.getTagHandle(tag.name)
                category_type = filter(lambda null: null != 0 , tag_handle[set])
                category_type = ''.join(chr(item) for item in category_type)
                if category_type == "Surface": surfs.append(set)

    print "There are " + str(len(surfs)) + " surfaces in this model."

    return surfs

def get_vols_by_group(volumes):
    """
    All group entity sets and name tag values
    in the current mesh instance are retrieved.
    
    The input list of volumes are then sorted into
    their appropriate group (if one exists).
    
    The grouped volumes, and group names are
    then returned as lists.
    """

    # Get all ent sets
    mesh_sets=mesh.getEntSets()

    group_vols=[]
    group_names=[]
    
    #For each set get all its tags
    for s in mesh_sets:
        tags=mesh.getAllTags(s)
        for t in tags:
            #if the tag label is NAME, we've found a group set
            if t.name == 'NAME':
                #save the group name
                group_name = tag_to_script(t[s])
                group_names.append(group_name)
                #add any volume in volumes that is 
                #contained by the group set
                vols = []
                for vol in volumes:
                    if s.contains(vol): vols.append(vol)

                group_vols.append(vols)

    return group_vols, group_names
                

#Borrowed this from Moataz
#https://github.com/moatazharb/DAGMC/tree/develop/tools/dagmc_get_materials.py
"""
function to transform the tags into strings
tag : string of the tag to add to tag_list
tag_list : vector of tags in the problem
returns tag_list
"""
def tag_to_script(tag):
    a = []
    # since we have a byte type tag loop over the 32 elements
    for part in tag:
        # if the byte char code is non 0
        if (part != 0):
            # convert to ascii
            a.append(str(unichr(part)))
            # join to end string
            test = ''.join(a)
            # the the string we are testing for is not in the list of found
            # tag values, add to the list of tag_values
    # if not already in list append to list
    #    if not any(test in s for s in tag_list):
    # the original code was incorrectly missing groups when one of the same
    # name with/without rho was added
    return test

    
def get_all_volumes():
    """
    Get all volumes in the current mesh instance
    and return them in a list.
    """
    #get all the entsets
    sets = mesh.getEntSets(0)
    
    vols = []
    #filter out all the volumes using the CATEGORY tag
    for set in sets:
        tags = mesh.getAllTags(set)
        
        for tag in tags:
            if tag.name == "CATEGORY":
                tag_handle = mesh.getTagHandle(tag.name)
                category_type = filter(lambda null: null != 0 , tag_handle[set])
                category_type = ''.join(chr(item) for item in category_type)
                if category_type == "Volume": vols.append(set)


    #print "There are " + str(len(vols)) + " volumes in this model."
    return vols


def get_vol_intersections(volume, intersect_dict):
    """
    Get's all child surface handles of the *volume*
    and returnes all ordered intersections from the
    input *intersect_dict*.
    """
    #get the surfaces for this volume
    surfs = volume.getChildren(0)

    intersections = []
    # return all intersections from the dictionary with
    # the surfaces of the volume as the key
    for surf in surfs:
        intersections += intersect_dict[surf]

    return intersections

def stitch(intersections):
    """
    Takes an input list of ordered *intersections* and
    arranges them into complete loops. The loops are then
    returned as a list.

    NOTE: THIS FUNCTION ASSUMES THAT EACH VOLUME-PLANE INTERSECTION
    WILL RESULT IN SOME NUMBER OF FINITE, COMPLETE LOOPS. IF THE 
    SLICED MODEL IS NOT CONSIDERED WATERTIGHT, THIS FUNCTION MAY
    NEVER FINISH.
    """
    colls = []
    #first, check for complete loops
    i=0
    while i < len(intersections):
        intersection = intersections[i]
        # add if this is full loop
        if point_match(intersection[0], intersection[-1]) and len(intersection) != 2:
            colls.append(intersection)
            del intersections[i]
            i=0
        # if this is a segment with a near-zerio length, remove it
        elif point_match(intersection[0], intersection[-1]) and len(intersection) == 2:
            del intersections[i]
            i=0
        i+=1
        
    if 0 == len(intersections):
        return colls
    
    #add the last intersection (arbitrary starting point)
    colls.append(intersections.pop())
    
    #do until all intersections are matched
    while len(intersections) != 0:
        
        #if the current collection is a loop, move to 
        #the next arb starting point
        if point_match(colls[-1][0], colls[-1][-1]):
            colls.append(intersections.pop())

        else:
            i = 0
            while i < len(intersections):
                intersection = intersections[i]
                if point_match(colls[-1][0], intersection[0]):
                    #reverse and attach to the front of colls[-1]
                    colls[-1] = np.append(intersection[::-1], colls[-1], axis=0)
                    del intersections[i]
                    i = 0
                elif point_match(colls[-1][0], intersection[-1]):
                    #attach to the front of colls[-1]
                    colls[-1] = np.append(intersection, colls[-1], axis=0)
                    #                   print type(intersection)
                    del intersections[i]
                    i = 0
                elif point_match(colls[-1][-1], intersection[0]):
                    #attach to the back of colls[-1]
                    colls[-1] = np.append(colls[-1], intersection, axis=0)
                    del intersections[i]
                    i = 0
                elif point_match(colls[-1][-1], intersection[-1]):
                    #reverse and attach to the back of colls[-1]
                    colls[-1] = np.append(colls[-1], intersection[::-1], axis=0)
                    del intersections[i]
                    i = 0
                #if no match is found, move to next intersection
                i+=1

    return colls

def parsing():
    """
    Sets up expected arguments when running the file as main. 
    """

    #Set help description
    description = 'This program is designed for the plotting of a watertight, DAGMC-ready .h5m file' 
    description += 'w/ triangular facets.'
    parser = argparse.ArgumentParser(description=description)
    
    #Filename argument
    parser.add_argument(
        '-f', action='store', dest='filename', required=True, help='The path to the .h5m file')

    # Indicates the axis along which to slice
    parser.add_argument('-axis', action='store', dest='axis', type=int,
                        help='Set axis along which to slice the model x=0, y=1, z=2 (default = x)')

    # Defines the coordinate of the plane along which to slice. 
    parser.add_argument('-coord', action='store', dest='coord', type=float,
                        help='Coordinate for the slice (default = 0)')
    
    # Indicates that the user wishes to plot groups as the same color
    parser.add_argument('--by-group', action='store_true', dest='by_group', 
                        help='Plot intersections by groups using the same color for each group')

    # Option for writing raw point data to file 
    parser.add_argument( '--write-pnts', action='store_true', dest='write_pnts', help = 'If set, the program will now write raw point data to file name "slicepnts.txt".')

    parser.set_defaults(axis = 0)
    parser.set_defaults(coord = 0)
    parser.set_defaults(by_group = False)
    parser.set_defaults(file_out = 'slicepnts.txt')
    parser.set_defaults(write_pnts = False)
    args = parser.parse_args()

    if not args.filename:
        raise Exception('h5m file path not specified!!. [-f] is not set')

    return args

def return_coding(ob):
    """
    The codes will be all "LINETO" commands, except for "MOVETO"s at the
    beginning of each subpath.
    """
    n = len(ob)
    codes = np.ones(n, dtype=Path.code_type) * Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def create_surface_intersections(surfs, axis, coord):
    """
    Generates and returns a dictionary with *surfs* as the keys
    and the lists ordered intersections returned from *surface_intersections*
    as the values.
    """
    intersection_dict={}
    for surf in surfs: 
        # get the surface's triangles
        surf_tris = surf.getEntities(iBase.Type.all, iMesh.Topology.triangle)
        #print "Retrieved " + str(len(surf_tris)) + " triangles from a surface set."
        # generate the surface intersections
        surf_intersections = surface_intersections(surf_tris, axis, coord)
        #add the surface's entry to the dictionary
        intersection_dict[surf] = surf_intersections
    return intersection_dict

def slice_faceted_model(filename, coord, axis, by_group=False):
    """
    Returns a list of the paths required for patches and group_names (if called for)
    for slice at *axis* and *coord* through a model contained by *filename*. 

    If *by_group* is true, the path information returned will be sets of volume
    patches sorted into their appropriate groups.
    """
    mesh.load(filename)

    #get all surfaces in the file
    surfs = get_surfaces()

    #get the intersections on each surface and store in dict with surfs as the keys
    intersection_dict = create_surface_intersections(surfs, axis, coord)


    #get all volumes in the model
    volumes = get_all_volumes()


    all_paths = []
    group_names = []
    if by_group:
        # if by group has been requested, sort volumes into their groups
        group_vols, group_names = get_vols_by_group(volumes)

        for vols, name in zip(group_vols, group_names):
            #get the coords and codes for each volume
            all_coords, all_codes = get_volume_paths(vols, axis, intersection_dict)
            #when doing this by group, concat all the volumes in a group into one path
            if all_coords == []:
                group_names.remove(name)
                continue
            else:
                group_path = np.concatenate(all_coords[:], axis=0)
                group_code = np.concatenate(all_codes[:], axis=0)
            all_paths.append([group_path, group_code])
            
    else:
        #if by group is not requested, return a path and code for each volume
        all_coords, all_codes  = get_volume_paths(volumes, axis, intersection_dict)

        all_paths=zip(all_coords, all_codes)

    return all_paths, group_names

def get_volume_paths(vols, axis, intersection_dict):
    """
    For a given set of volumes, this will return a list of 2D paths
    for each volume in *vols*.
    """
    all_coordinates = []
    all_codes = []
    orient_time = 0
    for vol in vols:
        
        #get the intersections for this volume based on its child surfaces
        intersects = get_vol_intersections(vol, intersection_dict)

        #if no intersections are returned, move on to the next volume
        if len(intersects) == 0: continue
        print "Retrieved "+str(len(intersects))+" intersections for this volume."

        #order the intersections into loops
        loops = stitch(intersects)
        print "Found "+str(len(loops))+" poly collections for this volume."

        #remove the axis of intersection from the points
        loops = [np.delete(loop, axis, 1) for loop in loops]


        if __name__ == "__main__":
            print "Re-orienting intersections..."
        start = time.clock()
        #orient the loops for proper plot fills
        loops = orient_loops(loops)
        orient_time += (time.clock()-start)

        #Reformat
        #rearrange coords into one long list and remove the coordinates for the slice
        coordinates = np.concatenate(loops[:], axis=0)

        #generate coding for the path that will allow for interior loops (see return_coding)
        codes = np.concatenate([return_coding(loop) for loop in loops])
        
        #add this volume's info to the global list
        all_coordinates.append(coordinates)
        all_codes.append(codes)
        coordinates = []
        codes = []

    if __name__ == "__main__":
        print "Took " + str(orient_time) + " seconds to reorient loops."
    return all_coordinates, all_codes


def orient_loops(loops):
    """
    Takes a set of loops generated by the stitch function and orients
    the windings of the loops such that the correct set of paths will
    be filled in a plot according to the non-zero fill rule. 
    """
    #first we will create a set of paths to represent the loops
    paths = []
    for loop in loops:
        paths.append(Path(loop))

    #now we will generate a containment matrix for the loops
    M = gen_containment(paths)
    #determine the current windings
    current_windings = get_windings(paths)
    #get the desired windings
    desired_windings = get_fill_windings(M)
    #alter the current windings to match the desired windings
    loops = set_windings(current_windings, desired_windings, loops)

    return loops

def set_windings(current_windings, desired_windings, loops):
    """
    Compares the *current_windings* to the *desired_winding(
    for each loop in *loops*.
    If the winding values do not match, the points are reversed.

    A list of the same size as *loops* is returned. 
    """
    n = len(current_windings)
    assert(len(current_windings) == len(desired_windings))
    assert(len(current_windings) == len(loops))

    for i in range(n):
        if current_windings[i] != desired_windings[i]:
            #reverse the loop
            loops[i] = loops[i][::-1]

    return loops

def find_winding(path):
    """
    Uses an by product of calculating a polygon's area
    using a simplified Green's Theorem. Using this 
    method, if the returned area is negative, the orientation
    of the vertices if CCW. If positive, CW.
    """
    area = 0
    verts = path.vertices
    j = len(verts) -1
    for i in range(len(verts)-1):
        area += (verts[j,0] + verts[i,0]) * (verts[j,1] - verts[i,1])
        j=i

    return CW if area >= 0 else CCW

def get_windings(paths):
    """
    Gets the windings of each path in path using find_winding()
    and return them in a list. 
    """
    windings=[]
    for path in paths:
        winding = find_winding(path)
        windings.append(winding)

    return windings

def gen_containment(paths):
    """
    Generates a binary containment matrix in which the 
    i,jth entry indicates whether or not path i is 
    contained by path j. 

    Uses the matplotlib path method .contains_path() to
    set the value.

    NOTE: NO PATHS SHOULD CROSS EACH OTHER, SO THIS WILL
    ALWAYS GIVE A DEFINITE NESTING OF THE PATHS.
    """
    n = len(paths)
    mat = np.empty([n, n])
    for i in range(n):
        for j in range(n):
            mat[i,j] = 1 if paths[j].contains_path(paths[i]) else 0

    return mat
    
def get_fill_windings(M):
    """
    Uses the sum of the rows in the binary containment matrix *M* to
    calculate the desired windings of the paths represented by the matrix.
    """
    a,b = M.shape
    assert(a == b)
    #use the sum of the rows of M to determine windings for path[i]
    windings = []
    # 1 indicates CCW; -1 indicates CW
    for i in range(a):
        wind = CCW if sum(M[i])%2 == 0 or sum(M[i]) == 0 else CW
        windings.append(wind)

    return windings
        
def main():

    #parse arguments and load the file
    args = parsing()

    print args.coord
    print args.axis
    all_paths, group_names = slice_faceted_model(args.filename, args.coord, args.axis, args.by_group)

    if __name__ == "__main__":
        print "Plotting..."

    file = open('slicepnts.txt', 'a')
    
    patches = []
    for coord, code in all_paths:
        if args.write_pnts: np.savetxt(file, coord, delimiter = ' ')
        path = Path(coord, code)
        color = np.random.rand(3, 1)
        patches.append(PathPatch(path, color=color, ec='black', lw=1, alpha=0.4))

        
    #create a new figure
    fig, ax = plt.subplots()

    #add the patches to the plot
    for patch in patches:
        ax.add_patch(patch)

    if args.by_group:
        ax.legend(patches, group_names, prop={'size':10})
    #show the plot!
    ax.autoscale_view()
    ax.set_aspect('equal')
    plt.show()  
   
    
if __name__ == "__main__":
    start = time.clock()
    main()
    print "Took " + str((time.clock()-start)) + " seconds total."




