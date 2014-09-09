#!/usr/bin/python 

from itaps import iMesh, iBase
import argparse
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
import numpy as np
import matplotlib.cm as colormap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
from pylab import * 
from random import random

import time
mesh = iMesh.Mesh()


CCW = 1
CW = -1

########

def point_match(pnt1, pnt2):

    b = False
    x1, y1, z1 = round(pnt1[0],6),round(pnt1[1],6),round(pnt1[2],6)
    x2, y2, z2 = round(pnt2[0],6),round(pnt2[1],6),round(pnt2[2],6)
    if x1 == x2 and y1 == y2 and z1 == z2:
        b = True
    return b

def surface_intersections(tris, axis, coord):
    
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

    #create an array for the coordinates
    triangle_verts=np.array([],ndmin=3)
    triangle_verts.shape =(0,3,3)

    verts = mesh.getEntAdj(triangle,iBase.Type.vertex)
    vert1 = mesh.getVtxCoords(verts[0])
    vert2 = mesh.getVtxCoords(verts[1])
    vert3 = mesh.getVtxCoords(verts[2])

    temp = np.vstack((vert1,vert2,vert3))
    #insert the coordinates into the array
    triangle_verts = np.append(triangle_verts,[temp],axis = 0)

    #check for an intersection                                                                                                                                                   
    line = triangle_plane_intersect(axis,coord,triangle_verts)

    #if a line is returned, indicate that we have an intersection
    intersect = True if line.size is not 0 else False

    return intersect, line

def insert_pnt( line, coll):

    
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
                

#########


def get_surfaces():

    
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

def get_vols_by_group():
    
    mesh_sets=mesh.getEntSets()


    all_vols = get_all_volumes()

    group_vols=[]
    group_names=[]
    for s in mesh_sets:
        tags=mesh.getAllTags(s)
        for t in tags:
            
            if t.name == 'NAME':
                group_name = tag_to_script(t[s])
                group_names.append(group_name)
                #assume only volumes in groups for now
                vols = []
                for vol in all_vols:
                    if s.contains(vol): vols.append(vol)

                group_vols.append(vols)

    return group_vols, group_names
                

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


    print "There are " + str(len(vols)) + " volumes in this model."
    return vols


def get_vol_intersections(volume, intersect_dict):

    #get the surfaces for this volume
    surfs = volume.getChildren(0)

    intersections = []
    # return all intersections from the dictionary with
    # the surfaces of the volume as the key
    for surf in surfs:
        intersections += intersect_dict[surf]

    return intersections

# this function assumes that each plane-volume intersection will result 
# in some number of complete loops
def stitch(intersections):

    colls = []
    #first, check for complete loops
    i=0
    while i < len(intersections):
        intersection = intersections[i]
        # add if this is full loop
        if point_match(intersection[0],intersection[-1]) and len(intersection) != 2:
            colls.append(intersection)
            del intersections[i]
            i=0
        # if this is a segment with a near-zerio length, remove it
        elif point_match(intersection[0],intersection[-1]) and len(intersection) == 2:
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
        if point_match(colls[-1][0],colls[-1][-1]):
            colls.append(intersections.pop())

        else:
            i=0
            while i < len(intersections):
                intersection=intersections[i]
                if point_match(colls[-1][0],intersection[0]):
                    #reverse and attach to the front of colls[-1]
                    colls[-1]=np.append(intersection[::-1],colls[-1],axis=0)
                    del intersections[i]
                    i=0
                elif point_match(colls[-1][0], intersection[-1]):
                    #attach to the front of colls[-1]
                    colls[-1]=np.append(intersection,colls[-1],axis=0)
                    #                   print type(intersection)
                    del intersections[i]
                    i=0
                elif point_match(colls[-1][-1], intersection[0]):
                    #attach to the back of colls[-1]
                    colls[-1]=np.append(colls[-1],intersection,axis=0)
                    del intersections[i]
                    i=0
                elif point_match(colls[-1][-1], intersection[-1]):
                    #reverse and attach to the back of colls[-1]
                    colls[-1]=np.append(colls[-1],intersection[::-1],axis=0)
                    del intersections[i]
                    i=0
                #if no match is found, move to next intersection
                i+=1

    return colls

def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='filename', help='The path to the .h5m file')

    parser.add_argument('-g', action='store_true', dest='by_group', 
                        help='Plot intersections by groups using the same color for each group')
    
    parser.set_defaults(by_group=False)

    args = parser.parse_args()

    if not args.filename:
        raise Exception('h5m file path not specified!!. [-f] is not set')

        
    return args

# This function sets up coding for the path such that interior loops
# will not be filled for a volume cross-section
def return_coding(ob):
    # The codes will be all "LINETO" commands, except for "MOVETO"s at the
    # beginning of each subpath
    n = len(ob)
    codes = np.ones(n, dtype=Path.code_type) * Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def create_surface_intersections(surfs, axis, coord):

    intersection_dict={}
    for surf in surfs: 
        # get the surface's triangles
        surf_tris = surf.getEntities(iBase.Type.all, iMesh.Topology.triangle)
        #print "Retrieved " + str(len(surf_tris)) + " triangles from a surface set."
        # generate the surface intersections
        surf_intersections = surface_intersections(surf_tris, axis, coord )
        #add the surface's entry to the dictionary
        intersection_dict[surf] = surf_intersections
    return intersection_dict

def slice_faceted_model(filename, coord, axis, by_group=False):

    mesh.load(filename)

    #get all surfaces in the file
    surfs = get_surfaces()

    #get the intersections on each surface and store in dict with surfs as the keys
    intersection_dict = create_surface_intersections(surfs, axis, coord)


    #get all volumes in the model
    volumes = get_all_volumes()


    all_paths=[]
    group_names=[]
    if by_group:
        # if by group has been requested, sort volumes into their groups
        group_vols, group_names =get_vols_by_group()

        for vols in group_vols:
            #get the coords and codes for each volume
            all_coords, all_codes = get_volume_paths(vols, axis, intersection_dict)
            #when doing this by group, concat all the volumes in a group into one path
            if all_coords == []:
                continue
            else:
                group_path = np.concatenate(all_coords[:],axis=0)
                group_code = np.concatenate(all_codes[:],axis=0)
            all_paths.append([group_path,group_code])
            
    else:
        #if by group is not requested, return a path and code for each volume
        all_coords, all_codes  = get_volume_paths(volumes, axis, intersection_dict)

        all_paths=zip(all_coords,all_codes)

    return all_paths, group_names

def get_volume_paths(vols, axis, intersection_dict):

    all_coordinates=[]
    all_codes=[]
    orient_time=0
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
        loops = [np.delete(loop,axis,1) for loop in loops]


        if __name__ == "__main__":
            print "Re-orienting intersections..."
        start = time.clock()
        #orient the loops for proper plot fills
        loops = orient_loops(loops)
        orient_time += (time.clock()-start)

        #Reformat
        #rearrange coords into one long list and remove the coordinates for the slice
        coordinates = np.concatenate(loops[:],axis=0)

        #generate coding for the path that will allow for interior loops (see return_coding)
        codes=np.concatenate([return_coding(loop) for loop in loops])
        
        #add this volume's info to the global list
        all_coordinates.append(coordinates)
        all_codes.append(codes)
        coordinates=[]
        codes=[]

    if __name__ == "__main__":
        print "Took " + str(orient_time) + " seconds to reorient loops."
    return all_coordinates, all_codes


def orient_loops(loops):
    
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

    n = len(current_windings)
    assert(len(current_windings)==len(desired_windings))
    assert(len(current_windings)==len(loops))

    for i in range(n):
        
        if current_windings[i] != desired_windings[i]:
            loops[i]=loops[i][::-1]

    return loops


def find_winding(path):

    area = 0
    verts = path.vertices
    j = len(verts) -1
    for i in range(len(verts)-1):
        area += (verts[j,0]+verts[i,0]) * (verts[j,1]-verts[i,1])
        j=i

    return CW if area >= 0 else CCW

def get_windings(paths):

    windings=[]
    for path in paths:
        winding = find_winding(path)
        windings.append(winding)
    return windings

def gen_containment(paths):

    n = len(paths)
    mat = np.empty([n,n])
    for i in range(n):
        for j in range(n):
            mat[i,j] = 1 if paths[j].contains_path(paths[i]) else 0

    return mat
    
def get_fill_windings(M):

    a,b = M.shape
    assert(a==b)
    #use the sum of the rows of M to determine windings for path[i]
    windings=[]
    # 1 indicates CCW; -1 indicates CW
    for i in range(a):
        wind = CCW if sum(M[i])%2 == 0 or sum(M[i]) == 0 else CW
        windings.append(wind)

    return windings
        

def main():

    #parse arguments and load the file
    args = parsing()

    axis = 1
    coord = 0.0
    all_paths, group_names = slice_faceted_model(args.filename, coord, axis, args.by_group)

    if __name__ == "__main__":
        print "Plotting..."

    patches=[]
    for coord, code in all_paths:
        path = Path(coord,code)
        color = np.random.rand(3,1)
        patches.append(PathPatch(path, color=color, ec='black', lw=1))

        
    #create a new figure
    fig, ax = plt.subplots()

    #add the patches to the plot
    for patch in patches:
        ax.add_patch(patch)

    #show the plot!
    ax.autoscale_view()    
    plt.show()  
       
if __name__ == "__main__":
    start = time.clock()
    main()
    print "Took " + str((time.clock()-start)) + " seconds total."




