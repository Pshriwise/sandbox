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

mesh = iMesh.Mesh()


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

def get_volumes():

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


def slice_faceted_model(filename):

    mesh.load(filename)
    
    #get all surfaces in the file
    surfs = get_surfaces()    
    intersection_dict={}

    for surf in surfs: 
        # get the surface's triangles
        surf_tris = surf.getEntities(iBase.Type.all, iMesh.Topology.triangle)
        #print "Retrieved " + str(len(surf_tris)) + " triangles from a surface set."
        # generate the surface intersections
        surf_intersections = surface_intersections(surf_tris, 2, 0 )
        #add the surface's entry to the dictionary
        intersection_dict[surf] = surf_intersections

    #get all the volumes
    vols = get_volumes()

    patches=[]
    for vol in vols:
        
        #get the intersections for this volume based on its child surfaces
        intersects = get_vol_intersections(vol, intersection_dict)
        #if no intersections are returned, move on to the next volume
        if len(intersects) == 0: continue
        print "Retrieved "+str(len(intersects))+" intersections for this volume."
        #order the intersections into loops
        collections = stitch(intersects)
        print "Found "+str(len(collections))+" poly collections for this volume."
        
        #PLOTTING
        #rearrange coords into one long list and remove the coordinates for the slice
        all_coords = np.delete(np.concatenate(collections[:],axis=0),2,1)   
        #generate coding for the path that will allow for interior loops (see return_coding)
        all_codes=np.concatenate([return_coding(collection) for collection in collections])
        #create a patch
        path = Path(all_coords, all_codes) 
        #make a patch for this path
        patches.append(PathPatch(path))
        #add the path to the plot

    return patches

def main():

    #parse arguments and load the file
    args = parsing()

    patches = slice_faceted_model(args.filename)

    colors = 100*np.random.rand(len(patches))
    p = PatchCollection(patches, alpha=0.4)
    p.set_array(np.array(colors))

    #create a new figure
    fig, ax = plt.subplots()

    ax.add_collection(p)
    #show the plot!
    ax.autoscale_view()    
    plt.show()  
       
if __name__ == "__main__":
    main()




