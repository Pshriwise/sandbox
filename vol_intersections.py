#!/usr/bin/python 

import h5py
import pytaps_get_tris as pytris
import numpy as np
from sets import Set
import time
from itaps import iMesh, iBase 
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from matplotlib import collections
import matplotlib.pyplot as pyplt
from mpl_toolkits.mplot3d import Axes3D
from pylab import * 

mesh = iMesh.Mesh()

def get_vol_tris(filename):

    #load the mesh file 
    mesh.load(filename)
    
    #get all of the entity sets
    sets = mesh.getEntSets(1)
    
    #filter out the voulme sets
    vol_sets=[]
    for set in sets:
        tags = mesh.getAllTags(set)
        for tag in tags:
            if tag.name == "CATEGORY":
                tag_handle = mesh.getTagHandle(tag.name)
                category_type = filter(lambda null: null != 0 , tag_handle[set])
                category_type = ''.join(chr(item) for item in category_type)
                if category_type == "Volume": vol_sets.append(set)

    #create a set to store triangles from multiple surface temporarily
    dummy_set = mesh.createEntSet(False)

    #start a list that will be a list of arrays
    #each array containing all the triangles for a surface
    all_vol_tris=[]
    vol_tris=[]
    # get the triangles for each volume and return them
    for set in vol_sets:
        #get the surface sets of the volume
        surfs= set.getChildren(0) 

        #for each surface, get the triangles and add them to the dummy_set
        for surf in surfs:
            tris = surf.getEntities(iBase.Type.all, iMesh.Topology.triangle)
            #dummy_set.add(tris)
            vol_tris.append(tris)
        all_vol_tris.append(vol_tris)
    
    return all_vol_tris


def intersection(axis,coord,triangle):

    triangle_verts=np.array([],ndmin=3)
    triangle_verts.shape =(0,3,3)

    verts = mesh.getEntAdj(triangle,iBase.Type.vertex)
    vert1 = mesh.getVtxCoords(verts[0])
    vert2 = mesh.getVtxCoords(verts[1])
    vert3 = mesh.getVtxCoords(verts[2])

    temp = np.vstack((vert1,vert2,vert3))
    triangle_verts = np.append(triangle_verts,[temp],axis = 0)

    #check for an intersection                                                                                                                                                   
    line = triangle_plane_intersect(axis,coord,triangle_verts)

    intersect = True if line.size is not 0 else False

    return intersect, line
    
def vol_plane_intersections(tris, axis, coord):

    pcolls = []
    print "Original length of tris: " + str(len(tris))
    while len(tris) is not 0:
       
        intersect, line = intersection(axis, coord, tris[0])
        #if we find an intersection, start looking for a poly collection
        if intersect:
            #print "Printing intersection"
            #print line
            checked_tris , pcoll = create_pcoll( tris[0], axis, coord, line)
            print "Length of checked tris is: " + str(len(checked_tris))
            pcolls.append(pcoll)
            print "Current length of tris is: " + str(len(tris))
            tris = remove_checked_tris(tris,checked_tris)
            print "After removing checked tris, tris length is: " + str(len(tris))
            line = []
        else:
            tris = tris[1:]
    return pcolls

def remove_checked_tris(tris,checked_tris):
    tris = list(tris)
    checked_tris = list(checked_tris)
    tris = list(set(tris)-set(checked_tris))
    return np.array(tris)
               
                
def create_pcoll(intersected_tri, axis, coord, pnts):
    
    pcoll = pnts[0]
    checked_tris = [intersected_tri]
    #get the all triangles adjacent to this one
    adj_tris = mesh.getEnt2ndAdj(intersected_tri, 1, iBase.Type.face)

    while len(adj_tris) is not 0:
        #print len(adj_tris)
        intersect, line = intersection(axis, coord, adj_tris[0])

        if intersect and not checked(adj_tris[0],checked_tris):

            #now we'll insert the point into the poly collection
            pcoll = insert_pnt(line[0], pcoll)

            # add this triangle to the checked_tris now that we've inserted the point
            checked_tris.append(adj_tris[0])
      
            # because an intersection was found on this tri, add its adj triangles to the stack
            new_tris = mesh.getEnt2ndAdj(adj_tris[0], 1, iBase.Type.face)
            adj_tris = np.concatenate((adj_tris,new_tris), axis=0)

        else:
            checked_tris.append(adj_tris[0])
            adj_tris=adj_tris[1:]




    return checked_tris, pcoll


def checked( tri, checked_tris):
    b = False
    for checked_tri in checked_tris:
        if checked_tri == tri: b = True
    return b

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

def point_match(pnt1, pnt2):

    b = False
    x1, y1, z1 = round(pnt1[0],6),round(pnt1[1],6),round(pnt1[2],6)
    x2, y2, z2 = round(pnt2[0],6),round(pnt2[1],6),round(pnt2[2],6)
    if x1 == x2 and y1 == y2 and z1 == z2:
        b = True
    return b

def main():
    volume_tris = get_vol_tris("cyl.h5m")

    for vol in volume_tris:

        for surface in vol:
            poly_collections=vol_plane_intersections(surface,0,0)
            print len(poly_collections)
            for coll in poly_collections:
                if len(coll) is not 0:
                    print "length of collection is: " + str(len(coll))
                    plot(coll[:,1],coll[:,2])
    show()

if __name__ == "__main__":
    main()
