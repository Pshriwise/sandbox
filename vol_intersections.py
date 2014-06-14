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
    # get the triangles for each volume and return them
    for set in vol_sets:
        #get the surface sets of the volume
        surfs= set.getChildren(0) 

        #for each surface, get the triangles and add them to the dummy_set
        for surf in surfs:
            tris = surf.getEntities(iBase.Type.all, iMesh.Topology.triangle)
            dummy_set.add(tris)
        
        #get all of the triangles for this volume from the dummy_set
        vol_tris = dummy_set.getEntities(iBase.Type.all, iMesh.Topology.triangle)
        #empty out our dummy_set before starting the next volume
        dummy_set.remove(vol_tris)
        #add this array of triangles to the total list
        all_vol_tris.append(vol_tris)

    #kill our dummy set to it isn't hanging around when we leave this function
    mesh.destroyEntSet(dummy_set)
    
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
    while len(tris) is not 0:
       
        intersect, line = intersection(axis, coord, tris[0])
        #if we find an intersection, start looking for a poly collection
        if intersect:
            #print "Printing intersection"
            #print line
            checked_tris , pcoll = create_pcoll( tris[0], axis, coord, line)
            pcolls = pcoll
            tris = list(Set(tris).difference(Set(checked_tris)))
            line = []
        tris = tris[1:]
    return pcolls

def create_pcoll(intersected_tri, axis, coord, pnts):
    
    pcoll = pnts[0]
    checked_tris = [intersected_tri]
    #get the all triangles adjacent to this one
    adj_tris = mesh.getEnt2ndAdj(intersected_tri, 1, iBase.Type.face)

    while len(adj_tris) is not 0:
        #print len(adj_tris)
        intersect, line = intersection(axis, coord, adj_tris[0])
        #print "Zero ID: " + str(id(adj_tris[0]))
        #print "Checked Tris IDs: "
        #for tri in checked_tris: print id(tri)
        #print not checked(adj_tris[0],checked_tris)
        if intersect and not checked(adj_tris[0],checked_tris):
            #print "Found new line!"
            #now we'll insert the point into the poly collection
            pcoll = insert_pnt(line[0], pcoll)
            #print pcoll
            #print "Next points: "
            #line.shape
            #print pcoll.shape
            checked_tris.append(adj_tris[0])
            adj_tris = mesh.getEnt2ndAdj(adj_tris[0],1, iBase.Type.face)
        else:
            checked_tris.append(adj_tris[0])
            adj_tris=adj_tris[1:]
            #print "Removed tri"


    return checked_tris, pcoll


def checked( tri, checked_tris):
    b = False
    for checked_tri in checked_tris:
        if checked_tri == tri: b = True
    return b

def insert_pnt( line, coll):

    if line[0].all() == coll[0].all():
        coll = np.insert(coll, 0, line[1], 0)
    if line[0].all() == coll[-1].all():
        coll = np.insert(coll, -1, line[1], 0)
    if line[1].all() == coll[0].all():
        coll = np.insert(coll, 0, line[0], 0)
    if line[1].all() == coll[-1].all():
        coll = np.insert(coll, -1, line[0], 0)
    else:
        print "Warning: no points coincide in this tri, moving on"

    return coll

if __name__ == "__main__":
    volume_tris = get_vol_tris("cyl.h5m")

for vol in volume_tris:

    poly_collection=vol_plane_intersections(vol,0,0.1)
    print poly_collection
    if len(poly_collection) is not 0:
        plot(poly_collection[:,0],poly_collection[:,1])
        show()
