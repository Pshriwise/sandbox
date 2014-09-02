#!/usr/bin/python 

from itaps import iMesh, iBase
import argparse
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
import numpy as np
from matplotlib.collections import PolyCollection
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


def checked( tri, checked_tris):
    b = False
    for checked_tri in checked_tris:
        if checked_tri == tri: b = True
    return b


def remove_checked_tris(tris,checked_tris):
    tris = list(tris)
    checked_tris = list(checked_tris)
    tris = list(set(tris)-set(checked_tris))
    return np.array(tris)
               


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

    
    #get all the surface entsets
    sets = mesh.getEntSets(1)
    
    surfs = []
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

    #get all the surface entsets
    sets = mesh.getEntSets(0)
    
    vols = []
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

    intersects = []
    for surf in surfs:
        intersects += intersect_dict[surf]

    return intersects

# this function assumes that each plane-volume intersection will result 
# in some number of complete loops
def stitch(intersections):

    colls = []
    #first, check for complete loops
    i=0
    while i < len(intersections):
        intersection = intersections[i]
        if point_match(intersection[0],intersection[-1]) and len(intersection) != 2:
            colls.append(intersection)
            del intersections[i]
            i=0
        elif point_match(intersection[0],intersection[-1]) and len(intersection) == 2:
            del intersections[i]
            i=0
        i+=1
        
    
    if 0 == len(intersections):
        return colls
    
    #add the last intersection (arbitrary starting point)
    colls.append(intersections.pop())
    counter = 0
    while len(intersections) != 0:
        if point_match(colls[-1][0],colls[-1][-1]):
            colls.append(intersections.pop())
        else:
            #            print type(intersections)
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
                i+=1
        counter+=1;
    return colls

def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', action='store', dest='filename', help='The path to the .h5m file')

    args = parser.parse_args()

    if not args.filename:
        raise Exception('h5m file path not specified!!. [-f] is not set')

        
    return args

def main():

    args = parsing()

    #load the mesh file
    mesh.load(args.filename)
    
    surfs = get_surfaces()
    
    intersection_dict={}
    for surf in surfs: 

        surf_tris = surf.getEntities(iBase.Type.all, iMesh.Topology.triangle)
        print "Retrieved " + str(len(surf_tris)) + " triangles from a surface set."

        surf_intersections = surface_intersections(surf_tris, 2, 0.2 )

        intersection_dict[surf] = surf_intersections

    #print intersection_dict 

    vols = get_volumes()

    colors = ['c','g','r','m','b']
    color = colors[0]

    fig, ax = plt.subplots()
            
    for vol in vols:
        
        color = colors.pop(0)
        colors.append(color)
        
        intersects = get_vol_intersections(vol, intersection_dict)
        print "Retrieved "+str(len(intersects))+" intersections for this volume."
        collections = stitch(intersects)
        print "Found "+str(len(collections))+" poly collections for this volume."
        
        for collection in collections:
            coll = PolyCollection([list(zip(collection[:,0],collection[:,1]))], alpha=0.4)
            ax.add_collection(coll)
    
    ax.autoscale_view()    
    plt.show()  
       
        #print collections 

if __name__ == "__main__":
    main()




