#!/usr/bin/python 

from itaps import iMesh, iBase
import argparse
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
mesh = iMesh.Mesh()

def get_surfaces(filename):

    #load the mesh file
    mesh.load(filename)
    
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
def surface_intersections( tris, axis, coord):

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
    
    surfs = get_surfaces(args.filename)

    for surf in surfs: 

        surf_tris = surf.getEntities(iBase.Type.all, iMesh.Topology.triangle)
        print "Retrieved " + str(len(surf_tris)) + " triangles from a surface set."

        


if __name__ == "__main__":
    main()
