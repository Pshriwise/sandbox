# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from itaps import iBase, iMesh
import numpy as np
mesh = iMesh.Mesh()

# <codecell>

def get_vol_ent_sets(sets):
    
    vol_sets = []
    # Filter out sets that aren't volumes using the category tag
    for set in sets:
        tags = mesh.getAllTags(set)
        for tag in tags:
            if tag.name == "CATEGORY":
                tag_handle = mesh.getTagHandle(tag.name)
                category_type = filter(lambda null: null != 0 , tag_handle[set])
                category_type = ''.join(chr(item) for item in category_type)
                if category_type == "Volume": vol_sets.append(set)
    return vol_sets
    


# <codecell>

def get_tri_points_for_vols(vol_sets):
    #Get the triangle points for a given volume
    all_vol_tris=[]
    for vol in vol_sets:
        all_vol_tris.append(get_tri_points_of_vol(vol))
        #append to list of tri_stacks
    return all_vol_tris

# <codecell>


def get_tri_points_of_vol(vol_set):
    # Get all the surfaces of the first volume
    surfs = vol_set.getChildren(0)
    vol_tri_stack = np.array([],ndmin=3)
    vol_tri_stack.shape =(0,3,3)
    for surf in surfs:
        surf_tris = get_surf_tri_points(surf)   
        vol_tri_stack = np.vstack((vol_tri_stack,get_surf_tri_points(surf)))
        #append to volume tri stack
    # Return a list of triangles of the volume
    return vol_tri_stack
    

# <codecell>

#Get the triangle points for a given surface
def get_surf_tri_points(surf_set):
    surf_tris = surf_set.getEntities(iBase.Type.all,iMesh.Topology.triangle)
    
    # Create a numpy array to hold the triangle points in the correct format
    surf_tri_stack=np.array([],ndmin=3)
    surf_tri_stack.shape = (0,3,3)
    
    # Add all vertex coordinates to a numpy array
    for tri in surf_tris:
        verts = mesh.getEntAdj(tri,iBase.Type.vertex)    
    
        vert1 = mesh.getVtxCoords(verts[0])
    
        vert2 = mesh.getVtxCoords(verts[1])
    
        vert3 = mesh.getVtxCoords(verts[2])
    
        triangle_verts = np.vstack((vert1,vert2,vert3))
        
        surf_tri_stack = np.append(surf_tri_stack, [triangle_verts], axis =0)
        
    return surf_tri_stack
    



    

# <codecell>

def get_tri_points_by_vol(filename):
    # Load the mesh file into the iMesh instance
    mesh.load(filename)
    # Get all of the entity sets
    sets = mesh.getEntSets(iBase.Type.all)
    vol_sets = get_vol_ent_sets(sets)
    vol_tri_points = get_tri_points_for_vols(vol_sets)

    return vol_tri_points
    

# <codecell>

#mesh.load("/home/patrick/cube.h5m")
#mesh.load("/home/patrick/scratch/moab_tools/make_watertight/test/iter_imprinted.h5m")
    



