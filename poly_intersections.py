#!/usr/bin/python 

from itaps import iMesh, iBase
import argparse

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


if __name__ == "__main__":
    main()
