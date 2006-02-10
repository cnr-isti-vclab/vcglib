######################################################################
# Hand made pro. 
######################################################################

TEMPLATE      = subdirs
SUBDIRS       = trimesh_base \
                trimesh_topology\
                trimesh_smooth \
                trimesh_refine \
                trimesh_isosurface \
                trimesh_join \
                trimesh_optional \
                trimesh_intersection \
                aabb_binary_tree
                

sources.files = *.pro
sources.path = .
INSTALLS += sources
