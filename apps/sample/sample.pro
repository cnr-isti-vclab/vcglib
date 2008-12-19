######################################################################
# Hand made pro. 
######################################################################

TEMPLATE      = subdirs
SUBDIRS       = trimesh_base  \
                trimesh_topology\
                trimesh_smooth \
                trimesh_refine \
                trimesh_clustering \
                trimesh_isosurface \
                trimesh_join \
                trimesh_optional \
                trimesh_intersection \
                trimesh_ball_pivoting \
                trimesh_hole \
                polygonmesh_base \
                aabb_binary_tree \
                edgemesh_grid \
		    trimesh_attribute

sources.files = *.pro
sources.path = .
INSTALLS += sources
