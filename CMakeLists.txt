# Copyright 2019, 2020, Collabora, Ltd.
# Copyright 2019, 2021, Visual Computing Lab, ISTI - Italian National Research Council
# SPDX-License-Identifier: BSL-1.0

cmake_minimum_required(VERSION 3.10)
project(VCGLib)

# Eigen options
option(ALLOW_BUNDLED_EIGEN "Allow use of bundled Eigen source" ON)
option(ALLOW_SYSTEM_EIGEN "Allow use of system-provided Eigen" ON)

# VCG options
option(VCG_HEADER_ONLY "Use VCG library in header only mode" ON)
option(VCG_BUILD_EXAMPLES "Build a set of examples of the library" OFF)

set (VCG_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR})
set (VCG_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)

### Build settings
set(CMAKE_CXX_STANDARD 11)

### Eigen
set(VCG_EIGEN_DIR ${CMAKE_CURRENT_LIST_DIR}/eigenlib)

if(ALLOW_SYSTEM_EIGEN AND EIGEN3_INCLUDE_DIR)
	message(STATUS "- Eigen - using system-provided library")
	set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
elseif(ALLOW_BUNDLED_EIGEN AND EXISTS "${VCG_EIGEN_DIR}/Eigen/Eigen")
	message(STATUS "- Eigen - using bundled source")
	set(EIGEN_INCLUDE_DIRS ${VCG_EIGEN_DIR})
else()
	message(
			FATAL_ERROR
			"Eigen is required - at least one of ALLOW_SYSTEM_EIGEN or ALLOW_BUNDLED_EIGEN must be enabled and found.")
endif()

### VCGLib headers and sources

set(VCG_HEADERS
		vcg/complex/append.h
		vcg/complex/all_types.h
		vcg/complex/complex.h
		vcg/complex/allocate.h
		vcg/complex/exception.h
		vcg/complex/algorithms/overlap_estimation.h
		vcg/complex/algorithms/dual_meshing.h
		vcg/complex/algorithms/intersection.h
		vcg/complex/algorithms/clip.h
		vcg/complex/algorithms/geodesic.h
		vcg/complex/algorithms/parametrization/poisson_solver.h
		vcg/complex/algorithms/parametrization/uv_utils.h
		vcg/complex/algorithms/parametrization/distortion.h
		vcg/complex/algorithms/parametrization/tangent_field_operators.h
		vcg/complex/algorithms/parametrization/voronoi_atlas.h
		vcg/complex/algorithms/edge_collapse.h
		vcg/complex/algorithms/hole.h
		vcg/complex/algorithms/align_pair.h
		vcg/complex/algorithms/closest.h
		vcg/complex/algorithms/tetra_implicit_smooth.h
		vcg/complex/algorithms/bitquad_support.h
		vcg/complex/algorithms/skeleton.h
		vcg/complex/algorithms/symmetry.h
		vcg/complex/algorithms/voronoi_volume_sampling.h
		vcg/complex/algorithms/polygon_polychord_collapse.h
		vcg/complex/algorithms/inside.h
		vcg/complex/algorithms/local_optimization/tri_edge_flip.h
		vcg/complex/algorithms/local_optimization/quad_diag_collapse.h
		vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h
		vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric_tex.h
		vcg/complex/algorithms/local_optimization/tri_edge_collapse.h
		vcg/complex/algorithms/local_optimization/tetra_edge_collapse.h
		vcg/complex/algorithms/polygonal_algorithms.h
		vcg/complex/algorithms/inertia.h
		vcg/complex/algorithms/mesh_assert.h
		vcg/complex/algorithms/occupancy_grid.h
		vcg/complex/algorithms/meshtree.h
		vcg/complex/algorithms/align_global.h
		vcg/complex/algorithms/cut_tree.h
		vcg/complex/algorithms/nring.h
		vcg/complex/algorithms/tetra/tetfuse_collapse.h
		vcg/complex/algorithms/stat.h
		vcg/complex/algorithms/ransac_matching.h
		vcg/complex/algorithms/refine.h
		vcg/complex/algorithms/outline_support.h
		vcg/complex/algorithms/convex_hull.h
		vcg/complex/algorithms/clean.h
		vcg/complex/algorithms/mesh_to_matrix.h
		vcg/complex/algorithms/quadrangulator.h
		vcg/complex/algorithms/isotropic_remeshing.h
		vcg/complex/algorithms/smooth.h
		vcg/complex/algorithms/autoalign_4pcs.h
		vcg/complex/algorithms/local_optimization.h
		vcg/complex/algorithms/curve_on_manifold.h
		vcg/complex/algorithms/clustering.h
		vcg/complex/algorithms/refine_loop.h
		vcg/complex/algorithms/cylinder_clipping.h
		vcg/complex/algorithms/pointcloud_normal.h
		vcg/complex/algorithms/bitquad_creation.h
		vcg/complex/algorithms/crease_cut.h
		vcg/complex/algorithms/implicit_smooth.h
		vcg/complex/algorithms/voronoi_remesher.h
		vcg/complex/algorithms/polygon_support.h
		vcg/complex/algorithms/point_sampling.h
		vcg/complex/algorithms/create/mc_lookup_table.h
		vcg/complex/algorithms/create/mc_trivial_walker.h
		vcg/complex/algorithms/create/extrude.h
		vcg/complex/algorithms/create/resampler.h
		vcg/complex/algorithms/create/ball_pivoting.h
		vcg/complex/algorithms/create/readme.txt
		vcg/complex/algorithms/create/zonohedron.h
		vcg/complex/algorithms/create/platonic.h
		vcg/complex/algorithms/create/marching_cubes.h
		vcg/complex/algorithms/create/plymc/voxel.h
		vcg/complex/algorithms/create/plymc/simplemeshprovider.h
		vcg/complex/algorithms/create/plymc/tri_edge_collapse_mc.h
		vcg/complex/algorithms/create/plymc/volume.h
		vcg/complex/algorithms/create/plymc/plymc.h
		vcg/complex/algorithms/create/plymc/svoxel.h
		vcg/complex/algorithms/create/tetramesh_support.h
		vcg/complex/algorithms/create/advancing_front.h
		vcg/complex/algorithms/textcoord_optimization.h
		vcg/complex/algorithms/bitquad_optimization.h
		vcg/complex/algorithms/halfedge_quad_clean.h
		vcg/complex/algorithms/voronoi_processing.h
		vcg/complex/algorithms/update/quality.h
		vcg/complex/algorithms/update/selection.h
		vcg/complex/algorithms/update/fitmaps.h
		vcg/complex/algorithms/update/component_ep.h
		vcg/complex/algorithms/update/texture.h
		vcg/complex/algorithms/update/curvature_fitting.h
		vcg/complex/algorithms/update/normal.h
		vcg/complex/algorithms/update/position.h
		vcg/complex/algorithms/update/halfedge_topology.h
		vcg/complex/algorithms/update/topology.h
		vcg/complex/algorithms/update/flag.h
		vcg/complex/algorithms/update/bounding.h
		vcg/complex/algorithms/update/halfedge_indexed.h
		vcg/complex/algorithms/update/color.h
		vcg/complex/algorithms/update/curvature.h
		vcg/complex/algorithms/point_outlier.h
		vcg/complex/algorithms/harmonic.h
		vcg/complex/algorithms/point_matching_scale.h
		vcg/complex/algorithms/attribute_seam.h
		vcg/complex/foreach.h
		vcg/complex/base.h
		vcg/complex/used_types.h
		vcg/container/entries_allocation_table.h
		vcg/container/container_allocation_table.h
		vcg/container/derivation_chain.h
		vcg/container/vector_occ.h
		vcg/container/simple_temporary_data.h
		vcg/space/segment2.h
		vcg/space/fitting3.h
		vcg/space/tetra3.h
		vcg/space/triangle2.h
		vcg/space/ray2.h
		vcg/space/deprecated_point2.h
		vcg/space/point4.h
		vcg/space/box2.h
		vcg/space/ray3.h
		vcg/space/planar_polygon_tessellation.h
		vcg/space/texcoord2.h
		vcg/space/deprecated_point3.h
		vcg/space/intersection/triangle_triangle3.h
		vcg/space/distance2.h
		vcg/space/point3.h
		vcg/space/deprecated_point.h
		vcg/space/space.h
		vcg/space/point.h
		vcg/space/colorspace.h
		vcg/space/rect_packer.h
		vcg/space/triangle3.h
		vcg/space/obox3.h
		vcg/space/point2.h
		vcg/space/smallest_enclosing.h
		vcg/space/color4.h
		vcg/space/polygon3.h
		vcg/space/line3.h
		vcg/space/index/octree.h
		vcg/space/index/grid_util2d.h
		vcg/space/index/grid_closest.h
		vcg/space/index/grid_static_ptr.h
		vcg/space/index/grid_util.h
		vcg/space/index/spatial_hashing.h
		vcg/space/index/closest2d.h
		vcg/space/index/grid_static_obj.h
		vcg/space/index/kdtree/kdtree.h
		vcg/space/index/kdtree/priorityqueue.h
		vcg/space/index/kdtree/kdtree_face.h
		vcg/space/index/kdtree/mlsutils.h
		vcg/space/index/octree_template.h
		vcg/space/index/aabb_binary_tree/kclosest.h
		vcg/space/index/aabb_binary_tree/closest.h
		vcg/space/index/aabb_binary_tree/ray.h
		vcg/space/index/aabb_binary_tree/frustum_cull.h
		vcg/space/index/aabb_binary_tree/aabb_binary_tree.h
		vcg/space/index/aabb_binary_tree/base.h
		vcg/space/index/grid_closest2d.h
		vcg/space/index/spatial_hashing2d.h
		vcg/space/index/space_iterators.h
		vcg/space/index/grid_static_ptr2d.h
		vcg/space/index/base2d.h
		vcg/space/index/base.h
		vcg/space/index/perfect_spatial_hashing.h
		vcg/space/index/space_iterators2d.h
		vcg/space/line2.h
		vcg/space/point_matching.h
		vcg/space/intersection3.h
		vcg/space/deprecated_point4.h
		vcg/space/rasterized_outline2_packer.h
		vcg/space/box.h
		vcg/space/plane3.h
		vcg/space/outline2_packer.h
		vcg/space/segment3.h
		vcg/space/intersection2.h
		vcg/space/sphere3.h
		vcg/space/box3.h
		vcg/space/distance3.h
		vcg/math/quadric5.h
		vcg/math/factorial.h
		vcg/math/eigen_matrix_addons.h
		vcg/math/quadric.h
		vcg/math/perlin_noise.h
		vcg/math/shot.h
		vcg/math/spherical_harmonics.h
		vcg/math/eigen_matrixbase_addons.h
		vcg/math/quaternion.h
		vcg/math/similarity.h
		vcg/math/disjoint_set.h
		vcg/math/random_generator.h
		vcg/math/camera.h
		vcg/math/linear.h
		vcg/math/matrix44.h
		vcg/math/eigen.h
		vcg/math/similarity2.h
		vcg/math/gen_normal.h
		vcg/math/polar_decomposition.h
		vcg/math/base.h
		vcg/math/histogram.h
		vcg/math/legendre.h
		vcg/math/matrix33.h
		vcg/simplex/edge/distance.h
		vcg/simplex/edge/topology.h
		vcg/simplex/edge/pos.h
		vcg/simplex/edge/component.h
		vcg/simplex/edge/base.h
		vcg/simplex/tetrahedron/tetrahedron.h
		vcg/simplex/tetrahedron/topology.h
		vcg/simplex/tetrahedron/pos.h
		vcg/simplex/tetrahedron/component.h
		vcg/simplex/tetrahedron/base.h
		vcg/simplex/face/component_occ.h
		vcg/simplex/face/component_ep.h
		vcg/simplex/face/jumping_pos.h
		vcg/simplex/face/distance.h
		vcg/simplex/face/component_polygon.h
		vcg/simplex/face/topology.h
		vcg/simplex/face/pos.h
		vcg/simplex/face/component.h
		vcg/simplex/face/component_ocf.h
		vcg/simplex/face/base.h
		vcg/simplex/vertex/component_occ.h
		vcg/simplex/vertex/component_sph.h
		vcg/simplex/vertex/distance.h
		vcg/simplex/vertex/component.h
		vcg/simplex/vertex/component_ocf.h
		vcg/simplex/vertex/base.h
		vcg/connectors/halfedge_pos.h
		vcg/connectors/hedge.h
		vcg/connectors/hedge_component.h

		#wrap
		wrap/callback.h
		)

set(SOURCES
		)

if (VCG_HEADER_ONLY)
	if (NOT TARGET vcglib) # to be sure that vcglib target is created just one time
		add_library(vcglib INTERFACE)
		target_include_directories(
				vcglib INTERFACE
				${CMAKE_CURRENT_LIST_DIR}
				${EIGEN_INCLUDE_DIRS})

		#just to show headers in ide
		add_custom_target(vcglib_ide SOURCES ${VCG_HEADERS})
	else()
		message(STATUS "- VCGLib - jumped - already included")
	endif()
else()
	#TODO make vcglib that includes all the wrap sources, checking everytime
	# if the the required targets (e.g. qt, gl, glew...) exists
endif()

if(VCG_BUILD_EXAMPLES)
	#TODO make the list of samples to build
	add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/apps)
endif()
