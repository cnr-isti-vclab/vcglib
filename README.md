The **_Visualization and Computer Graphics Library_** (VCGlib) is an open source, portable, and templated library written in C++, with no external dependencies, for manipulation, processing, cleaning, and simplifying triangle meshes.

The library, composed by more than 100k lines of code, is released under the GPL license, and it is the base of most of the software tools of the [Visual Computing Lab](http://vcg.isti.cnr.it) of the Italian National Research Council Institute - ISTI, like [MeshLab](http://www.meshlab.net/), [Metro](http://vcg.isti.cnr.it/vcglib/metro.html), and many others.

The VCG library is tailored to mostly manage triangular meshes. The library is fairly large and offers many state-of-the-art capabilities for processing meshes, such as:

- high quality quadric-error edge-collapse based simplfication
- efficient spatial query structures (uniform grids, hashed grids, kdtree, etc)
- advanced smoothing and fairing algorithms
- computation of curvature
- optimization of texture coordinates
- Hausdorff distance computation
- geodesic paths
- mesh repairing capabilities
- isosurface extraction and advancing front meshing algorithms
- Poisson disk sampling and other tools to sample point distributions over meshes
- subdivision surfaces

## Notable Applications

A number of applications have been developed using VCGlib:

- MeshLab: the renowed open source mesh processing program
- Metro, the tool for measuring differences between meshes
- The first high quality out-of-core mesh simplifier that was used by the Stanford Digital Michelangelo project to process their huge 3D scanned models.

## Contacts

For any info about licensing (portions of) the library please contact us:
Paolo Cignoni (p.cignoni@isti.cnr.it) 
Visual Computing Lab of the Italian National Research Council - ISTI

In case of bugs please report them [here](https://github.com/cnr-isti-vclab/vcglib/issues).
