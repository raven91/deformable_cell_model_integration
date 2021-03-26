# deformable_cell_model_integration

Implementation of the deformable cell model that describes the life cycle of a number of interacting cells, 
each of which is modelled as a viscoelastic mesh. Its mathematical description is largely based 
on the work of [Van Liedekerke et al. 2020](https://link.springer.com/article/10.1007/s10237-019-01204-7).

Minimum C++ Version: C++17.

Dependencies: [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).

Visualization of the cell dynamics is implemented in [@deformable_cell_model_rendering](https://github.com/raven91/deformable_cell_model_rendering).

Initial cell mesh can be generated as in [@deformable_cell_model_mesh_generation](https://github.com/raven91/deformable_cell_model_mesh_generation).
