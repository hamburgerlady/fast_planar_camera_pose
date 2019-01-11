# fast_planar_camera_pose
A solver for camera pose estimation with unknown radial distortion and unknown focal length from four coplanar points as described in "A fast minimal solver for absolute camera pose with unknown focal length and radial distortion from four planar points", M. Oskarsson, arXiv:1805.10705v2, 2018.

The Matlab solver is called by

[Rsol,tsol,fsol,ksol] = solver_planar_p4pfr_fast(X,U)

Input:

X: 2x4 3D-coordinates (the 3D plane is assumed to be z=0) 

U: 2x4 image coordinates

Output:

returns n solutions, 

Rsol: contains cell 1xn of rotations

tsol: 3xn translations

fsol: 1xn focal lengths

ksol: 1xn radial distortion 
