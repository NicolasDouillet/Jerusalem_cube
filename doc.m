%% Jerusalem_cube
%
% Function to compute, display, and save the
% Jerusalem cube at any iteration / depth level.
%
% Author & support : nicolas (dot) douillet (at) free.fr, 2021.
%
%% Syntax
%
% Jerusalem_cube;
% 
% Jerusalem_cube(nb_iterations);
%
% Jerusalem_cube(nb_iterations, printable_ready);
%
% Jerusalem_cube(nb_iterations, printable_ready, option_display);
%
% [V, T] = Jerusalem_cube(nb_iterations, printable_ready, option_display);
%
%% Description
%
% Jerusalem_cube computes and displays the 2-
% Jerusalem cube included in the unit sphere.
% 
% Jerusalem_cube(nb_it) computes and displays
% the nb_it Jerusalem cube.
%
% Jerusalem_cube(nb_iterations, printable_ready) prevents from
% creating non manifold edges when printable_ready is set to *true /
% logical 1, and remove duplicated vertices and faces when it is set to
% false / logical 0. In this latter case, the model is lighter (less
% vertices, less faces), but at the cost of non manifoldness.
%
% Jerusalem_cube(nb_iterations, printable_ready, option_display)
% displays it when option_display is set to logical *true/1 (default),
% and doesn't when it is set to  logical false/0.
%
% [V, T] = Jerusalem_cube(nb_iterations, printable_ready, option_display) saves
% the resulting vertex coordinates in the array V, and the triangulation in the array T.
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/73447-sierpinski-menger-sponge?s_tid=prof_contriblnk
% Sierpinski-Menger sponge> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73216-cubic-based-3d-koch-snowflake?s_tid=prof_contriblnk cubic_based_3D_koch_snowflake>
%
%% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - printable_ready : either logical, true/*false or numeric 1/*0.
%
% - option_display : either logical, *true/false or numeric *1/0.
%
%% Output arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - T = [i1 i2 i3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Computes and displays the Jerusalem cube at iteration 2, with minimum vertex and face numbers

Jerusalem_cube;

%% Example #2
% Computes and saves the 3D printable ready Jerusalem cube at iteration 3

[V,T] = Jerusalem_cube(3,true,false);