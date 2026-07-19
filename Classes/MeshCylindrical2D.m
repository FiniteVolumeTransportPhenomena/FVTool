classdef MeshCylindrical2D < Mesh2D
    %MeshCylindrical2D 2D axisymmetric cylindrical mesh (r, z).
    % SEE ALSO: createMeshCylindrical2D, MeshStructure
    methods
        function m = MeshCylindrical2D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {dims, cellsize, cellcenters, facecenters, corners, edges};
            end
            m = m@Mesh2D(args{:});
            if nargin>0
                m.coordsystem = 'cylindrical';
            end
        end
    end
end
