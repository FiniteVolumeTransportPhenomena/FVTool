classdef MeshCylindrical3D < Mesh3D
    %MeshCylindrical3D 3D cylindrical mesh (r, theta, z).
    % SEE ALSO: createMeshCylindrical3D, MeshStructure
    methods
        function m = MeshCylindrical3D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {dims, cellsize, cellcenters, facecenters, corners, edges};
            end
            m = m@Mesh3D(args{:});
            if nargin>0
                m.coordsystem = 'cylindrical';
            end
        end
    end
end
