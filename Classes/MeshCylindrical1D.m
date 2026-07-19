classdef MeshCylindrical1D < Mesh1D
    %MeshCylindrical1D 1D cylindrical (radial) mesh.
    % SEE ALSO: createMeshCylindrical1D, MeshStructure
    methods
        function m = MeshCylindrical1D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {dims, cellsize, cellcenters, facecenters, corners, edges};
            end
            m = m@Mesh1D(args{:});
            if nargin>0
                m.coordsystem = 'cylindrical';
            end
        end
    end
end
