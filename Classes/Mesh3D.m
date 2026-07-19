classdef Mesh3D < MeshStructure
    %Mesh3D 3D Cartesian mesh. SEE ALSO: createMesh3D, MeshStructure
    methods
        function m = Mesh3D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {3, dims, cellsize, cellcenters, facecenters, corners, edges, 'cartesian'};
            end
            m = m@MeshStructure(args{:});
        end
    end
end
