classdef Mesh2D < MeshStructure
    %Mesh2D 2D Cartesian mesh. SEE ALSO: createMesh2D, MeshStructure
    methods
        function m = Mesh2D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {2, dims, cellsize, cellcenters, facecenters, corners, edges, 'cartesian'};
            end
            m = m@MeshStructure(args{:});
        end
    end
end
