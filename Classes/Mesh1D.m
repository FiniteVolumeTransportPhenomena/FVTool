classdef Mesh1D < MeshStructure
    %Mesh1D 1D Cartesian mesh. SEE ALSO: createMesh1D, MeshStructure
    methods
        function m = Mesh1D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {1, dims, cellsize, cellcenters, facecenters, corners, edges, 'cartesian'};
            end
            m = m@MeshStructure(args{:});
        end
    end
end
