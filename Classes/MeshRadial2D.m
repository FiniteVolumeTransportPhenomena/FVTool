classdef MeshRadial2D < Mesh2D
    %MeshRadial2D 2D polar/radial mesh (r, theta).
    % SEE ALSO: createMeshRadial2D, MeshStructure
    methods
        function m = MeshRadial2D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {dims, cellsize, cellcenters, facecenters, corners, edges};
            end
            m = m@Mesh2D(args{:});
            if nargin>0
                m.coordsystem = 'radial';
            end
        end
    end
end
