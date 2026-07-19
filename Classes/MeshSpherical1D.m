classdef MeshSpherical1D < Mesh1D
    %MeshSpherical1D 1D spherical (radial) mesh.
    % SEE ALSO: createMeshSpherical1D, MeshStructure
    methods
        function m = MeshSpherical1D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {dims, cellsize, cellcenters, facecenters, corners, edges};
            end
            m = m@Mesh1D(args{:});
            if nargin>0
                m.coordsystem = 'spherical';
            end
        end
    end
end
