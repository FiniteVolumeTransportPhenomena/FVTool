classdef MeshSpherical3D < Mesh3D
    %MeshSpherical3D 3D spherical mesh (r, theta, phi).
    % SEE ALSO: createMeshSpherical3D, MeshStructure
    methods
        function m = MeshSpherical3D(dims, cellsize, cellcenters, facecenters, corners, edges)
            if nargin==0
                args = {};
            else
                args = {dims, cellsize, cellcenters, facecenters, corners, edges};
            end
            m = m@Mesh3D(args{:});
            if nargin>0
                m.coordsystem = 'spherical';
            end
        end
    end
end
