classdef MeshStructure
    %MeshStructure base class for all FVTool meshes.
    % Holds the geometric description of the domain and the mesh. The
    % spatial dimension and the coordinate system are stored as two
    % explicit fields (an integer `dimension` and a string `coordsystem`)
    % rather than the old floating-point dimension codes. Concrete meshes
    % are the subclasses Mesh1D/Mesh2D/Mesh3D and their curvilinear
    % variants (MeshCylindrical*, MeshRadial2D, MeshSpherical*).
    %
    % SEE ALSO: createMesh1D, createMesh2D, createMesh3D,
    %   createMeshCylindrical1D, createMeshCylindrical2D,
    %   createMeshCylindrical3D, createMeshRadial2D,
    %   createMeshSpherical1D, createMeshSpherical3D, geometryTag

    properties
        dimension    % integer spatial dimension: 1, 2 or 3
        coordsystem  % 'cartesian' | 'cylindrical' | 'radial' | 'spherical'
        dims
        cellsize
        cellcenters
        facecenters
        corners
        edges
    end

    methods
        function meshVar = MeshStructure(dimension, dims, cellsize, ...
          cellcenters, facecenters, corners, edges, coordsystem)
            if nargin>0
                meshVar.dimension = dimension;
                meshVar.dims = dims;
                meshVar.cellsize = cellsize;
                meshVar.cellcenters = cellcenters;
                meshVar.facecenters = facecenters;
                meshVar.corners= corners;
                meshVar.edges= edges;
                if nargin>7
                    meshVar.coordsystem = coordsystem;
                else
                    meshVar.coordsystem = 'cartesian';
                end
            end
        end

        function tag = geometryTag(m)
            % GEOMETRYTAG canonical dispatch key for a mesh, e.g. '1D',
            % '2D', 'Cylindrical2D', 'Radial2D', 'Spherical3D'. Cartesian
            % meshes use the bare dimension suffix; curvilinear meshes are
            % prefixed with the capitalised coordinate system. The tag
            % matches the suffix of the geometry-specific implementation
            % files so dispatchers can switch on it directly.
            if strcmp(m.coordsystem, 'cartesian')
                tag = sprintf('%dD', m.dimension);
            else
                coord = [upper(m.coordsystem(1)) m.coordsystem(2:end)];
                tag = sprintf('%s%dD', coord, m.dimension);
            end
        end
    end
end
