%This function extracts information from the .obj file and organize them
%in 3 data structures:
%Vertices: is a nx2 matrix where n represents the number of vertices composing the matrix
%so that each row is a different vertex and the two columns contain the x and y coordinates
%for the corresponding vertex.
%Edges: is a nx2 matrix where n represents the number of edges in the mesh and the
%columns contain the index of the two vertices in the Vertices matrix that compose the
%edge.
%Faces: is a nx3 matrix where n represents the number of triangles in the mesh. Each
%column contains the index of the 3 vertices at the 3 corners of the triangle.
function [triangles, vertices] = loadFile(file)
mesh = fopen(file,'r');
verticesI = 0;
trianglesI = 0;

while ~feof(mesh)
    %Read new line from file
    line = fgetl(mesh);
    switch line(1)
        case '#' % number of vertices
            %memory preallocation for the structures
            if (strcmp(line(2:9),'vertices'))
                numVerteces = line(11:end);
                numVerteces = str2num(numVerteces);
                %creating a matrix with x,y,z column for each vertex
                vertices = zeros(numVerteces, 3);
            elseif(strcmp(line(2:6),'faces'))
                numtriangles = line(9:end);
                numtriangles = str2num(numtriangles);
                %creating a matrix with references to the three points part
                %of each triengle
                triangles = zeros(numtriangles, 3);
            end
        case 'v' % vertices
            verticesI = verticesI + 1;
            vertices(verticesI,:) = sscanf(line, 'v %f %f %f');
        case 'f' % triangles
            trianglesI = trianglesI + 1;
            triangles(trianglesI,:) = sscanf(line, 'f %d %d %d');
        otherwise
            fprintf('comment');
    end
end
fclose(mesh);

end