%Display the mesh
function plot = plotMesh(V, F) 
    %2D view
    view(2);
    axis off
    hold on
    plot = trisurf(F, V(:,1), V(:,2), zeros(size(V,1),1),'FaceColor','interp');
end