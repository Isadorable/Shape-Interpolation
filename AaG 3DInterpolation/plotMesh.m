%Display the mesh
function plot = plotMesh(V, F)  
    axis off
    hold on
    plot = trisurf(F, V(:,1), V(:,2), V(:,3),'FaceColor','interp');
end