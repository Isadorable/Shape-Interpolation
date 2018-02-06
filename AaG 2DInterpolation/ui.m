%%%%%%%%%%%%%%%2D as-rigid-as-possible interpolation%%%%%%%%%%%%%%%%%%%%%%%
%Author: Isadora Sanna is430@bath.ac.uk
%The following program in an implementation of two famous methods used
%for shape interpolation and image morphing: Linear interpolation and
%As-rigid-as-possible interpolation from the paper "As-rigid-as-possible 
%shape interpolation by Alexa et al.
%Formulas listed in the comments are related to the report.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = ui(varargin)
%UI initialization
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ui_OpeningFcn, ...
                   'gui_OutputFcn',  @ui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code 


% --- Executes just before ui is made visible.
function ui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ui (see VARARGIN)

%update the "global" variables
handles.output = hObject;

%Load source, target and intermediate meshes
%the number of meshes is set in "handles.poses".
%In case we want to use more or less than 3 meshes we need to modify the
%value of "handles.poses" and add another sequence of the following 5 lines
handles.poses = 3;

%%Load and display the SOURCE mesh
fileName = 'man.obj';
[handles.F1, handles.V1] = loadFile(fileName);
plot = plotMesh(handles.V1,handles.F1);
set(handles.sourceAxes,'Visible','off');
set(plot,'Parent', handles.sourceAxes);

%%Load and display the INTERMEDIATE mesh
fileName = 'man2.obj';
[handles.F2, handles.V2] = loadFile(fileName);
plot = plotMesh(handles.V2,handles.F2);
set(handles.second,'Visible','off');
set(plot,'Parent', handles.second);

%%Load and display the TARGET mesh
fileName = 'man3.obj';
[handles.F3, handles.V3] = loadFile(fileName);
plot = plotMesh(handles.V3,handles.F3);
set(handles.targetAxes,'Visible','off');
set(plot,'Parent', handles.targetAxes);

%we use a cell to save the vertices set retrived from "loadFile()" in
%order to access their values more easily
handles.V{1} = handles.V1;
handles.V{2} = handles.V2;
handles.V{3} = handles.V3;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startNormalInterpolation(hObject);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startARAPInterpolation(hObject);


%%%%%%%%%%%%%%%%%%%% LINEAR INTERPOLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startNormalInterpolation(hObject)
handles = guidata(hObject);
%set the title of the window
figure('Name','Linear interpolation','NumberTitle','off');

%For each pair of poses, we set the first one as the source pose and the
%second one as the target pose. We use the values of their respetives vertices
%to compute the position of the vertices v' in between. t is set to 0.01.
for p = 1:handles.poses-1
newV = zeros(size(handles.V{p},1),3);
    for t = 0:0.01:1
        for i = 1 : size(handles.V{p},1)
            %apply the standard formula (1-t)*v0 + t*v1
            v =  (1-t)*[handles.V{p}(i,1);handles.V{p}(i,2);handles.V{p}(i,3)]+...
                t*[handles.V{p+1}(i,1);handles.V{p+1}(i,2);handles.V{p+1}(i,3)];
            newV(i,:) = v';
        end

        %Plot each frame in between
        plotInBetween(handles.F1,newV)
    end
end

%%%%%%%%%%%%%%%%%% As-Rigid-As-Possible INTERPOLATION %%%%%%%%%%%%%%%%%%%%%
function startARAPInterpolation(hObject)
handles = guidata(hObject);
%set the title of the window
figure('Name','As-Rigid-As-Possible interpolation','NumberTitle','off');

%For each pair of poses, we set the first one as the source pose and the
%second one as the target pose. We use the values of their respetives vertices
%to compute the position of the vertices v' in between. t is set to 0.05.
for p = 1:handles.poses-1
    
    %Taking in account the report:
    %-Faces = triangles in the mesh. Info about the links between the vertices
    %-P = set of source vertices
    %-Q = set of target vertices
    %-pLeq and qLeq = respectively the matrix of known coordinates from the
    %source mesh and the column vector made with the known coordinates of
    %the target mesh. Description : Chapter 2.3 Formula 2.7
    %-b and u are respectively the matrix and column vector in chapter 2.3
    %Formula 2.18
    Faces = handles.F1;
    P = handles.V{p};
    Q = handles.V{p+1};
    b = zeros(size(P,1)*4,1);
    u = zeros(size(P,1)*4,size(P,1));
    pLeq = zeros(6);
    
    nTriangles = size(handles.F1,1);
    
    tempSolution = zeros(3);
    for t = 0:0.05:1
        %known vertex with interpolation
        vFirst = (1-t)*[P(1,1);P(1,2)]+ t*[Q(1,1);Q(1,2)];
                
        for i = 1:nTriangles
            %Formula 2.7
            pLeq(1,1:2) = P(Faces(i,1),1:2);
            pLeq(2,3:4) = P(Faces(i,1),1:2);
            pLeq(3,1:2) = P(Faces(i,2),1:2);
            pLeq(4,3:4) = P(Faces(i,2),1:2);
            pLeq(5,1:2) = P(Faces(i,3),1:2);
            pLeq(6,3:4) = P(Faces(i,3),1:2);
            pLeq([1,3,5],5) = 1;
            pLeq([2,4,6],6) = 1;
            qLeq = [Q(Faces(i,1),1);Q(Faces(i,1),2);Q(Faces(i,2),1);Q(Faces(i,2),2);Q(Faces(i,3),1);Q(Faces(i,3),2)];
            RS = pLeq\qLeq;
            
            %Formula 2.8 and 2.9
            tempSolution(1:2,1:2) = [RS(1:2,1)';RS(3:4,1)'];
            [U,D,V] = svd(tempSolution);
            Rgamma = U*V';
            S = V*D*V';
            
            %Formula 2.11
            %A = Rgamma*((1-t)*eye(size(Rgamma,1)) + t*S);
            A = ((1-t)*eye(3) + t*Rgamma) * ((1-t)*eye(3) + t*S);
            A = A(1:2,1:2);
            
            %Creation of 2.18 Chapter 2.3
            %We use the values of inv(pLeq) as entries for the matrix u
            pInv = inv(pLeq);
            
            u(i*4-3,Faces(i,1)*2-1) = pInv(1,1);
            u(i*4-3,Faces(i,2)*2-1) = pInv(1,3);
            u(i*4-3,Faces(i,3)*2-1) = pInv(1,5);
            
            u(i*4-2,Faces(i,1)*2-1) = pInv(2,1);
            u(i*4-2,Faces(i,2)*2-1) = pInv(2,3);
            u(i*4-2,Faces(i,3)*2-1) = pInv(2,5);
            
            u(i*4-1,Faces(i,1)*2) = pInv(3,2);
            u(i*4-1,Faces(i,2)*2) = pInv(3,4);
            u(i*4-1,Faces(i,3)*2) = pInv(3,6);
            
            u(i*4,Faces(i,1)*2) = pInv(4,2);
            u(i*4,Faces(i,2)*2) = pInv(4,4);
            u(i*4,Faces(i,3)*2) = pInv(4,6);
            %We use the values of A as entries for the column vector b
            %considering that the first vertex is fixed (Formula 2.19)
            b(i*4-3:i*4) = [A(1,1)-dot(vFirst,u(i*4-3,1:2));A(1,2)-dot(vFirst,u(i*4-2,1:2));...
                A(2,1)-dot(vFirst,u(i*4-1,1:2));A(2,2)-dot(vFirst,u(i*4,1:2))];
            
        end
        %We remove the first two columns of u and then we compute the
        %solution of the linear system ux = b. We set as first value of the
        %solution the previously fixed first vertex.
        ui = u(:,3:end);
        newV = ui\b;
        newV = [vFirst;newV];
        Vx(:,1) = newV(1:2:end,:);
        Vy(:,1) = newV(2:2:end,:);
        newV = [Vx,Vy];
        
        %Plot each frame in between
        plotInBetween(handles.F1,newV)
 
    end
end

function plotInBetween(F,V)
    clf;
    hold on;
    %2D view
    view(2);
    axis off
    %[xmin xmax ymin ymax]
    axis([-2 2 -1.5 2]);
    trisurf(F, V(:,1), V(:,2), zeros(size(V,1),1),'FaceColor','interp');
    pause(0.01);
    drawnow nocallbacks

