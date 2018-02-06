%%%%%%%%%%%%%%%3D as-rigid-as-possible interpolation%%%%%%%%%%%%%%%%%%%%%%%
%Author: Isadora Sanna is430@bath.ac.uk
%The following program in an implementation of two famous methods used
%for shape interpolation and image morphing: Linear interpolation and
%As-rigid-as-possible interpolation from the paper "As-rigid-as-possible 
%shape interpolation by Alexa et al.
%Formulas listed in the comments are related to the report.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = ui3d(varargin)
%UI initialization
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ui3d_OpeningFcn, ...
                   'gui_OutputFcn',  @ui3d_OutputFcn, ...
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


% --- Executes just before ui3d is made visible.
function ui3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ui3d (see VARARGIN)

% Choose default command line output for ui3d
handles.output = hObject;

handles.poses = 2;

%%Load and display the SOURCE mesh
fileName = 'cube_source.obj';
[handles.F1, handles.V1] = loadFile(fileName);
plot = plotMesh(handles.V1,handles.F1);
set(handles.sourceAxes,'Visible','off');
set(plot,'Parent', handles.sourceAxes);

%%Load and display the TARGET mesh
fileName = 'cube_target.obj';
[handles.F2, handles.V2] = loadFile(fileName);
plot = plotMesh(handles.V2,handles.F2);
set(handles.targetAxes,'Visible','off');
set(plot,'Parent', handles.targetAxes);

%we use a cell to save the vertices set retrived from "loadFile()" in
%order to access their values more easily
handles.V{1} = handles.V1;
handles.V{2} = handles.V2;

%axes limits
handles.limXtarget = get(handles.targetAxes,'xlim');
handles.limYtarget = get(handles.targetAxes,'ylim');

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ui3d_OutputFcn(hObject, eventdata, handles) 
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
        plotInBetween(handles.F1,newV,handles.limXtarget,handles.limYtarget)
    end
end

%%%%%%%%%%%%%%%%%% As-Rigid-As-Possible INTERPOLATION %%%%%%%%%%%%%%%%%%%%%
function startARAPInterpolation(hObject)
handles = guidata(hObject);
figure('Name','As-Rigid-As-Possible interpolation','NumberTitle','off');
%For each pair of poses, we set the first one as the source pose and the
%second one as the target pose. We use the values of their respetives vertices
%to compute the position of the vertices v' in between. t is set to 0.01.
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
    u = zeros(size(P,1)*9,size(P,1)*3);
    pLeq = zeros(9);
    
    nTriangles = size(Faces,1);
    
    tempSolution = zeros(3);
    for t = 0:0.01:1
        %known vertex with interpolation
        vFirst = (1-t)*P(1,1:3)'+ t*Q(1,1:3)';
        
        for i = 1:nTriangles
            pLeq(1,1:3) = P(Faces(i,1),1:3);
            pLeq(2,4:6) = P(Faces(i,1),1:3);
            pLeq(3,7:9) = P(Faces(i,1),1:3);
            pLeq(4,1:3) = P(Faces(i,2),1:3);
            pLeq(5,4:6) = P(Faces(i,2),1:3);
            pLeq(6,7:9) = P(Faces(i,2),1:3);
            pLeq(7,1:3) = P(Faces(i,3),1:3);
            pLeq(8,4:6) = P(Faces(i,3),1:3);
            pLeq(9,7:9) = P(Faces(i,3),1:3);
            pLeq([1,4,7],10) = 1;
            pLeq([2,5,8],11) = 1;
            pLeq([3,6,9],12) = 1;
            qLeq = [Q(Faces(i,1),1);Q(Faces(i,1),2);Q(Faces(i,1),3);...
                Q(Faces(i,2),1);Q(Faces(i,2),2);Q(Faces(i,2),3);...
                Q(Faces(i,3),1);Q(Faces(i,3),2);Q(Faces(i,3),3)];
            RS = pLeq\qLeq;
            
            %Formula 2.8 and 2.9
            tempSolution(1:3,1:3) = [RS(1:3)';RS(4:6)';RS(7:9)'];
            [U,D,V] = svd(tempSolution);
            Rgamma = U*V';
            S = V*D*V';
            
            %Formula 2.10
            A = Rgamma*((1-t)*eye(size(Rgamma,1)) + t*S);
            A = A(1:3,1:3);
            
            %Creation of 2.18 Chapter 2.3
            %We use the values of inv(pLeq) as entries for the matrix u
            pInv = pLeq\eye(size(pLeq));
            
            u(i*9-8:i*9-6,Faces(i,1)*3-2) = pInv(1:3,1);
            u(i*9-8:i*9-6,Faces(i,2)*3-2) = pInv(1:3,4);
            u(i*9-8:i*9-6,Faces(i,3)*3-2) = pInv(1:3,7);
            
            u(i*9-5:i*9-3,Faces(i,1)*3-1) = pInv(4:6,2);
            u(i*9-5:i*9-3,Faces(i,2)*3-1) = pInv(4:6,5);
            u(i*9-5:i*9-3,Faces(i,3)*3-1) = pInv(4:6,8);
            
            u(i*9-2:i*9,Faces(i,1)*3) = pInv(7:9,3);
            u(i*9-2:i*9,Faces(i,2)*3) = pInv(7:9,6);
            u(i*9-2:i*9,Faces(i,3)*3) = pInv(7:9,9);
            
            b(i*9-8:i*9) = [A(1,1)-dot(vFirst,u(i*9-8,1:3));A(1,2)-dot(vFirst,u(i*9-7,1:3));A(1,3)-dot(vFirst,u(i*9-6,1:3));...
                A(2,1)-dot(vFirst,u(i*9-5,1:3));A(2,2)-dot(vFirst,u(i*9-4,1:3));A(2,3)-dot(vFirst,u(i*9-3,1:3));...
                A(3,1)-dot(vFirst,u(i*9-2,1:3));A(3,2)-dot(vFirst,u(i*9-1,1:3));A(3,3)-dot(vFirst,u(i*9,1:3))];
        end
        %We remove the first three columns of u and then we compute the
        %solution of the linear system ux = b. We set as first value of the
        %solution the previously fixed first vertex.
        ui = u(:,4:end);
        ui = sparse(ui);
        newV = ui\b;
        newV = [vFirst;newV];
        Vx(:,1) = newV(1:3:end,:);
        Vy(:,1) = newV(2:3:end,:);
        Vz(:,1) = newV(3:3:end,:);
        newV = [Vx,Vy,Vz];
        
        %Plot each frame in between
        plotInBetween(handles.F1,newV,handles.limXtarget,handles.limYtarget)
    end
end

function plotInBetween(F,V,limX,limY)
    clf;
    hold on;
    view(2);
    axis off
    axis([limX(1) limX(2) limY(1) limY(2)]);
    trisurf(F, V(:,1), V(:,2), V(:,3),'FaceColor','interp');
    pause(0.05);
    drawnow nocallbacks
