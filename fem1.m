%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM for  - div . (u_x,u_y) + u  = f    in  \Omega             %
%                   (u_x, u_y). n = u_N  on  \partial\Omega_N   %
%                               u = g    on  \partial\Omega_D   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
format long;

% Geometry of finite elements - triangulation, local global node numbering,
% coordinates, boundary nodes
%%%%%

[Coord,Elem,Nb,Db]=InitialMesh(1);
for j=1:1
%% Edge-Node-Element Connections
[n2ed,ed2el]=edge(Elem,Coord);
%% Element Redrefine
[Coord,Elem,Db,Nb]=redrefine(Coord,Elem,n2ed,ed2el,Db,Nb);
end
% No of degrees of freedom (initially solution at all the nodes
%are assumed as unknowns, dirichlet boundary conditions 
% if any are to be incorporated later on )
FullNodes=[1:size(Coord,1)];                        % define the components of uh to be present or not
FreeNodes=setdiff(FullNodes, unique(Db));
%g=0
FullNodes
FreeNodes
% Intializing the matrices
A=sparse(size(Coord,1),size(Coord,1)); % A is global stiffness matrix
b=sparse(size(Coord,1),1); % global load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly of A 
% stima is element stiffness matrices
for j=1:size(Elem,1),
    A(Elem(j,:),Elem(j,:))=A(Elem(j,:),Elem(j,:))+...
                                    stima(Coord(Elem(j,:),:));                                                                
end
% Assembly of load vector b
for j=1:size(Elem,1),
    b(Elem(j,:))=b(Elem(j,:))+det([1 1 1; Coord(Elem(j,:),:)'])*f(sum(Coord(Elem(j,:),:))/3)/6;
end
% Neumann Conditions
if (~isempty(Nb))
   for j=1:size(Nb,1)
       b(Nb(j,:))= b(Nb(j,:))+norm(Coord(Nb(j,1),:)-Coord(Nb(j,2),:))*...
                                     u_N(Coord(Nb(j,1),:),Coord(Nb(j,2),:))/2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uh=zeros(length(FullNodes),1);
uh
% Dirichlet Conditions
 if (~isempty(Db))
     Dbnodes=unique(Db);                                 %Db contains the information of g
     for j=1:size(Dbnodes,1)
        uh(Dbnodes(j),1)=ue(Coord(Dbnodes(j),:));
     end
 end
b=b-A*uh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the linear system
uh(FreeNodes)=A(FreeNodes,FreeNodes)\b(FreeNodes);
uh

% Exact solution at the nodes
u=u_nodes(Coord);
% Display the computed solution
 figure(1)
 show(Coord,Elem,uh,u)

