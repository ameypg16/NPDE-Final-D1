
%%%function for calculating the various delta q errors%%%%%%

function ufinal=fem3(q1in,q2in)
[Coord,Elem,Nb,Db]=InitialMesh2(1);
epsilon=0.02;


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
% U1=sparse(size(Coord,1),size(Coord,1));
% qin=fem2();
% q1in=zeros(1,13);
% q2in=zeros(1,13);
% for i=1:13
%     q1in(i)=qin(i);
%     q2in(i)=qin(13+i);
% end
    

%%%%%%%% assembled matrix2 %%%%%%%%%%%%%
% for j=1:size(Elem,1)
%        U(Elem(j,:),Elem(j,:))=U(Elem(j,:),Elem(j,:))+3*stima2(Coord(Elem(j,:),:),q1in,Elem(j,:));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% alpha        
% qin is the input vector which will be updated as the iteration goes on.
U1=sparse(size(Coord,1),size(Coord,1));
V1=sparse(size(Coord,1),size(Coord,1));
W1=sparse(size(Coord,1),size(Coord,1));
X1=sparse(size(Coord,1),size(Coord,1));
U2=sparse(size(Coord,1),size(Coord,1));
V2=sparse(size(Coord,1),size(Coord,1));
W2=sparse(size(Coord,1),size(Coord,1));
X2=sparse(size(Coord,1),size(Coord,1));
%%%%%%%% assembled matrix3 %%%%%%%%%%%%%
  for j=1:size(Elem,1)
       V1(Elem(j,:),Elem(j,:))=V1(Elem(j,:),Elem(j,:))+stima3(Coord(Elem(j,:),:),q2in,Elem(j,:));
       U1(Elem(j,:),Elem(j,:))=U1(Elem(j,:),Elem(j,:))+3*stima2(Coord(Elem(j,:),:),q1in,Elem(j,:));
       W1(Elem(j,:),Elem(j,:))=W1(Elem(j,:),Elem(j,:))+stima4(Coord(Elem(j,:),:));
       X1(Elem(j,:),Elem(j,:))=X1(Elem(j,:),Elem(j,:))+2*stima5(Coord(Elem(j,:),:),q1in,q2in,Elem(j,:));
       V2(Elem(j,:),Elem(j,:))=V2(Elem(j,:),Elem(j,:))+stima3(Coord(Elem(j,:),:),q1in,Elem(j,:));
       U2(Elem(j,:),Elem(j,:))=U2(Elem(j,:),Elem(j,:))+3*stima2(Coord(Elem(j,:),:),q2in,Elem(j,:));
       W2(Elem(j,:),Elem(j,:))=W2(Elem(j,:),Elem(j,:))+stima4(Coord(Elem(j,:),:));
       X2(Elem(j,:),Elem(j,:))=X2(Elem(j,:),Elem(j,:))+2*stima5(Coord(Elem(j,:),:),q2in,q1in,Elem(j,:));       
  end
  
  B1=(2/(epsilon)^2)*(U1+V1-W1);
  B2=(2/(epsilon)^2)*(U2+V2-W2);
  C=(2/(epsilon)^2)*2*X1;
%   C2=(2/(epsilon)^2)*2*X2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5     beta     
% qin is the input vector
%%%%%%%% assembled matrix4 %%%%%%%%%%%%%
% for j=1:size(Elem,1)
%      W(Elem(j,:),Elem(j,:))=W(Elem(j,:),Elem(j,:))+stima4(Coord(Elem(j,:),:));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5          
% qin is the input vector
%%%%%%%% assembled matrix5 %%%%%%%%%%%%%
% for j=1:size(Elem,1)
%      X(Elem(j,:),Elem(j,:))=X(Elem(j,:),Elem(j,:))+2*stima5(Coord(Elem(j,:),:),q1in,q2in,Elem(j,:));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5          
% qin=[q1in,q2in] is the input vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
% % % % % % % % % % final assembly
% stima1+(2/e^2)*(U+V-W+X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555


Y1=sparse(size(Coord,1),1);
Y2=sparse(size(Coord,1),1);
Zeq1=sparse(size(Coord,1),1);
Zeq2=sparse(size(Coord,1),1);
Z1eq1=sparse(size(Coord,1),1);
Z2eq1=sparse(size(Coord,1),1);
Z1eq2=sparse(size(Coord,1),1);
Z2eq2=sparse(size(Coord,1),1);
%%%%%assembly of rhs
for j=1:size(Elem,1),
    Y1(Elem(j,:))=Y1(Elem(j,:))+rhs1(Coord(Elem(j,:),:),q1in,Elem(j,:));
    Y2(Elem(j,:))=Y2(Elem(j,:))+rhs1(Coord(Elem(j,:),:),q2in,Elem(j,:));
    Zeq1(Elem(j,:))=Zeq1(Elem(j,:))+rhs2(Coord(Elem(j,:),:),q2in,q1in,Elem(j,:));
    Zeq2(Elem(j,:))=Zeq2(Elem(j,:))+rhs2(Coord(Elem(j,:),:),q1in,q2in,Elem(j,:));
    Z1eq1(Elem(j,:))=Z1eq1(Elem(j,:))+rhs3(Coord(Elem(j,:),:),q1in,Elem(j,:));
    Z2eq1(Elem(j,:))=Z2eq1(Elem(j,:))+rhs4(Coord(Elem(j,:),:),q1in,Elem(j,:));
    Z1eq2(Elem(j,:))=Z1eq2(Elem(j,:))+rhs3(Coord(Elem(j,:),:),q2in,Elem(j,:));
    Z2eq2(Elem(j,:))=Z2eq2(Elem(j,:))+rhs4(Coord(Elem(j,:),:),q2in,Elem(j,:));
end

D1=-Z2eq1-(2/epsilon^2)*(Y1+Zeq1-Z1eq1);
D2=-Z2eq2-(2/epsilon^2)*(Y2+Zeq2-Z1eq2);


%%%%%%%%%%%%%%
% final rhs=-Z2-(2/e^2)*(Y+Z-Z1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembly of load vector b
% % % % % for j=1:size(Elem,1),
% % %     b(Elem(j,:))=b(Elem(j,:))+det([1 1 1; Coord(Elem(j,:),:)'])*f(sum(Coord(Elem(j,:),:))/3)/6;
% % % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann Conditions
if (~isempty(Nb))
   for j=1:size(Nb,1)
       b(Nb(j,:))= b(Nb(j,:))+norm(Coord(Nb(j,1),:)-Coord(Nb(j,2),:))*...
                                     u_N(Coord(Nb(j,1),:),Coord(Nb(j,2),:))/2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uh=zeros(length(FullNodes),1);
uh1=zeros(size(Coord,1),1);
uh2=zeros(size(Coord,1),1);
%uh
% Dirichlet Conditions
% % % if (~isempty(Db))
% % %     Dbnodes=unique(Db);                                 %Db contains the information of g
% % %     for j=1:size(Dbnodes,1)
% % %        uh(Dbnodes(j),1)=ue(Coord(Dbnodes(j),:));
% % %     end
% % % end
% % % b=b-A*uh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Solving the linear system
% uh(FreeNodes)=A(FreeNodes,FreeNodes)\b(FreeNodes);
% uh

uh2=(C*D1-(A+B1)*D2)\(C*C-(A+B2)*(A+B1));
uh1=(C*D2-(A+B2)*D1)\(C*C-(A+B1)*(A+B2));
uh1
uh2
% Exact solution at the nodes
% u=u_nodes(Coord);

% Display the computed solution
% figure(1)
% show(Coord,Elem,uh,u)
ufinal=[uh1 uh2];
% uhh1=zeros(1,size(Coord,1));
% uhh2=zeros(1,size(Coord,1));
% for i=1:size(Coord,1)
%     uhh1(1,i)=uh1(i,1);
%     uhh2(1,i)=uh2(i,1);
% end

ufinal
size(ufinal)

