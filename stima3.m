function V=stima3(vertices,qin,nodes)
mk=1/2*det([ones(1,3);vertices']);
qin=qin(nodes);
S=zeros(6,1);
S(1,1)=qin(1,1)^2; S(4,1)=qin(1,1)*qin(1,2)*2;
S(2,1)=qin(1,2)^2; S(5,1)=qin(1,2)*qin(1,3)*2;
S(3,1)=qin(1,3)^2; S(6,1)=qin(1,1)*qin(1,3)*2;
T11=[24 4 4 2*6 2*2 2*6]*2*mk/720;
T12=[6 6 2 2*2*2 2*2 2*2]*2*mk/720;
T13=[6 2 6 2*2 2*2 2*4]*2*mk/720;
T22=[4 24 4 2*6 2*6 2*2]*2*mk/720;
T23=[2 6 6 2*2 2*4 2*2]*2*mk/720;
T33=[4 4 24 2*2 2*6 2*6]*2*mk/720;
T21=T12;T32=T23;T31=T13;
V=zeros(3,3);
V(1,1)=T11*S;V(1,2)=T12*S;
V(1,3)=T13*S;V(2,2)=T22*S;
V(2,3)=T23*S;V(3,3)=T33*S;
V(2,1)=V(1,2);V(3,2)=V(2,3);V(3,1)=V(1,3);

%%%%%%%% assembled matrix3 %%%%%%%%%%%%%
%%%  for j=1:size(Elem,1)
%       V(Elem(j,:),Elem(j,:))=V(Elem(j,:),Elem(j,:))+stima3(Coord(Elem(j,:),:),q2in,Elem(j,:));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5          
% q2in is the second component of the input vector
