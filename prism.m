%%%%prism form-finding
clc;clear all;close all
%material properties
E_bar=1e3*eye(12);
A_bar=diag([1000*ones(3,1);ones(9,1)]);

%geometry
R=5;
r=5;
h=5;
theta=pi/6;
N1=[R*cos(pi/6),R*sin(pi/6) 0
   -R*cos(pi/6),R*sin(pi/6) 0
   0 -R 0]';
T=[cos(theta) -sin(theta) 0
    sin(theta) cos(theta) 0
    0 0 1];
N2=T*[r*cos(pi/6),r*sin(pi/6) h
   -r*cos(pi/6),r*sin(pi/6) h
   0 -R h]';
N=[N1,N2];

%connectivity matrix
C_b_in = [1 5;2 6;3 4];  % Bar in circle
C_s_in=[1 4;2 5;3 6;1 2;2 3;3 1;4 5;5 6;6 4];

% Convert above index notation into actual connectivity matrices
C_b = tenseg_ind2C(C_b_in,N);
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
n=size(N,2); %number of node
ne=size(C,1);%number of element
%plot
tenseg_plot(N,C_b,C_s)
%% self-stress design
A=kron(C',eye(3))*diag(kron(C,eye(3))*N(:))*...
    kron(eye(ne),ones(3,1));  %equilibrium matrix
[U,S,V] = svd(A);
r=rank(A); %rank of (A_hat*Gp)

V1=V(:,1:r);V2=V(:,r+1:end);

t=1e1*V2; % force density
%%length and restlength design
L=sqrt(sum((N*C').^2))';    %length
l_bar=diag(L);
f=L.*t;           %force
f_bar=diag(f);

l0_bar=E_bar*A_bar*l_bar/(f_bar+E_bar*A_bar);% rest length
disp('length is')
disp(L')
disp('rest length is')
disp(diag(l0_bar)')
disp('l-l0 is')
disp(L'-diag(l0_bar)')

%% calculate the connecting point in plate
S=sqrt(2-sqrt(3))*R
alpha=acos((3*R^2+S^2-2*R^2)/(2*sqrt(3)*R*S))
P=S*sin(alpha)/sin(2*alpha+pi/3)
Q=P*(sqrt(2)*R-S)/(sqrt(2)*R)    %Q is the distance from pinned node to top or bottom node
