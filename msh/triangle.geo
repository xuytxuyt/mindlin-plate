
a = 10;
lc = a/3^0.5/1.5;

Point(1) = {2*a/3,           0, 0, lc};
Point(2) = {- a/3,  a/3^0.5, 0, lc};
Point(3) = {- a/3, -a/3^0.5, 0, lc};

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,1} ;

Line Loop(1) = {1,2,3};

Plane Surface(1) = {1};

Physical Line("Γ₁") = {1};
Physical Line("Γ₂") = {2};
Physical Line("Γ₃") = {3};

Physical Surface("Ω") = {1};


Mesh.Algorithm = 6;
// Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RefineMesh;
// RecombineMesh;
