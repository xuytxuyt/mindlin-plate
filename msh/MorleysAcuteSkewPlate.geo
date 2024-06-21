
a = 100.0;
n = 8;
b = 3^0.5*a/2;
Point(1) = {0.0, 0.0, 0.0};
Point(2) = {a/2, 0.0, 0.0};
Point(3) = {  a, 0.0, 0.0};
Point(4) = {b+a, a/2 , 0.0};
Point(5) = {b+a/2, a/2 , 0.0};
Point(6) = {b, a/2 , 0.0};
Point(7) = {(b+a)/2, a/4 , 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,7};
Line(8) = {7,5};


Curve Loop(1) = {1,2,3,4,5,6};

Plane Surface(1) = {1};

Transfinite Curve{3,6} = n;
Transfinite Curve{1,2,4,5,7,8} = n/2+1;
Transfinite Surface{1};

Curve{7} In Surface{1};
Curve{8} In Surface{1};

Physical Curve("Γᵇ") = {1,2};
Physical Curve("Γʳ") = {3};
Physical Curve("Γᵗ") = {4,5};
Physical Curve("Γˡ") = {6};
Physical Surface("Ω") = {1};
Physical Point("𝐴") = {7};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;
