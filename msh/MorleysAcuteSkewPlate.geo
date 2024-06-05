
a = 100.0;
n = 64;
b = 3^0.5*a/2;
Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {b+a, a/2 , 0.0};
Point(4) = {b, a/2 , 0.0};
Point(5) = {(b+a)/2, a/4 , 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {4,5};
Line(6) = {5,2};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,2,3,4,5,6} = n;
Transfinite Curve{5,6} = n/4+1;
//Transfinite Surface{1};

Point{5} In Surface{1};
Curve{5} In Surface{1};
Curve{6} In Surface{1};

Physical Curve("Γᵇ") = {1};
Physical Curve("Γʳ") = {2};
Physical Curve("Γᵗ") = {3};
Physical Curve("Γˡ") = {4};
Physical Surface("Ω") = {1};
Physical Point("𝐴") = {5};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;
