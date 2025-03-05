n = 64;
//a = 1/n;
//b = 1.0-1/n;


//Point(1) = {a, a, 0.0};
//Point(2) = {  b, a, 0.0};
//Point(3) = {  b,   b, 0.0};
//Point(4) = {a,   b, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,2,3,4} = n-1;
Transfinite Surface{1};
Physical Curve("Γᵇ") = {1};
Physical Curve("Γᵗ") = {2};
Physical Curve("Γˡ") = {3};
Physical Curve("Γʳ") = {4};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
//Mesh.SecondOrderIncomplete = 1;
Mesh 2;
//RecombineMesh;
//SetOrder 2;
