
a = 5.0;
n = 27;
b = a/(2*32);

Point(1) = {b, b, 0.0};
Point(2) = {(a-b)/2, b, 0.0};
Point(3) = {  a-b, b, 0.0};
Point(4) = {b+(a-2*b)/2^0.5, b+(a-2*b)/2^0.5, 0.0};
Point(5) = {b,   a-b, 0.0};
Point(6) = {b, (a-b)/2, 0.0};
Point(7) = {2*b+a/2^0.5/2, 2*b+a/2^0.5/2, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Circle(3) = {3,1,4};
Circle(4) = {4,1,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {7,2};
Line(8) = {7,6};
Line(9) = {7,4};

Curve Loop(1) = {1,-7,8,6};
Curve Loop(2) = {2,3,-9,7};
Curve Loop(3) = {9,4,5,-8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Transfinite Curve{1,2,3,4,5,6,7,8,9} = n;
Transfinite Surface{2} ;
Transfinite Surface{1,3} Right;

Physical Curve("Γᵇ") = {1,2};
Physical Curve("Γᵉ") = {3,4};
Physical Curve("Γˡ") = {5,6};
Physical Surface("Ω") = {1,2,3};
Physical Point("𝐴") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
//Mesh.SecondOrderIncomplete = 1;
Mesh 2;
RecombineMesh;
//SetOrder 2;
