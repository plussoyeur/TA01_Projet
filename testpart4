h1 = 0.2;
h2 = 0.1;

Point(1) = {0,0,0,h1};
Point(2) = {12,0,0,h1};
Point(3) = {12,4,0,h1};
Point(4) = {0,4,0,h1};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};




// Generation d'un cercle
Point(5) = {1,2,0,h1};
Point(6) = {3,2,0,h1};
Point(7) = {2,2,0,h1};

Circle(5) = {5,7,6};
Circle(6) = {6,7,5};

Line Loop(2) = {5,6};




// Generation d'une ellipse
Point(8) = {4,2,0,h1};
Point(9) = {5,1.5,0,h1};
Point(10) = {5,2,0,h1};
Point(11) = {6,2,0,h1};
Point(12) = {5,2.5,0,h1};

Ellipse(7) = {8,10,8,9};
Ellipse(8) = {9,10,11,11};
Ellipse(9) = {11,10,11,12};
Ellipse(10) = {12,10,8,8};

Line Loop(3) = {7,8,9,10};



// Generation d'un cercle plus petit
Point(13) = {9,3,0,h2};
Point(14) = {9.5,3,0,h2};
Point(15) = {9.25,3,0,h2};

Circle(11) = {13,15,14};
Circle(12) = {14,15,13};

Line Loop(4) = {11,12};



// Generation d'un cercle plus petit
Point(16) = {9,1,0,h2};
Point(17) = {9.5,1,0,h2};
Point(18) = {9.25,1,0,h2};

Circle(13) = {16,18,17};
Circle(14) = {17,18,16};

Line Loop(5) = {13,14};





// Generation d'un rectangle avec un cercle au bout
Point(19) = {8,1.5,0,h1};
Point(20) = {11,1.5,0,h1};
Point(21) = {11,2.5,0,h1};
Point(22) = {8,2.5,0,h1};
Point(23) = {11,2,0,h1};
Point(24) = {8,2,0,h1};

Line(15) = {21,22};
Circle(16) = {22,24,19};
Line(17) = {19,20};
Circle(18) = {20,23,21};

Line Loop(6) = {15,16,17,18};





// Generation d'un rectangle fin
Point(25) = {6.8,3.5,0,h1};
Point(26) = {6.8,0.5,0,h1};
Point(27) = {7.0,0.5,0,h1};
Point(28) = {7.0,3.5,0,h1};

Line(19) = {25,26};
Line(20) = {26,27};
Line(21) = {27,28};
Line(22) = {28,25};


Line Loop(7) = {19,20,21,22};






Plane Surface(1) = {1,2,3,4,5,6,7};


Physical Line(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
Physical Surface(2) = {1};
