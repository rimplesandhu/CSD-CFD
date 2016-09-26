mfac = 0.08;
cr = 0.557784;
tp = 0.8093964;
ct = 0.3681984;
s = 0.762;
t = 0.0075;

// define the points
Point(1) = {0, 0, t, mfac};
Point(2) = {cr, 0, t, mfac};
Point(3) = {tp+ct, s, t, mfac};
Point(4) = {tp, s, t, mfac};
Point(5) = {0, 0, -t, mfac};
Point(6) = {cr, 0, -t, mfac};
Point(7) = {tp+ct, s, -t, mfac};
Point(8) = {tp, s, -t, mfac};

// defines the lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// define the surfaces
Line Loop(1) = {1, 10, -5, -9};
Plane Surface(1) = {1};
Line Loop(2) = {2, 11, -6, -10};
Plane Surface(2) = {2};
Line Loop(3) = {3, 12, -7, -11};
Plane Surface(3) = {3};
Line Loop(4) = {-4, 12, 8, -9};
Plane Surface(4) = {4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {5};
Line Loop(6) = {5, 6, 7, 8};
Plane Surface(6) = {6};

// volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
