clear
close all
clc
clear classes
rehash

fprintf('Running MVector tests...\n');
tol = 1e-12;

%% Constructors
v0 = MVector();                 assertEqual(v0.get(), [0 0], tol);
v2 = MVector(1,2);              assertEqual(v2.get(), [1 2], tol);
v3 = MVector(1,2,3);            assertEqual(v3.get(), [1 2 3], tol);
v4 = MVector([4 5]);            assertEqual(v4.get(), [4 5], tol);
v5 = MVector([4 5 6]);          assertEqual(v5.get(), [4 5 6], tol);
v6 = MVector(v5);               assertEqual(v6.get(), [4 5 6], tol);

%% get with template shape
row3 = v2.get([0 0 0]);         assertEqual(row3, [1 2 0], tol);
col3 = v2.get(zeros(3,1));      assert(all(size(col3)==[3 1]));  assertEqual(col3.', [1 2 0], tol);

%% set: functional vs in-place
u0 = MVector(9,9);
u1 = u0.set(7,8);               % functional: u0 unchanged
assertEqual(u0.get(), [9 9], tol);
assertEqual(u1.get(), [7 8], tol);

u1.set(1,2,3);                  % in-place
assertEqual(u1.get(), [1 2 3], tol);
u1.set([9 10]);                 % in-place, demote to 2D
assertEqual(u1.get(), [9 10], tol);
u1.set(MVector([2 3 4]));       % in-place, promote to 3D
assertEqual(u1.get(), [2 3 4], tol);

%% add / sub: functional vs in-place
v = MVector(1,2);
w = v.add(5,6);                 % functional
assertEqual(v.get(), [1 2], tol);
assertEqual(w.get(), [6 8], tol);

v.add([1 1 1]);                 % in-place (promote)
assertEqual(v.get([0 0 0]), [2 3 1], tol);
v2 = MVector(1,1,1);
v.add(v2);                      % in-place
assertEqual(v.get(), [3 4 2], tol);

v = MVector(10,10);
x = v.sub(1,2);                 % functional
assertEqual(v.get(), [10 10], tol);
assertEqual(x.get(), [9 8], tol);
v.sub([1 1 1]);                 % in-place
assertEqual(v.get([0 0 0]), [9 9 -1], tol);
v.sub(MVector(1,2,3));          % in-place
assertEqual(v.get(), [8 7 -4], tol);

%% mult/div, mag, magSq, dist, dot, cross (instance)
v = MVector(1,2).mult(3);       % functional
assertEqual(v.get(), [3 6], tol);
z = MVector(1,2,3); z.mult(2);  % in-place
assertEqual(z.get(), [2 4 6], tol);

d = MVector(2,4,6).div(2);      % functional
assertEqual(d.get(), [1 2 3], tol);

v = MVector(3,4);               assert(abs(v.mag() - 5) < tol);
assert(abs(v.magSq() - 25) < tol);
v = MVector(1,2,2);             assert(abs(v.mag() - 3) < tol);
assert(abs(v.magSq() - 9) < tol);

a = MVector(0,0); b = MVector(3,4);                 assert(abs(a.dist(b) - 5) < tol);
a3 = MVector(0,0,0); b3 = MVector(1,2,2);           assert(abs(a3.dist(b3) - 3) < tol);

a = MVector(1,2); b = MVector(3,4);                 assert(abs(a.dot(b) - 11) < tol);
a = MVector(1,2,3); b = MVector(4,5,6);             assert(abs(a.dot(b) - 32) < tol);

x = MVector(1,0,0); y = MVector(0,1,0);
z = x.cross(y);                                      assertEqual(z.get(), [0 0 1], tol);
t = MVector(); x.cross(y, t);                        assertEqual(t.get([0 0 0]), [0 0 1], tol);

%% normalize, limit, setMag, heading, rotate, lerp, array, toString, equals
u = MVector(3,4).normalize();                        assert(abs(u.mag() - 1) < 1e-12);
assertEqual(u.get(), [0.6 0.8], 1e-12);

v = MVector(3,4); v.limit(2);                        assert(abs(v.mag()-2) < 1e-12);
v = MVector(1,2,2); v.limit(2);                      assert(abs(v.mag()-2) < 1e-12);
v = MVector(0,0); v.limit(0);                        assertEqual(v.get(), [0 0], tol);

v = MVector(3,4).setMag(10);                         assert(abs(v.mag()-10) < 1e-12);
src = MVector(3,4); tar = MVector();
src.setMag(tar, 5);                                  assert(abs(tar.mag()-5) < 1e-12); assertEqual(src.get(), [3 4], tol);

assert(abs(MVector(1,0).heading() - 0) < tol);
assert(abs(MVector(0,1).heading() - pi/2) < tol);
assert(abs(MVector(-1,0).heading() - pi) < tol);

v = MVector(1,0).rotate(pi/2);                       assertEqual(v.get(zeros(1,2)), [0 1], 1e-12);
v = MVector(1,0,7).rotate(pi/2);  r = v.get([0 0 0]); assert(abs(r(1)) < 1e-12 && abs(r(2)-1) < 1e-12 && abs(r(3)-7) < 1e-12);

v = MVector(0,0).lerp(MVector(10,0), 0.3);           assertEqual(v.get(), [3 0], tol);
v = MVector(0,0,0).lerp(MVector(0,0,10), 0.5);       assertEqual(v.get(), [0 0 5], tol);

v = MVector(7,8); aarr = v.array();                  assertEqual(aarr, [7 8], tol);
s = v.toString();                                     assert(ischar(s) || isstring(s)); assert(contains(string(s),"MVector"));
assert(MVector(1,2).equals(MVector(1,2)));
assert(~MVector(1,2).equals(MVector(1,2,0)));
assert(MVector(1,2,3).equals([1 2 3]));

%% Operator overloads
a = MVector(1,2); b = MVector(3,4);
c = a + b;                    assertEqual(c.get(), [4 6], tol);
d = b - a;                    assertEqual(d.get(), [2 2], tol);
e = -a;                       assertEqual(e.get(), [-1 -2], tol);

f = MVector(1,2,3) .* 2;      assertEqual(f.get(), [2 4 6], tol);
g = 3 .* MVector(1,2);        assertEqual(g.get(), [3 6], tol);
h = MVector(2,4) ./ 2;        assertEqual(h.get(), [1 2], tol);

% equality operators
assert(MVector(1,2) == MVector(1,2));
assert(~(MVector(1,2) == MVector(1,2,0)));
assert(MVector(1,2,3) == [1 2 3]);
assert(MVector(1,2) ~= [1 3]);

%% random2D / random3D / fromAngle / angleBetween (static utilities retained)
r2 = MVector.random2D();         assert(numel(r2.get())==2);  assert(abs(r2.mag()-1) < 1e-10);
rt = MVector(); MVector.random2D(rt);  assert(numel(rt.get())==2 && abs(rt.mag()-1) < 1e-10);

r3 = MVector.random3D();         assert(numel(r3.get([0 0 0]))==3); assert(abs(r3.mag()-1) < 1e-10);
MVector.random3D(rt);            assert(numel(rt.get([0 0 0]))==3 && abs(rt.mag()-1) < 1e-10);

fa = MVector.fromAngle(pi/3);    assertEqual(fa.get(), [cos(pi/3) sin(pi/3)], 1e-12);
MVector.fromAngle(-pi/2, rt);    assertEqual(rt.get(), [0 -1], 1e-12);

ang = MVector.angleBetween(MVector(1,0), MVector(0,1));     assert(abs(ang - pi/2) < 1e-12);
ang = MVector.angleBetween(MVector(1,0,0), MVector(-1,0,0)); assert(abs(ang - pi) < 1e-12);
ang = MVector.angleBetween(MVector(0,0,0), MVector(1,0,0));  assert(abs(ang - 0) < 1e-12);

fprintf('All MVector tests passed.\n');

%% helper
function assertEqual(a,b,tolerance)
    assert(all(abs(a - b) <= tolerance), ...
        sprintf('Expected [%s] but got [%s]', num2str(b), num2str(a)));
end
