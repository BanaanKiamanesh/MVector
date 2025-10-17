# MVector - A fast MATLAB implementation of Processing's PVector

---

**MVector** is a high-performance, mutable 2D/3D vector class for MATLAB, inspired by **Processing's PVector** (Processing 3/4). It mirrors the core PVector behavior while following MATLAB best practices for speed, readability, and clear error handling.

* 2D and 3D vectors (auto-promote to 3D when `z` is used)
* In-place **or** functional style (same method: mutates when no output is captured; returns a new vector when output **is** captured)
* Operator overloads: `+  -  .*  ./  uminus  ==  ~=`
* No toolboxes required
* Clear, identifier-based errors

> **Attribution**: This is a MATLAB implementation inspired by **Processing's PVector** (Processing 3/4). It is not affiliated with or endorsed by the Processing Foundation.

---

## Folder Layout

```
<repo-root>/
install.m        % installer (adds src/test to the MATLAB path)
src/
MVector.m      % class definition
test/
Test.m         % test script
README.md
```

---

## Requirements

* MATLAB R2016b or newer
* No additional toolboxes

---

## Installation

From the **repo root** (where `install.m` lives):

### Temporary (current session)

```matlab
>> install
```

This adds:

* `<root>` (kept so you can still call `install`)
* `<root>/src` (recursively)
* `<root>/test` (recursively)

### Permanent

```matlab
>> install('save')
```

Adds the same paths and attempts to **save** them for future sessions.

> If saving fails (permissions), run MATLAB with elevated privileges or use **Home → Set Path → Save**.

### Reset (remove paths)

```matlab
>> install('reset')
```

Removes `/src` and `/test` from the current session; keeps `<root>` so `install` remains callable.

---

## Quickstart

```matlab
% After install or install('save'):

% Functional (returns new vectors)
u = MVector(1,2) + MVector(3,4);     % [4 6]
w = MVector(1,2,3) .* 2;             % [2 4 6]
d = MVector(2,4) ./ 2;               % [1 2]

% In-place vs. functional (same methods)
v  = MVector(1,2);
v.add(3,4);                          % in-place: v -> [4 6]
z  = v.add(1,1);                     % functional: z = [5 7], v remains [4 6]

% Dot / cross / magnitude
s  = v.dot([1 0]);                   % 5
cx = MVector(1,0,0).cross(MVector(0,1,0));   % [0 0 1]
m  = v.mag();                        % magnitude

% Factory utilities
r2 = MVector.random2D();             % unit vector in 2D
r3 = MVector.random3D();             % unit vector in 3D
a  = MVector.fromAngle(pi/4);        % [cos,sin]
ang = MVector.angleBetween(MVector(1,0), MVector(0,1));  % pi/2
```

---

## API Overview

### Constructors

```matlab
MVector()                 % -> [0 0] (2D)
MVector(x, y)             % -> [x y] (2D)
MVector(x, y, z)          % -> [x y z] (3D)
MVector([x y])            % -> 2D
MVector([x y z])          % -> 3D
MVector(otherMVector)     % copy
```

### Properties

* `x`, `y`, `z` (dependent). Setting `z` promotes to 3D automatically.

```matlab
v = MVector(1,2);
v.z = 5;   % now 3D: [1 2 5]
```

### Operator Overloads (always return a **new** vector)

* `+`  vector addition: `MVector + MVector` or `MVector + [x y (z)]`
* `-`  vector subtraction
* `uminus`  unary minus (`-v`)
* `.*` scalar multiplication: `vector .* scalar` or `scalar .* vector`
* `./` scalar division: `vector ./ scalar` (scalar on right)
* `==`, `~=` exact component equality (dims must match)

> Element-wise `vector .* vector` is **not supported** (use `dot` or `cross` as appropriate).

### Instance Methods

**Mutate in place** if you **don't** capture an output; **return a new vector** if you **do** capture an output.

* Construction / access
* `copy()`
* `get(template)` — respects length (2/3) and row/column shape of `template`
* `set(x,y)` / `set(x,y,z)` / `set(MVector)` / `set([..])`
* Algebra
* `add(...)` / `sub(...)` / `mult(n)` / `div(n)`
* Norm & limits
* `mag()` / `magSq()` (return scalars)
* `normalize([target])`
* `limit(max)`
* `setMag(len)` / `setMag(target, len)`
* Angles & transforms
* `heading()` — `atan2(y,x)` (2D)
* `rotate(theta)` — 2D rotation (preserves `z`)
* `lerp(v,amt)` / `lerp(x,y[,z],amt)`
* Utilities
* `array()` — `[x y]` or `[x y z]`
* `toString()`
* `equals(obj)` — exact equality
* Other scalar returns
* `dist(MVector)` — scalar distance
* `dot(MVector | [x y (z)] | x,y[,z])` — scalar dot product
* `cross(MVector[, target])` — cross product (returns new or writes into `target`)

### Static Utilities

* `random2D([target[, parent]])` — parent arg ignored (kept for signature parity)
* `random3D([target[, parent]])`
* `fromAngle(angle[, target])`
* `angleBetween(v1, v2)`

---

## Examples

### Dual semantics (in-place vs. return)

```matlab
v = MVector(1,2);

% In-place:
v.add(3,4).mult(2);          % v -> ([1 2] + [3 4]) * 2 = [8 12]

% Functional (non-mutating):
w = v.add(1,1);              % w -> [9 13], v remains [8 12]
u = v.mult(0.5);             % u -> [4 6],  v remains [8 12]
```

### Operator style

```matlab
a = MVector(1,2);
b = MVector(3,4);
c = a + b;                   % [4 6]
d = b - a;                   % [2 2]
e = -a;                      % [-1 -2]
f = a .* 3;                  % [3 6]
g = MVector(2,4) ./ 2;       % [1 2]
```

### Geometry

```matlab
p = MVector(3,4);
p.normalize();               % in-place -> unit length
p = p.setMag(10);            % functional -> returns scaled copy

th  = MVector(1,0).heading();            % 0
rot = MVector(1,0).rotate(pi/2);         % [0 1]
ang = MVector.angleBetween(MVector(1,0), MVector(0,1));  % pi/2
```

---

## Testing

From the repo root (after `install`):

```matlab
>> Test
Running MVector tests...
All MVector tests passed.
```

`test/Test.m` exercises constructors, operator overloads, in-place/functional behavior, and all methods/utilities.

---

## Performance Tips

* Use **in-place** methods (`add/sub/mult/div/normalize/limit/setMag/rotate/lerp`) in tight loops to avoid allocations.
* Prefer `magSq()` when only relative magnitudes matter (avoids `sqrt`).
* Stay in 2D when possible (z=0) to avoid unnecessary 3D work.

---

## Troubleshooting

* **`MVector` not found** → run `install` from repo root; verify `src/MVector.m` exists.
* **Path not saved** → `install('save')` with elevated privileges, or use **Set Path → Save**.
* **Class changes not picked up** → `clear classes; rehash;` (or restart MATLAB).
* **Error IDs** → all `error()` calls use identifiers + messages; match this style when extending.
* **Floating-point equality** → `equals()` is exact; for tolerances:

```matlab
tol = 1e-12;
isequal = norm(a.array() - b.array()) <= tol;
```

---

## License

(Add your license text here, e.g., MIT.)

---

## Acknowledgements

* Inspired by **Processing's PVector** (Processing 3/4).
* Thanks to the Processing community for a simple, pragmatic vector API.

---

Happy vectoring!
