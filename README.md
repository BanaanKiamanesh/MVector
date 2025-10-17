# MVector — A fast MATLAB implementation of Processing's PVector (3/4)

**MVector** is a high-performance, mutable 2D/3D vector class for MATLAB, modeled closely after **Processing's PVector** (Processing 3/4). It mirrors the PVector API and behavior wherever possible while embracing MATLAB best practices for speed, readability, and robust error handling.

* ✅ 2D and 3D vectors (auto-promote to 3D when `z` is used)
* ✅ In-place, chainable operations (handle class)
* ✅ No toolboxes required
* ✅ Clear errors & input validation
* ✅ Static operations (PVector-style) with MATLAB-safe names

> **Attribution**: This is a MATLAB implementation inspired by **Processing's PVector** class (Processing 3/4). It is not affiliated with or endorsed by the Processing Foundation.

---

## Table of Contents

* [Folder Layout](#folder-layout)
* [Requirements](#requirements)
* [Installation](#installation)

  * [Temporary (current session)](#temporary-current-session)
  * [Permanent](#permanent)
  * [Reset (remove paths)](#reset-remove-paths)
* [Quickstart](#quickstart)
* [API Overview](#api-overview)

  * [Constructors](#constructors)
  * [Properties](#properties)
  * [Instance Methods](#instance-methods)
  * [Static Methods](#static-methods)
  * [Compatibility Notes](#compatibility-notes)
* [Examples](#examples)
* [Testing](#testing)
* [Performance Tips](#performance-tips)
* [Troubleshooting](#troubleshooting)
* [License](#license)
* [Acknowledgements](#acknowledgements)

---

## Folder Layout

```
<repo-root>/
  install.m        % installer (add src/test to MATLAB path)
  src/
    MVector.m      % class definition
  test/
    Test.m         % unit test script covering functionality
  README.md
```

---

## Requirements

* MATLAB R2016b or newer (tested on 2025b releases)
* No additional toolboxes

---

## Installation

Place your MATLAB current directory at the **repo root** (where `install.m` lives).

### Temporary (current session)

```matlab
>> install
```

Adds:

* `<root>` (so you can still call `install`)
* `<root>/src` (recursively)
* `<root>/test` (recursively)

### Permanent

```matlab
>> install('save')
```

Does the same as temporary install, then attempts to **save the path** so `MVector` is available in future sessions.

> If MATLAB cannot save the path (permissions), run MATLAB as administrator / with elevated privileges, or use the **Set Path** UI to save.

### Reset (remove paths)

```matlab
>> install('reset')
```

Removes only the `/src` and `/test` paths from the current session (keeps the repo root so `install` remains callable).

---

## Quickstart

```matlab
% After running install or install('save'):

v = MVector(1, 2);      % 2D vector [1 2]
v.add(3, 4).normalize();% in-place, chainable
len = v.mag();          % magnitude

w = MVector.fromAngle(pi/4);  % 2D unit vector at 45 degrees
u = MVector.addS(v, w);       % static add -> new MVector

disp(v.toString())
disp(u.array())
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

* `x`, `y`, `z` (dependent properties).
  Setting `z` automatically promotes the vector to 3D.

```matlab
v = MVector(1,2);
v.z = 5;       % now 3D: [1 2 5]
```

### Instance Methods

All instance methods are **in-place** (they modify the object) and **chainable**.

* `copy()` → returns a new `MVector` copy
* `get(template)` → returns `[x y]` or `[x y z]`. If `template` is a length-2 or length-3 vector, the output matches its length and shape (row/column).
* `set(x,y)` / `set(x,y,z)` / `set(MVector)` / `set([..])`
* `add(MVector|[..])` / `add(x,y)` / `add(x,y,z)`
* `sub(MVector|[..])` / `sub(x,y)` / `sub(x,y,z)`
* `mult(n)` / `div(n)` (scalar; `div` requires nonzero)
* `mag()` / `magSq()` (magnitude and squared magnitude)
* `dist(MVector)` / `dot(MVector|components)` / `cross(MVector[, target])`
* `normalize([target])` — if a `target` is provided, writes normalized values there
* `limit(max)` — clamps magnitude
* `setMag(len)` / `setMag(target, len)`
* `heading()` — 2D heading angle `atan2(y,x)` (z ignored)
* `rotate(theta)` — 2D rotation by `theta` (radians), preserves `z`
* `lerp(v, amt)` / `lerp(x,y,amt)` / `lerp(x,y,z,amt)` — in place
* `array()` — `[x y]` or `[x y z]`
* `toString()` — pretty string
* `equals(obj)` — exact equality (see note below)

> **Note**: `hashCode()` is intentionally **not implemented** (omitted by request).

### Static Methods

To avoid MATLAB name clashes (MATLAB disallows static and instance methods with the same name), PVector static ops are provided with a **`*S` suffix**:

* `addS(v1, v2)` / `addS(v1, v2, target)`
* `subS(v1, v2)` / `subS(v1, v2, target)`
* `multS(v, n)` / `multS(v, n, target)`
* `divS(v, n)` / `divS(v, n, target)`
* `distS(v1, v2)`
* `dotS(v1, v2)`
* `crossS(v1, v2, target)` (**target required**, returns 3D cross)
* `lerpS(v1, v2, amt)`
* `random2D([target[, parent]])` — returns/writes a 2D unit vector (parent ignored)
* `random3D([target[, parent]])` — returns/writes a 3D unit vector (parent ignored)
* `fromAngle(angle[, target])` — 2D unit vector
* `angleBetween(v1, v2)` — angle in radians

---

## Compatibility Notes

* **Static names**: Because MATLAB cannot have both instance and static methods with the same name in one class, static PVector methods are
  exposed as `addS`, `subS`, `multS`, `divS`, `dotS`, `distS`, `crossS`, `lerpS`.

  > If you need *exact* PVector names, create a tiny wrapper class (e.g., `MVectorOps`) that forwards to these `*S` methods.

* **`get(target)`**: Processing's `get(float[] target)` is emulated as `get(template)` where `template`'s **length and shape** (2 or 3, row/column) determines the returned shape.

* **`cross`**:

  * Instance: `cross(v)` returns a **new** 3D vector; `cross(v, target)` writes into `target`.
  * Static: `crossS(v1, v2, target)` requires `target` (to avoid ambiguity in MATLAB).

* **`equals`**: **Exact** component equality. For floating-point comparisons with tolerance, use your own check (e.g., `norm(a.array()-b.array()) < tol`).

* **`heading` / `rotate`**: 2D only (z is ignored/preserved), matching typical PVector usage.

* **Random**: `random2D`/`random3D` accept an optional `parent` arg (ignored) for Processing signature compatibility.

* **`hashCode`**: intentionally **not implemented**.

---

## Examples

### Basic algebra (chainable)

```matlab
v = MVector(1, 2).add(3, 4).mult(0.5);  % => [2 3]
u = MVector(0, 0, 5).sub(0, 0, 2);      % => [0 0 3]
```

### Normalization, limiting, setting magnitude

```matlab
v = MVector(3, 4).normalize();           % => [0.6 0.8]
v.limit(0.5);                            % magnitude <= 0.5

src = MVector(3, 4);
dst = MVector();
src.setMag(dst, 10);                     % writes normalized*10 into dst
```

### Dot, cross, distance

```matlab
a = MVector(1, 2, 3);
b = MVector(4, 5, 6);
sd = a.dot(b);                           % 32
tx = MVector();
MVector.crossS(a, b, tx);                % tx = cross(a,b)

d2 = MVector(0,0).dist(MVector(3,4));    % 5
```

### Lerp & angles

```matlab
p = MVector(0,0).lerp(10, 0, 0.25);      % => [2.5 0]
q = MVector.fromAngle(pi/3);             % [cos, sin]
ang = MVector.angleBetween(MVector(1,0), MVector(0,1));  % pi/2
```

### Random unit vectors

```matlab
r2 = MVector.random2D();                 % 2D unit vector
r3 = MVector.random3D();                 % 3D unit vector

t = MVector();
MVector.random2D(t);                     % write into target
```

---

## Testing

From the repo root (after `install`):

```matlab
>> Test
Running MVector tests...
All MVector tests passed.
```

`test/Test.m` covers constructors, all instance methods (except `hashCode()` by design), and static methods. It also validates random 2D/3D outputs are unit length.

---

## Performance Tips

* **Prefer in-place ops**: The class is a handle; operations like `add`, `sub`, `mult`, etc. **mutate** and are **chainable**.
* **Avoid creating temporaries** in tight loops. Reuse target vectors for static ops (`crossS`, `addS`, etc.) to minimize allocations.
* **Use `magSq()`** when you only need relative magnitudes—no `sqrt`, faster in hot paths.
* **Stay in 2D when possible**: z = 0 vectors remain 2D; promoting to 3D adds minor overhead.

---

## Troubleshooting

* **`MVector` not found**
  Run `install` from the repo root. If you moved files, verify `src/MVector.m` exists and is readable.

* **Path not saved**
  Use `install('save')`. If it warns, run MATLAB as administrator/elevated or save the path via **Set Path** > **Save**.

* **Class changes not picked up**
  Run:

  ```matlab
  clear classes; rehash;  % or restart MATLAB
  ```

* **Error IDs require message strings**
  All thrown errors use `error('ID','Message')`. If you add/modify errors, ensure both **ID** and **message** are provided.

* **Floating-point equality**
  `equals()` uses exact equality. For tolerance checks:

  ```matlab
  tol = 1e-12;
  isequal = norm(a.array() - b.array()) <= tol;
  ```

---

## License

Add your license text here (e.g., MIT). If including third-party assets, ensure their licenses are respected.

---

## Acknowledgements

* Inspired by **Processing's PVector** (Processing 3/4).
* Thanks to the Processing community for the original API design that makes vector math friendly and consistent.

---

Happy vectoring!
