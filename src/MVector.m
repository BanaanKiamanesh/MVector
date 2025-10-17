classdef MVector < handle
%MVector Fast 2D/3D mutable vector (Processing PVector-style) for MATLAB.
% 
% MVector is a high-performance, mutable vector with Processing PVector-like
% API. It automatically handles 2D and 3D usage (z=0 for 2D).
%
% Syntax (construction)
%   v = MVector()                 % -> [0 0] (2D)
%   v = MVector(x,y)              % -> [x y] (2D)
%   v = MVector(x,y,z)            % -> [x y z] (3D)
%   v = MVector([x y]) | MVector([x y z])
%   v = MVector(otherMVector)     % copy
%
% Key behavior
%   • Mutable & chainable (handle class): v.add(...).normalize().mult(5)
%   • 2D/3D automatic; setting/using z promotes to 3D
%   • Robust input validation & clear error messages
%   • No toolboxes required
%
% Instance methods (PVector-equivalent)
%   copy, get, set, add, sub, mult, div, mag, magSq, dist, dot, cross,
%   normalize, limit, setMag, heading, rotate, lerp, array, toString, equals
%
% Static methods
%   addS, subS, multS, divS, distS, dotS, crossS, lerpS   % static ops (*S naming)
%   random2D, random3D, fromAngle, angleBetween
%
% NOTE on static names:
%   MATLAB disallows static and instance methods with the same name. To
%   mirror PVector's static operations (add/sub/mult/div/dot/dist/cross/lerp),
%   this class provides them with a "*S" suffix (e.g., MVector.addS(v1,v2)).
%
% Compatibility notes (Processing PVector)
%   • get(targetArray) is emulated via get(template): returns values in the
%     shape/length of 'template' (length 2 or 3, row or column).
%   • random2D / random3D accept an optional target MVector. Any extra
%     argument (e.g., a PApplet 'parent') is accepted and ignored.
%
% Example
%   v = MVector(1,2).add(3,4).normalize().mult(10);
%   w = MVector.fromAngle(pi/4);
%   a = MVector.addS(v,w);
%   t = MVector(); MVector.subS(a,w,t);

    properties (Access = private)
        v   (1,3) double = [0 0 0];
        dim (1,1) double {mustBeMember(dim,[2 3])} = 2;
    end

    properties (Dependent)
        x; y; z;
    end

    %% Dependent property accessors
    methods
        function val = get.x(obj), val = obj.v(1); end
        function set.x(obj,val), obj.v(1) = obj.assertFiniteScalar(val,'x must be finite.'); end
        function val = get.y(obj), val = obj.v(2); end
        function set.y(obj,val), obj.v(2) = obj.assertFiniteScalar(val,'y must be finite.'); end
        function val = get.z(obj), val = obj.v(3); end
        function set.z(obj,val), obj.v(3) = obj.assertFiniteScalar(val,'z must be finite.'); obj.dim=3; end
    end

    %% Constructors
    methods
        function obj = MVector(varargin)
            if nargin==0
                obj.v=[0 0 0]; obj.dim=2; 
                return
            end
            if nargin==1
                a = varargin{1};
                if isa(a,'MVector')
                    obj.v=a.v; obj.dim=a.dim; 
                    return
                end
                a = obj.ensureRowVec(a,'Constructor input must be 1x2 or 1x3 numeric or MVector.');
                n=numel(a);
                if n==2
                    obj.v=[a(1) a(2) 0]; obj.dim=2;
                elseif n==3
                    obj.v=[a(1) a(2) a(3)]; obj.dim=3;
                else
                    error('MVector:Constructor:BadLength','Length must be 2 or 3.');
                end
                return
            end
            if nargin==2 || nargin==3
                nums = cellfun(@(x) obj.assertFiniteScalar(x,'Constructor components must be finite.'), varargin);
                if nargin==2
                    obj.v=[nums(1) nums(2) 0]; obj.dim=2;
                else
                    obj.v=[nums(1) nums(2) nums(3)]; obj.dim=3;
                end
                return
            end
            error('MVector:Constructor:Arity','Accepts 0, 1, 2 or 3 args.');
        end
    end

    %% Instance methods
    methods
        function out = copy(obj)
        %COPY  Deep copy of this vector.
            out = MVector(obj);
        end

        function out = get(obj, template)
        %GET  Return coordinates; if TEMPLATE given (length 2 or 3), return
        %      values matching its shape (row/col) and length.
            if nargin==1
                out = obj.v(1:obj.dim); 
                return
            end
            if ~isnumeric(template) || ~isvector(template)
                error('MVector:get:TemplateType','Template must be numeric vector.');
            end
            len = numel(template);
            if len~=2 && len~=3
                error('MVector:get:TemplateLen','Template length must be 2 or 3.');
            end
            vals = obj.v(1:min(obj.dim,len));
            if len==3 && obj.dim==2
                vals = [vals 0];
            end
            if iscolumn(template)
                out = vals(:);
            else
                out = vals;
            end
        end

        function obj = set(obj, varargin)
        %SET  Overloaded: set(x,y) | set(x,y,z) | set(MVector) | set([..]).
            n = nargin-1;
            switch n
                case 1
                    a = varargin{1};
                    if isa(a,'MVector')
                        obj.v=a.v; obj.dim=a.dim; 
                        return
                    end
                    a = obj.ensureRowVec(a,'set(array) expects 1x2 or 1x3 numeric.');
                    obj.applySetFromArray(a);
                case 2
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    obj.v=[x y 0]; obj.dim=2;
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    obj.v=[x y z]; obj.dim=3;
                otherwise
                    error('MVector:set:Arity','set accepts 1,2,3 args.');
            end
        end

        function obj = add(obj, varargin)
        %ADD  In-place add: add(MVector|[...]) | add(x,y) | add(x,y,z).
            n=nargin-1;
            switch n
                case 1
                    a=varargin{1};
                    if isa(a,'MVector')
                        if a.dim==3 || obj.dim==3
                            obj.dim=3; obj.v=obj.v+[a.v(1) a.v(2) a.v(3)];
                        else
                            obj.v(1:2)=obj.v(1:2)+a.v(1:2);
                        end
                    else
                        a=obj.ensureRowVec(a,'add(array) expects 1x2 or 1x3 numeric.');
                        obj.addFromArray(a);
                    end
                case 2
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    obj.v(1)=obj.v(1)+x; obj.v(2)=obj.v(2)+y;
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    obj.v=obj.v+[x y z]; obj.dim=3;
                otherwise
                    error('MVector:add:Arity','add accepts 1,2,3 args.');
            end
        end

        function obj = sub(obj, varargin)
        %SUB  In-place subtract: sub(MVector|[...]) | sub(x,y) | sub(x,y,z).
            n=nargin-1;
            switch n
                case 1
                    a=varargin{1};
                    if isa(a,'MVector')
                        if a.dim==3 || obj.dim==3
                            obj.dim=3; obj.v=obj.v-[a.v(1) a.v(2) a.v(3)];
                        else
                            obj.v(1:2)=obj.v(1:2)-a.v(1:2);
                        end
                    else
                        a=obj.ensureRowVec(a,'sub(array) expects 1x2 or 1x3 numeric.');
                        obj.subFromArray(a);
                    end
                case 2
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    obj.v(1)=obj.v(1)-x; obj.v(2)=obj.v(2)-y;
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    obj.v=obj.v-[x y z]; obj.dim=3;
                otherwise
                    error('MVector:sub:Arity','sub accepts 1,2,3 args.');
            end
        end

        function obj = mult(obj, n)
        %MULT  In-place multiply by scalar n.
            n=obj.assertFiniteScalar(n,'n must be finite.');
            if obj.dim==3, obj.v=obj.v.*n; else, obj.v(1:2)=obj.v(1:2).*n; end
        end

        function obj = div(obj, n)
        %DIV  In-place divide by nonzero scalar n.
            n=obj.assertFiniteNonzeroScalar(n,'n must be nonzero finite.');
            invn=1.0/n;
            if obj.dim==3, obj.v=obj.v.*invn; else, obj.v(1:2)=obj.v(1:2).*invn; end
        end

        function m = mag(obj)
        %MAG  Euclidean length (2D or 3D as applicable).
            if obj.dim==3, m = hypot(obj.v(1),hypot(obj.v(2),obj.v(3)));
            else,           m = hypot(obj.v(1),obj.v(2)); end
        end

        function m2 = magSq(obj)
        %MAGSQ  Squared magnitude (fast, no sqrt).
            vv=obj.v;
            if obj.dim==3, m2 = vv(1)^2 + vv(2)^2 + vv(3)^2;
            else,          m2 = vv(1)^2 + vv(2)^2; end
        end

        function d = dist(obj, other)
        %DIST  Distance to another MVector.
            if ~isa(other,'MVector'), error('MVector:dist:Type','dist expects an MVector.'); end
            if obj.dim==3 || other.dim==3
                dv=obj.v-[other.v(1) other.v(2) other.v(3)];
                d=hypot(dv(1),hypot(dv(2),dv(3)));
            else
                dv=obj.v(1:2)-other.v(1:2); d=hypot(dv(1),dv(2));
            end
        end

        function s = dot(obj, varargin)
        %DOT  Dot product with MVector or components: dot(v) | dot(x,y[,z]).
            n=nargin-1;
            switch n
                case 1
                    a=varargin{1};
                    if isa(a,'MVector')
                        if obj.dim==3 || a.dim==3
                            s = obj.v(1)*a.v(1)+obj.v(2)*a.v(2)+obj.v(3)*a.v(3);
                        else
                            s = obj.v(1)*a.v(1)+obj.v(2)*a.v(2);
                        end
                    else
                        a=obj.ensureRowVec(a,'dot(array) expects 1x2 or 1x3 numeric.');
                        if numel(a)==2
                            s=obj.v(1)*a(1)+obj.v(2)*a(2);
                        elseif numel(a)==3
                            s=obj.v(1)*a(1)+obj.v(2)*a(2)+obj.v(3)*a(3);
                        else
                            error('MVector:dot:ArrayLen','dot(array) length must be 2 or 3.');
                        end
                    end
                case 2
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    s=obj.v(1)*x+obj.v(2)*y;
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    s=obj.v(1)*x+obj.v(2)*y+obj.v(3)*z;
                otherwise
                    error('MVector:dot:Arity','dot accepts 1,2,3 args.');
            end
        end

        function out = cross(obj, other, target)
        %CROSS  3D cross product; returns new vector or writes into TARGET.
            if ~isa(other,'MVector'), error('MVector:cross:Type','cross expects an MVector.'); end
            ax=obj.v(1); ay=obj.v(2); az=obj.v(3);
            bx=other.v(1); by=other.v(2); bz=other.v(3);
            cx=ay*bz-az*by; cy=az*bx-ax*bz; cz=ax*by-ay*bx;
            if nargin<3
                out=MVector([cx cy cz]);
            else
                if ~isa(target,'MVector'), error('MVector:cross:TargetType','target must be an MVector.'); end
                target.v=[cx cy cz]; target.dim=3; out=target;
            end
        end

        function out = normalize(obj, varargin)
        %NORMALIZE  Make unit length (in-place); or write into TARGET.
            if nargin==1
                m=obj.mag();
                if m>0
                    invm=1.0/m;
                    if obj.dim==3, obj.v=obj.v.*invm; else, obj.v(1:2)=obj.v(1:2).*invm; end
                end
                out=obj;
            elseif nargin==2
                target=varargin{1}; 
                if ~isa(target,'MVector'), error('MVector:normalize:TargetType','target must be an MVector.'); end
                m=obj.mag();
                if m>0
                    invm=1.0/m;
                    if obj.dim==3, target.v=obj.v.*invm; target.dim=3;
                    else, target.v=[obj.v(1)*invm obj.v(2)*invm 0]; target.dim=2; end
                else
                    target.v=[0 0 0]; target.dim=obj.dim;
                end
                out=target;
            else
                error('MVector:normalize:Arity','normalize accepts 0 or 1 argument (target).');
            end
        end

        function obj = limit(obj, maxMag)
        %LIMIT  Clamp magnitude to maxMag (in-place).
            maxMag=obj.assertFiniteNonnegativeScalar(maxMag,'max must be >=0 finite.');
            if maxMag==0
                obj.v(:)=0; 
                return
            end
            if obj.magSq()>maxMag*maxMag
                obj.normalize();
                obj.mult(maxMag);
            end
        end

        function out = setMag(obj, varargin)
        %SETMAG  Set magnitude: setMag(len) in-place or setMag(target,len).
            if nargin==2
                len=obj.assertFiniteNonnegativeScalar(varargin{1},'len must be >=0 finite.');
                m=obj.mag(); if m>0, obj.mult(len/m); end
                out=obj;
            elseif nargin==3
                target=varargin{1};
                len=obj.assertFiniteNonnegativeScalar(varargin{2},'len must be >=0 finite.');
                if ~isa(target,'MVector'), error('MVector:setMag:TargetType','target must be an MVector.'); end
                m=obj.mag();
                if m>0
                    sc=len/m;
                    if obj.dim==3, target.v=obj.v.*sc; target.dim=3;
                    else, target.v=[obj.v(1)*sc obj.v(2)*sc 0]; target.dim=2; end
                else
                    target.v=[0 0 0]; target.dim=obj.dim;
                end
                out=target;
            else
                error('MVector:setMag:Arity','setMag accepts (len) or (target,len).');
            end
        end

        function a = heading(obj)
        %HEADING  2D heading angle (atan2(y,x)); z ignored.
            a = atan2(obj.v(2), obj.v(1));
        end

        function obj = rotate(obj, theta)
        %ROTATE  2D rotation by theta (radians). z preserved.
            theta=obj.assertFiniteScalar(theta,'theta must be finite.');
            c=cos(theta); s=sin(theta); x=obj.v(1); y=obj.v(2);
            obj.v(1)=x*c-y*s; obj.v(2)=x*s+y*c;
        end

        function obj = lerp(obj, varargin)
        %LERP  In-place linear interpolation:
        %      lerp(v,amt) | lerp(x,y,amt) | lerp(x,y,z,amt)
            n=nargin-1;
            switch n
                case 2
                    v2=varargin{1};
                    if ~isa(v2,'MVector'), error('MVector:lerp:Type','lerp(v,amt) expects v as MVector.'); end
                    amt=obj.assertFiniteScalar(varargin{2},'amt must be finite.');
                    if obj.dim==3 || v2.dim==3
                        obj.dim=3; obj.v=obj.v + amt.*([v2.v(1) v2.v(2) v2.v(3)]-obj.v);
                    else
                        obj.v(1:2)=obj.v(1:2)+amt.*(v2.v(1:2)-obj.v(1:2));
                    end
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    amt=obj.assertFiniteScalar(varargin{3},'amt must be finite.');
                    obj.v(1:2)=obj.v(1:2)+amt.*([x y]-obj.v(1:2));
                case 4
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    amt=obj.assertFiniteScalar(varargin{4},'amt must be finite.');
                    obj.v=obj.v+amt.*([x y z]-obj.v); obj.dim=3;
                otherwise
                    error('MVector:lerp:Arity','lerp accepts (v,amt) or (x,y,amt) or (x,y,z,amt).');
            end
        end

        function a = array(obj)
        %ARRAY  Return [x y] or [x y z] as numeric row vector.
            a = obj.v(1:obj.dim);
        end

        function s = toString(obj)
        %TOSTRING  Human-readable string of the vector.
            if obj.dim==3
                s = sprintf('MVector[%.15g, %.15g, %.15g]',obj.v(1),obj.v(2),obj.v(3));
            else
                s = sprintf('MVector[%.15g, %.15g]',obj.v(1),obj.v(2));
            end
        end

        function tf = equals(obj, other)
        %EQUALS  True if same size and same components (exact equality).
            if isa(other,'MVector')
                tf = (obj.dim==other.dim) && all(obj.v(1:obj.dim)==other.v(1:other.dim));
            elseif isnumeric(other) && isvector(other)
                other = obj.ensureRowVec(other,'equals: numeric must be 1x2 or 1x3.');
                n=numel(other);
                if n~=2 && n~=3, error('MVector:equals:ArrayLen','equals(array) length must be 2 or 3.'); end
                tf = (obj.dim==n) && all(obj.v(1:n)==other(1:n));
            else
                tf=false;
            end
        end
    end

    %% Static methods (*S to avoid name clash with instance names)
    methods (Static)
        function out = addS(v1, v2, target)
        %ADDS  Static add: out = addS(v1,v2) or addS(v1,v2,target).
            [a,b]=MVector.assertTwoVectors(v1,v2,'addS');
            if nargin<3
                if a.dim==3||b.dim==3
                    out=MVector([a.v(1)+b.v(1), a.v(2)+b.v(2), a.v(3)+b.v(3)]);
                else
                    out=MVector([a.v(1)+b.v(1), a.v(2)+b.v(2)]);
                end
            else
                MVector.assertTarget(target,'addS');
                if a.dim==3||b.dim==3
                    target.v=[a.v(1)+b.v(1) a.v(2)+b.v(2) a.v(3)+b.v(3)]; target.dim=3;
                else
                    target.v=[a.v(1)+b.v(1) a.v(2)+b.v(2) 0]; target.dim=2;
                end
                out=target;
            end
        end

        function out = subS(v1, v2, target)
        %SUBS  Static subtract: out = subS(v1,v2) or subS(v1,v2,target).
            [a,b]=MVector.assertTwoVectors(v1,v2,'subS');
            if nargin<3
                if a.dim==3||b.dim==3
                    out=MVector([a.v(1)-b.v(1), a.v(2)-b.v(2), a.v(3)-b.v(3)]);
                else
                    out=MVector([a.v(1)-b.v(1), a.v(2)-b.v(2)]);
                end
            else
                MVector.assertTarget(target,'subS');
                if a.dim==3||b.dim==3
                    target.v=[a.v(1)-b.v(1) a.v(2)-b.v(2) a.v(3)-b.v(3)]; target.dim=3;
                else
                    target.v=[a.v(1)-b.v(1) a.v(2)-b.v(2) 0]; target.dim=2;
                end
                out=target;
            end
        end

        function out = multS(v, n, target)
        %MULTS  Static multiply: out = multS(v,n) or multS(v,n,target).
            MVector.assertVector(v,'multS');
            n=MVector.assertFiniteScalarStatic(n,'multS: n must be finite.');
            if nargin<3
                out=MVector(v); out.mult(n);
            else
                MVector.assertTarget(target,'multS');
                if v.dim==3, target.v=v.v.*n; target.dim=3;
                else, target.v=[v.v(1)*n v.v(2)*n 0]; target.dim=2; end
                out=target;
            end
        end

        function out = divS(v, n, target)
        %DIVS  Static divide: out = divS(v,n) or divS(v,n,target).
            MVector.assertVector(v,'divS');
            n=MVector.assertFiniteNonzeroScalarStatic(n,'divS: n must be nonzero finite.');
            invn=1.0/n;
            if nargin<3
                out=MVector(v); out.mult(invn);
            else
                MVector.assertTarget(target,'divS');
                if v.dim==3, target.v=v.v.*invn; target.dim=3;
                else, target.v=[v.v(1)*invn v.v(2)*invn 0]; target.dim=2; end
                out=target;
            end
        end

        function d = distS(v1, v2)
        %DISTS  Static distance between two vectors.
            [a,b]=MVector.assertTwoVectors(v1,v2,'distS');
            if a.dim==3||b.dim==3
                dv=[a.v(1)-b.v(1) a.v(2)-b.v(2) a.v(3)-b.v(3)];
                d=hypot(dv(1),hypot(dv(2),dv(3)));
            else
                dv=[a.v(1)-b.v(1) a.v(2)-b.v(2)];
                d=hypot(dv(1),dv(2));
            end
        end

        function s = dotS(v1, v2)
        %DOTS  Static dot product of two vectors.
            [a,b]=MVector.assertTwoVectors(v1,v2,'dotS');
            if a.dim==3||b.dim==3
                s = a.v(1)*b.v(1)+a.v(2)*b.v(2)+a.v(3)*b.v(3);
            else
                s = a.v(1)*b.v(1)+a.v(2)*b.v(2);
            end
        end

        function out = crossS(v1, v2, target)
        %CROSSS  Static cross product into TARGET (required).
            [a,b]=MVector.assertTwoVectors(v1,v2,'crossS'); 
            MVector.assertTarget(target,'crossS');
            ax=a.v(1); ay=a.v(2); az=a.v(3);
            bx=b.v(1); by=b.v(2); bz=b.v(3);
            target.v=[ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx]; target.dim=3; 
            out=target;
        end

        function out = lerpS(v1, v2, amt)
        %LERPS  Static linear interpolation: out = (1-amt)*v1 + amt*v2.
            [a,b]=MVector.assertTwoVectors(v1,v2,'lerpS');
            amt=MVector.assertFiniteScalarStatic(amt,'lerpS: amt must be finite.');
            if a.dim==3||b.dim==3
                out=MVector([a.v(1)+amt*(b.v(1)-a.v(1)), ...
                             a.v(2)+amt*(b.v(2)-a.v(2)), ...
                             a.v(3)+amt*(b.v(3)-a.v(3))]);
            else
                out=MVector([a.v(1)+amt*(b.v(1)-a.v(1)), ...
                             a.v(2)+amt*(b.v(2)-a.v(2))]);
            end
        end

        function out = random2D(varargin)
        %RANDOM2D  Unit vector in random 2D direction.
        %   out = random2D()                 % returns new MVector
        %   out = random2D(target)           % writes into target MVector
        %   out = random2D(target,parent)    % 'parent' ignored (compat)
        %   out = random2D(parent)           % parent ignored (compat)
            target=[];
            if nargin>=1 && isa(varargin{1},'MVector'), target=varargin{1}; end
            theta=2*pi*rand(); c=cos(theta); s=sin(theta);
            if isempty(target), out=MVector([c s]);
            else
                MVector.assertTarget(target,'random2D');
                target.v=[c s 0]; target.dim=2; out=target;
            end
        end

        function out = random3D(varargin)
        %RANDOM3D  Unit vector in random 3D direction.
        %   out = random3D()                 % returns new MVector
        %   out = random3D(target)           % writes into target MVector
        %   out = random3D(target,parent)    % 'parent' ignored (compat)
        %   out = random3D(parent)           % parent ignored (compat)
            target=[];
            if nargin>=1 && isa(varargin{1},'MVector'), target=varargin{1}; end
            z=2*rand()-1; t=2*pi*rand(); r=sqrt(max(0,1-z*z)); x=r*cos(t); y=r*sin(t);
            if isempty(target), out=MVector([x y z]);
            else
                MVector.assertTarget(target,'random3D');
                target.v=[x y z]; target.dim=3; out=target;
            end
        end

        function out = fromAngle(angle, target)
        %FROMANGLE  2D unit vector from angle (radians).
            a=MVector.assertFiniteScalarStatic(angle,'fromAngle: angle must be finite.');
            if nargin<2
                out=MVector([cos(a) sin(a)]);
            else
                MVector.assertTarget(target,'fromAngle');
                target.v=[cos(a) sin(a) 0]; target.dim=2; out=target;
            end
        end

        function ang = angleBetween(v1, v2)
        %ANGLEBETWEEN  Angle between two vectors (radians).
            [a,b]=MVector.assertTwoVectors(v1,v2,'angleBetween');
            ma=a.mag(); mb=b.mag();
            if ma==0 || mb==0, ang=0.0; return, end
            d=MVector.dotS(a,b)/(ma*mb);
            d=max(-1.0,min(1.0,d));
            ang=acos(d);
        end
    end

    %% Private helpers
    methods (Access = private)
        function applySetFromArray(obj, a)
            n=numel(a);
            if n==2
                obj.v=[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite') 0]; obj.dim=2;
            elseif n==3
                obj.v=[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite') obj.assertFiniteScalar(a(3),'z finite')]; obj.dim=3;
            else
                error('MVector:set:ArrayLen','set(array) length must be 2 or 3.');
            end
        end
        function addFromArray(obj,a)
            n=numel(a);
            if n==2
                obj.v(1:2)=obj.v(1:2)+[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite')];
            elseif n==3
                obj.v=obj.v+[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite') obj.assertFiniteScalar(a(3),'z finite')]; obj.dim=3;
            else
                error('MVector:add:ArrayLen','add(array) length must be 2 or 3.');
            end
        end
        function subFromArray(obj,a)
            n=numel(a);
            if n==2
                obj.v(1:2)=obj.v(1:2)-[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite')];
            elseif n==3
                obj.v=obj.v-[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite') obj.assertFiniteScalar(a(3),'z finite')]; obj.dim=3;
            else
                error('MVector:sub:ArrayLen','sub(array) length must be 2 or 3.');
            end
        end
        function v = ensureRowVec(~, a, msg)
            if ~isnumeric(a) || ~isvector(a), error('MVector:VectorType',msg); end
            if iscolumn(a), v = a.'; else, v = a; end
        end
        function x = assertFiniteScalar(~, x, msg)
            if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x), error('MVector:FiniteScalar',msg); end
            x=double(x);
        end
        function x = assertFiniteNonzeroScalar(~, x, msg)
            if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x) || x==0, error('MVector:NonzeroScalar',msg); end
            x=double(x);
        end
        function x = assertFiniteNonnegativeScalar(~, x, msg)
            if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x) || x<0, error('MVector:NonnegativeScalar',msg); end
            x=double(x);
        end
    end

    %% Static private helpers
    methods (Static, Access = private)
        function [a,b]=assertTwoVectors(v1,v2,fname)
            if ~isa(v1,'MVector') || ~isa(v2,'MVector')
                error(sprintf('MVector:%s:Type',fname),'Both inputs must be MVector.');
            end
            a=v1; b=v2;
        end
        function assertVector(v,fname)
            if ~isa(v,'MVector')
                error(sprintf('MVector:%s:Type',fname),'Input must be an MVector.');
            end
        end
        function assertTarget(target,fname)
            if ~isa(target,'MVector')
                error(sprintf('MVector:%s:TargetType',fname),'target must be an MVector.');
            end
        end
        function x=assertFiniteScalarStatic(x,msg)
            if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x), error('MVector:FiniteScalar',msg); end
            x=double(x);
        end
        function x=assertFiniteNonzeroScalarStatic(x,msg)
            if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x) || x==0, error('MVector:NonzeroScalar',msg); end
            x=double(x);
        end
    end
end
