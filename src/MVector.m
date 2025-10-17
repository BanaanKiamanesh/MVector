classdef MVector < handle
%MVector Fast 2D/3D mutable vector (Processing PVector-style) for MATLAB.
%
% Key features
%   • 2D/3D vectors (z=0 for 2D; touching z promotes to 3D)
%   • In-place OR functional style:
%       - No output  -> in-place (mutates this)
%       - With output -> returns new MVector (does NOT mutate this)
%   • Operator overloads: +  -  .*  ./  uminus  ==  ~=
%   • No toolboxes required
%
% Constructors
%   v = MVector()                 % -> [0 0] (2D)
%   v = MVector(x,y)              % -> [x y] (2D)
%   v = MVector(x,y,z)            % -> [x y z] (3D)
%   v = MVector([x y]) | MVector([x y z])
%   v = MVector(otherMVector)     % copy
%
% Instance methods (mutate-in-place if no output; otherwise return new)
%   copy            get(template)        set(x,y|x,y,z|MVector|[...])
%   add(...)        sub(...)             mult(n)             div(n)
%   normalize()     limit(max)           setMag(len|target,len)
%   rotate(theta)   lerp(v,amt | x,y[,z],amt)
%
% Query / non-mutating returns (always return outputs)
%   mag()           magSq()              dist(v)             dot(...)
%   cross(v [,target])                  heading()           array()
%   toString()      equals(obj)
%
% Static utilities (kept static for factory-style usage)
%   random2D([target[, parent]])        random3D([target[, parent]])
%   fromAngle(angle[, target])          angleBetween(v1, v2)

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

    %% Operator overloads (always return NEW MVector; non-mutating)
    methods
        function out = plus(a,b)
            % a + b : MVector ⊕ (MVector | [x y] | [x y z])
            [va, da] = MVector.asRowWithDim(a);
            [vb, db] = MVector.asRowWithDim(b);
            dim = max(da, db);
            if dim==2, out = MVector([va(1)+vb(1), va(2)+vb(2)]);
            else,       out = MVector([va(1)+vb(1), va(2)+vb(2), va(3)+vb(3)]); end
        end

        function out = minus(a,b)
            % a - b : MVector ⊖ (MVector | [x y] | [x y z])
            [va, da] = MVector.asRowWithDim(a);
            [vb, db] = MVector.asRowWithDim(b);
            dim = max(da, db);
            if dim==2, out = MVector([va(1)-vb(1), va(2)-vb(2)]);
            else,       out = MVector([va(1)-vb(1), va(2)-vb(2), va(3)-vb(3)]); end
        end

        function out = uminus(a)
            % -a : unary negate
            if a.dim==2, out = MVector([-a.v(1), -a.v(2)]);
            else,         out = MVector([-a.v(1), -a.v(2), -a.v(3)]); end
        end

        function out = times(a,b)
            % a .* b : scalar .* MVector OR MVector .* scalar
            if isa(a,'MVector') && isnumeric(b) && isscalar(b)
                s = double(b); dim=a.dim;
                if dim==2, out = MVector([a.v(1)*s, a.v(2)*s]);
                else,       out = MVector([a.v(1)*s, a.v(2)*s, a.v(3)*s]); end
            elseif isa(b,'MVector') && isnumeric(a) && isscalar(a)
                s = double(a); dim=b.dim;
                if dim==2, out = MVector([b.v(1)*s, b.v(2)*s]);
                else,       out = MVector([b.v(1)*s, b.v(2)*s, b.v(3)*s]); end
            elseif isa(a,'MVector') && isa(b,'MVector')
                error('MVector:times:VectorVector','Element-wise vector.*vector not supported; use dot() or cross().');
            else
                error('MVector:times:Types','times only supports scalar with MVector.');
            end
        end

        function out = rdivide(a,b)
            % a ./ b : MVector ./ scalar   (scalar ./ MVector not supported)
            if isa(a,'MVector') && isnumeric(b) && isscalar(b)
                b = double(b);
                if b==0, error('MVector:rdivide:Zero','Division by zero.'); end
                s = 1.0/b; dim=a.dim;
                if dim==2, out = MVector([a.v(1)*s, a.v(2)*s]);
                else,       out = MVector([a.v(1)*s, a.v(2)*s, a.v(3)*s]); end
            else
                error('MVector:rdivide:Types','Supported form is MVector ./ scalar.');
            end
        end

        function tf = eq(a,b)
            % a == b : exact equality (dims must match)
            [va, da] = MVector.asRowWithDim(a);
            [vb, db] = MVector.asRowWithDim(b);
            if da ~= db, tf = false; return; end
            tf = all(va(1:da) == vb(1:db));
        end

        function tf = ne(a,b)
            tf = ~eq(a,b);
        end
    end

    %% Instance methods (mutate in-place if no output, else return NEW)
    methods
        function out = copy(obj)
            out = MVector(obj);
        end

        function out = get(obj, template)
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
            if iscolumn(template), out = vals(:); else, out = vals; end
        end

        function varargout = set(obj, varargin)
            % set(x,y) | set(x,y,z) | set(MVector) | set([..])
            [nv, nd] = obj.computeSetResult(varargin{:});
            if nargout>0
                varargout{1} = MVector(nv(1:nd));
            else
                obj.v = [nv(1) nv(2) (nd==3)*nv(3)]; obj.dim = nd;
            end
        end

        function varargout = add(obj, varargin)
            [nv, nd] = obj.computeAddResult(varargin{:});
            if nargout>0
                varargout{1} = MVector(nv(1:nd));
            else
                obj.v = [nv(1) nv(2) (nd==3)*nv(3)]; obj.dim = nd;
            end
        end

        function varargout = sub(obj, varargin)
            [nv, nd] = obj.computeSubResult(varargin{:});
            if nargout>0
                varargout{1} = MVector(nv(1:nd));
            else
                obj.v = [nv(1) nv(2) (nd==3)*nv(3)]; obj.dim = nd;
            end
        end

        function varargout = mult(obj, n)
            n=obj.assertFiniteScalar(n,'n must be finite.');
            dim=obj.dim; vv=obj.v;
            if dim==3, nv=vv.*n; else, nv=[vv(1)*n vv(2)*n 0]; end
            if nargout>0
                varargout{1} = MVector(nv(1:dim));
            else
                obj.v = nv; obj.dim = dim;
            end
        end

        function varargout = div(obj, n)
            n=obj.assertFiniteNonzeroScalar(n,'n must be nonzero finite.');
            invn=1.0/n; dim=obj.dim; vv=obj.v;
            if dim==3, nv=vv.*invn; else, nv=[vv(1)*invn vv(2)*invn 0]; end
            if nargout>0
                varargout{1} = MVector(nv(1:dim));
            else
                obj.v = nv; obj.dim = dim;
            end
        end

        function m = mag(obj)
            if obj.dim==3, m = hypot(obj.v(1),hypot(obj.v(2),obj.v(3)));
            else,           m = hypot(obj.v(1),obj.v(2)); end
        end

        function m2 = magSq(obj)
            vv=obj.v;
            if obj.dim==3, m2 = vv(1)^2 + vv(2)^2 + vv(3)^2;
            else,          m2 = vv(1)^2 + vv(2)^2; end
        end

        function d = dist(obj, other)
            if ~isa(other,'MVector'), error('MVector:dist:Type','dist expects an MVector.'); end
            if obj.dim==3 || other.dim==3
                dv=obj.v-[other.v(1) other.v(2) other.v(3)];
                d=hypot(dv(1),hypot(dv(2),dv(3)));
            else
                dv=obj.v(1:2)-other.v(1:2); d=hypot(dv(1),dv(2));
            end
        end

        function s = dot(obj, varargin)
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

        function varargout = normalize(obj, varargin)
            if nargin==1
                m=obj.mag(); vv=obj.v; dim=obj.dim;
                if m>0
                    invm=1.0/m;
                    if dim==3, nv=vv.*invm; else, nv=[vv(1)*invm vv(2)*invm 0]; end
                else
                    nv=[0 0 0];
                end
                if nargout>0, varargout{1}=MVector(nv(1:dim));
                else, obj.v=nv; obj.dim=dim; end
            elseif nargin==2
                target=varargin{1}; 
                if ~isa(target,'MVector'), error('MVector:normalize:TargetType','target must be an MVector.'); end
                m=obj.mag(); dim=obj.dim;
                if m>0
                    invm=1.0/m;
                    if dim==3, target.v=obj.v.*invm; target.dim=3;
                    else, target.v=[obj.v(1)*invm obj.v(2)*invm 0]; target.dim=2; end
                else
                    target.v=[0 0 0]; target.dim=dim;
                end
                varargout{1}=target;
            else
                error('MVector:normalize:Arity','normalize accepts 0 or 1 argument (target).');
            end
        end

        function varargout = limit(obj, maxMag)
            maxMag=obj.assertFiniteNonnegativeScalar(maxMag,'max must be >=0 finite.');
            dim=obj.dim; vv=obj.v;
            if maxMag==0
                nv=[0 0 0];
            else
                if (vv(1)^2+vv(2)^2 + (dim==3)*vv(3)^2) > maxMag*maxMag
                    m = obj.mag(); sc = maxMag/max(m,eps);
                    if dim==3, nv=vv.*sc; else, nv=[vv(1)*sc vv(2)*sc 0]; end
                else
                    nv=vv;
                end
            end
            if nargout>0, varargout{1}=MVector(nv(1:dim));
            else, obj.v=nv; obj.dim=dim; end
        end

        function varargout = setMag(obj, varargin)
            if nargin==2
                len=obj.assertFiniteNonnegativeScalar(varargin{1},'len must be >=0 finite.');
                dim=obj.dim; vv=obj.v; m=obj.mag();
                if m>0, sc=len/m;
                else,    sc=0; end
                if dim==3, nv=vv.*sc; else, nv=[vv(1)*sc vv(2)*sc 0]; end
                if nargout>0, varargout{1}=MVector(nv(1:dim));
                else, obj.v=nv; obj.dim=dim; end
            elseif nargin==3
                target=varargin{1};
                len=obj.assertFiniteNonnegativeScalar(varargin{2},'len must be >=0 finite.');
                if ~isa(target,'MVector'), error('MVector:setMag:TargetType','target must be an MVector.'); end
                dim=obj.dim; m=obj.mag();
                if m>0, sc=len/m; else, sc=0; end
                if dim==3, target.v=obj.v.*sc; target.dim=3;
                else, target.v=[obj.v(1)*sc obj.v(2)*sc 0]; target.dim=2; end
                varargout{1}=target;
            else
                error('MVector:setMag:Arity','setMag accepts (len) or (target,len).');
            end
        end

        function a = heading(obj)
            a = atan2(obj.v(2), obj.v(1));
        end

        function varargout = rotate(obj, theta)
            theta=obj.assertFiniteScalar(theta,'theta must be finite.');
            c=cos(theta); s=sin(theta); x=obj.v(1); y=obj.v(2);
            nv=[x*c-y*s, x*s+y*c, obj.v(3)];
            dim=obj.dim;
            if nargout>0, varargout{1}=MVector(nv(1:dim));
            else, obj.v=nv; obj.dim=dim; end
        end

        function varargout = lerp(obj, varargin)
            n=nargin-1;
            switch n
                case 2
                    v2=varargin{1};
                    if ~isa(v2,'MVector'), error('MVector:lerp:Type','lerp(v,amt) expects v as MVector.'); end
                    amt=obj.assertFiniteScalar(varargin{2},'amt must be finite.');
                    dim=max(obj.dim, v2.dim);
                    if dim==3
                        goal=[v2.v(1) v2.v(2) v2.v(3)];
                        base=[obj.v(1) obj.v(2) obj.v(3)];
                    else
                        goal=v2.v(1:2);
                        base=obj.v(1:2);
                    end
                    res = base + amt*(goal - base);
                    if dim==2, nv=[res(1) res(2) 0]; else, nv=[res(1) res(2) res(3)]; end
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    amt=obj.assertFiniteScalar(varargin{3},'amt must be finite.');
                    base=obj.v(1:2); goal=[x y]; dim=2;
                    res = base + amt*(goal - base);
                    nv=[res(1) res(2) 0];
                case 4
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    amt=obj.assertFiniteScalar(varargin{4},'amt must be finite.');
                    base=[obj.v(1) obj.v(2) obj.v(3)]; goal=[x y z]; dim=3;
                    res = base + amt*(goal - base);
                    nv=[res(1) res(2) res(3)];
                otherwise
                    error('MVector:lerp:Arity','lerp accepts (v,amt) or (x,y,amt) or (x,y,z,amt).');
            end
            if nargout>0, varargout{1}=MVector(nv(1:dim));
            else, obj.v=nv; obj.dim=dim; end
        end

        function a = array(obj)
            a = obj.v(1:obj.dim);
        end

        function s = toString(obj)
            if obj.dim==3
                s = sprintf('MVector[%.15g, %.15g, %.15g]',obj.v(1),obj.v(2),obj.v(3));
            else
                s = sprintf('MVector[%.15g, %.15g]',obj.v(1),obj.v(2));
            end
        end

        function tf = equals(obj, other)
            if isa(other,'MVector')
                tf = (obj.dim==other.dim) && all(obj.v(1:obj.dim)==other.v(1:other.dim));
            elseif isnumeric(other) && isvector(other)
                other = obj.ensureRowVec(other,'equals: numeric must be 1x2 or 1x3.');
                n=numel(other); if n~=2 && n~=3, error('MVector:equals:ArrayLen','equals(array) length must be 2 or 3.'); end
                tf = (obj.dim==n) && all(obj.v(1:n)==other(1:n));
            else
                tf=false;
            end
        end
    end

    %% Static utilities (factories retained)
    methods (Static)
        function out = random2D(varargin)
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
            a=MVector.assertFiniteScalarStatic(angle,'fromAngle: angle must be finite.');
            if nargin<2
                out=MVector([cos(a) sin(a)]);
            else
                MVector.assertTarget(target,'fromAngle');
                target.v=[cos(a) sin(a) 0]; target.dim=2; out=target;
            end
        end

        function ang = angleBetween(v1, v2)
            [a,b]=MVector.assertTwoVectors(v1,v2,'angleBetween');
            ma=a.mag(); mb=b.mag();
            if ma==0 || mb==0, ang=0.0; return, end
            d=(a.v(1)*b.v(1)+a.v(2)*b.v(2)+a.v(3)*b.v(3))/(ma*mb);
            d=max(-1.0,min(1.0,d));
            ang=acos(d);
        end
    end

    %% Private helpers
    methods (Access = private)
        function [nv, nd] = computeSetResult(obj, varargin)
            n = numel(varargin);
            switch n
                case 1
                    a = varargin{1};
                    if isa(a,'MVector')
                        nv=a.v; nd=a.dim;
                    else
                        a = obj.ensureRowVec(a,'set(array) expects 1x2 or 1x3 numeric.');
                        nd=numel(a);
                        if nd==2
                            nv=[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite') 0];
                        elseif nd==3
                            nv=[obj.assertFiniteScalar(a(1),'x finite') obj.assertFiniteScalar(a(2),'y finite') obj.assertFiniteScalar(a(3),'z finite')];
                        else
                            error('MVector:set:ArrayLen','set(array) length must be 2 or 3.');
                        end
                    end
                case 2
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    nv=[x y 0]; nd=2;
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    nv=[x y z]; nd=3;
                otherwise
                    error('MVector:set:Arity','set accepts 1,2,3 args.');
            end
        end

        function [nv, nd] = computeAddResult(obj, varargin)
            % return NEW coords without mutating obj
            base=obj.v; dim=obj.dim;
            n=nargin-1;
            switch n
                case 1
                    a=varargin{1};
                    if isa(a,'MVector')
                        nd = max(dim, a.dim);
                        if nd==3, nv=[base(1)+a.v(1) base(2)+a.v(2) base(3)+a.v(3)];
                        else,      nv=[base(1)+a.v(1) base(2)+a.v(2) 0]; end
                    else
                        a=obj.ensureRowVec(a,'add(array) expects 1x2 or 1x3 numeric.');
                        if numel(a)==2, nd=2; nv=[base(1)+obj.assertFiniteScalar(a(1),'x finite') base(2)+obj.assertFiniteScalar(a(2),'y finite') 0];
                        elseif numel(a)==3, nd=3; nv=[base(1)+obj.assertFiniteScalar(a(1),'x finite') base(2)+obj.assertFiniteScalar(a(2),'y finite') base(3)+obj.assertFiniteScalar(a(3),'z finite')];
                        else, error('MVector:add:ArrayLen','add(array) length must be 2 or 3.');
                        end
                    end
                case 2
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    nd=2; nv=[base(1)+x base(2)+y 0];
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    nd=3; nv=[base(1)+x base(2)+y base(3)+z];
                otherwise
                    error('MVector:add:Arity','add accepts 1,2,3 args.');
            end
        end

        function [nv, nd] = computeSubResult(obj, varargin)
            base=obj.v; dim=obj.dim;
            n=nargin-1;
            switch n
                case 1
                    a=varargin{1};
                    if isa(a,'MVector')
                        nd = max(dim, a.dim);
                        if nd==3, nv=[base(1)-a.v(1) base(2)-a.v(2) base(3)-a.v(3)];
                        else,      nv=[base(1)-a.v(1) base(2)-a.v(2) 0]; end
                    else
                        a=obj.ensureRowVec(a,'sub(array) expects 1x2 or 1x3 numeric.');
                        if numel(a)==2, nd=2; nv=[base(1)-obj.assertFiniteScalar(a(1),'x finite') base(2)-obj.assertFiniteScalar(a(2),'y finite') 0];
                        elseif numel(a)==3, nd=3; nv=[base(1)-obj.assertFiniteScalar(a(1),'x finite') base(2)-obj.assertFiniteScalar(a(2),'y finite') base(3)-obj.assertFiniteScalar(a(3),'z finite')];
                        else, error('MVector:sub:ArrayLen','sub(array) length must be 2 or 3.');
                        end
                    end
                case 2
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    nd=2; nv=[base(1)-x base(2)-y 0];
                case 3
                    x=obj.assertFiniteScalar(varargin{1},'x must be finite.');
                    y=obj.assertFiniteScalar(varargin{2},'y must be finite.');
                    z=obj.assertFiniteScalar(varargin{3},'z must be finite.');
                    nd=3; nv=[base(1)-x base(2)-y base(3)-z];
                otherwise
                    error('MVector:sub:Arity','sub accepts 1,2,3 args.');
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
        function [row, dim] = asRowWithDim(x)
            if isa(x,'MVector')
                row = x.v;
                dim = x.dim;
            elseif isnumeric(x) && isvector(x)
                if iscolumn(x), x=x.'; end
                n=numel(x);
                if n==2
                    row=[double(x(1)) double(x(2)) 0]; dim=2;
                elseif n==3
                    row=[double(x(1)) double(x(2)) double(x(3))]; dim=3;
                else
                    error('MVector:asRowWithDim:Len','Vector must have length 2 or 3.');
                end
            else
                error('MVector:asRowWithDim:Type','Operand must be MVector or numeric vector.');
            end
        end
        function [a,b]=assertTwoVectors(v1,v2,fname)
            if ~isa(v1,'MVector') || ~isa(v2,'MVector')
                error(sprintf('MVector:%s:Type',fname),'Both inputs must be MVector.');
            end
            a=v1; b=v2;
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
    end
end
