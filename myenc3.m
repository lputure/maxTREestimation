function [x,z] = myenc3(B,w)
% weighted minimum bound circle for exactly 3 points
%
% 要求点数大于2或者3

% 只有两个点，返回加权值
if size(B,1) == 2    
    x1 = (w(1)*B(1,:)+w(2)*B(2,:))/(w(1)+w(2));
    z1 = w(1)*norm(x1-B(1,:));    
    x = x1; z = z1;
    return    
end

if w(1)==w(2)&&w(2)==w(3) %权重全部相等的情况,调用正常的函数处理
    [c,r] = enc3(B(:,1),B(:,2));
    z = w(1)*r;
    x = c;
    return;
end

% 权重只有一组相等的情况
if w(1)==w(2)&&w(2)~=w(3)
    wn = w.*w;
    
    x1 = (B(1,:)+B(2,:))/2; %特殊处理的部分
    z1 = w(1)*norm(x1-B(1,:));
    if z1>w(3)*norm(x1-B(3,:))
        x = x1; z = z1;
        return;
    end
    
    c2 = (wn(1)*B(1,:)-wn(3)*B(3,:))/(wn(1)-wn(3));
    x2 = (w(1)*B(1,:)+w(3)*B(3,:))/(w(1)+w(3));
    z2 = w(1)*norm(x2-B(1,:));
    if z2>w(2)*norm(x2-B(2,:))
        x = x2; z = z2;
        return
    end
    
    c3 = (wn(3)*B(3,:)-wn(2)*B(2,:))/(wn(3)-wn(2));
    x3 = (w(3)*B(3,:)+w(2)*B(2,:))/(w(3)+w(2));
    z3 = w(3)*norm(x3-B(3,:));
    if z3>w(1)*norm(x3-B(1,:))
        x = x3; z = z3;
        return
    end
    
    r2 = norm(x2-c2);
    r3 = norm(x3-c3);
    
    [ix,iy] = circcirc(c2(1),c2(2),r2,c3(1),c3(2),r3);
    A = [B(2,:)'-B(1,:)',B(3,:)'-B(1,:)'];
    b = [[ix(1);iy(1)]-B(1,:)', [ix(2);iy(2)]-B(1,:)'];
    y = A\b;
    
    if(all(y(:,1)>=0 & y(:,1)<=1))
        x = [ix(1),iy(1)];
        z = w(1)*norm(x-B(1,:));
        return
    end
    if(all(y(:,2)>=0 & y(:,2)<=1))
        x = [ix(2),iy(2)];
        z = w(1)*norm(x-B(1,:));
        return
    end
end

if w(1)==w(3)&&w(3)~=w(2)
    wn = w.*w;
    
    c1 = (wn(1)*B(1,:)-wn(2)*B(2,:))/(wn(1)-wn(2));
    x1 = (w(1)*B(1,:)+w(2)*B(2,:))/(w(1)+w(2));
    z1 = w(1)*norm(x1-B(1,:));
    f1 = z1>w(3)*norm(x1-B(3,:));
    if f1
        x = x1; z = z1;
        return
    end
    
    x2 = (B(1,:)+B(3,:))/2; %特殊处理的部分
    z2 = w(1)*norm(x1-B(1,:));
    if z2>w(3)*norm(x2-B(2,:))
        x = x2; z = z2;
        return;
    end
    
    c3 = (wn(3)*B(3,:)-wn(2)*B(2,:))/(wn(3)-wn(2));
    x3 = (w(3)*B(3,:)+w(2)*B(2,:))/(w(3)+w(2));
    z3 = w(3)*norm(x3-B(3,:));
    if z3>w(1)*norm(x3-B(1,:))
        x = x3; z = z3;
        return
    end
    
    r1 = norm(x1-c1);
    r3 = norm(x3-c3);
    
    [ix,iy] = circcirc(c1(1),c1(2),r1,c3(1),c3(2),r3);
    A = [B(2,:)'-B(1,:)',B(3,:)'-B(1,:)'];
    b = [[ix(1);iy(1)]-B(1,:)', [ix(2);iy(2)]-B(1,:)'];
    y = A\b;
    
    if(all(y(:,1)>=0 & y(:,1)<=1))
        x = [ix(1),iy(1)];
        z = w(1)*norm(x-B(1,:));
        return
    end
    if(all(y(:,2)>=0 & y(:,2)<=1))
        x = [ix(2),iy(2)];
        z = w(1)*norm(x-B(1,:));
        return
    end
end

if w(2)==w(3)&&w(1)~=w(2)
    wn = w.*w;
    
    c1 = (wn(1)*B(1,:)-wn(2)*B(2,:))/(wn(1)-wn(2));
    x1 = (w(1)*B(1,:)+w(2)*B(2,:))/(w(1)+w(2));
    z1 = w(1)*norm(x1-B(1,:));
    f1 = z1>w(3)*norm(x1-B(3,:));
    if f1
        x = x1; z = z1;
        return
    end
    
    c2 = (wn(1)*B(1,:)-wn(3)*B(3,:))/(wn(1)-wn(3));
    x2 = (w(1)*B(1,:)+w(3)*B(3,:))/(w(1)+w(3));
    z2 = w(1)*norm(x2-B(1,:));
    f2 = z2>w(2)*norm(x2-B(2,:));
    if f2
        x = x2; z = z2;
        return
    end
    
    x3 = (B(2,:)+B(3,:))/2; %特殊处理的部分
    z3 = w(2)*norm(x3-B(2,:));
    if z3>w(1)*norm(x3-B(1,:))
        x = x3; z = z3;
        return;
    end
    
    r1 = norm(x1-c1);
    r2 = norm(x2-c2);
    
    [ix,iy] = circcirc(c1(1),c1(2),r1,c2(1),c2(2),r2);
    
    A = [B(2,:)'-B(1,:)',B(3,:)'-B(1,:)'];
    b = [[ix(1);iy(1)]-B(1,:)', [ix(2);iy(2)]-B(1,:)'];
    y = A\b;
    
    if(all(y(:,1)>=0 & y(:,1)<=1))
        x = [ix(1),iy(1)];
        z = w(1)*norm(x-B(1,:));
        return
    end
    if(all(y(:,2)>=0 & y(:,2)<=1))
        x = [ix(2),iy(2)];
        z = w(1)*norm(x-B(1,:));
        return
    end
    
end
% 权重都不相等的情况

wn = w.*w;
c1 = (wn(1)*B(1,:)-wn(2)*B(2,:))/(wn(1)-wn(2));
x1 = (w(1)*B(1,:)+w(2)*B(2,:))/(w(1)+w(2));
z1 = w(1)*norm(x1-B(1,:));
f1 = z1>w(3)*norm(x1-B(3,:));
if f1
    x = x1; z = z1;
    return
end

c2 = (wn(1)*B(1,:)-wn(3)*B(3,:))/(wn(1)-wn(3));
x2 = (w(1)*B(1,:)+w(3)*B(3,:))/(w(1)+w(3));
z2 = w(1)*norm(x2-B(1,:));
f2 = z2>w(2)*norm(x2-B(2,:));
if f2
    x = x2; z = z2;
    return
end

% c3 = (wn(3)*B(3,:)-wn(2)*B(2,:))/(wn(3)-wn(2));
x3 = (w(3)*B(3,:)+w(2)*B(2,:))/(w(3)+w(2));
z3 = w(3)*norm(x3-B(3,:));
f3 = z3>w(1)*norm(x3-B(1,:));
if f3
    x = x3; z = z3;
    return
end

r1 = norm(x1-c1);
r2 = norm(x2-c2);

[ix,iy] = circcirc(c1(1),c1(2),r1,c2(1),c2(2),r2);
A = [B(2,:)'-B(1,:)',B(3,:)'-B(1,:)'];
b = [[ix(1);iy(1)]-B(1,:)', [ix(2);iy(2)]-B(1,:)'];
y = A\b;

if(all(y(:,1)>=0 & y(:,1)<=1))
    x = [ix(1),iy(1)];
    z = w(1)*norm(x-B(1,:));
    return
end
if(all(y(:,2)>=0 & y(:,2)<=1))
    x = [ix(2),iy(2)];
    z = w(1)*norm(x-B(1,:));
    return
end

if 1
    fprintf('error\n');
end

end

function [center,radius] = enc3(X,Y)
% minimum radius enclosing circle for exactly 3 points
%
% x, y are 3x1 vectors

% convert to complex
xy = X + sqrt(-1)*Y;

% just in case the points are collinear or nearly so, get
% the interpoint distances, and test the farthest pair
% to see if they work.
Dij = @(XY,i,j) abs(XY(i) - XY(j));
D12 = Dij(xy,1,2);
D13 = Dij(xy,1,3);
D23 = Dij(xy,2,3);

% Find the most distant pair. Test if their circumcircle
% also encloses the third point.
if (D12>=D13) && (D12>=D23)
    center = (xy(1) + xy(2))/2;
    radius = D12/2;
    if abs(center - xy(3)) <= radius
        center = [real(center),imag(center)];
        return
    end
elseif (D13>=D12) && (D13>=D23)
    center = (xy(1) + xy(3))/2;
    radius = D13/2;
    if abs(center - xy(2)) <= radius
        center = [real(center),imag(center)];
        return
    end
elseif (D23>=D12) && (D23>=D13)
    center = (xy(2) + xy(3))/2;
    radius = D23/2;
    if abs(center - xy(1)) <= radius
        center = [real(center),imag(center)];
        return
    end
end

% if we drop down to here, then the points cannot
% be collinear, so the resulting 2x2 linear system
% of equations will not be singular.
A = 2*[X(2)-X(1), Y(2)-Y(1); X(3)-X(1), Y(3)-Y(1)];
rhs = [X(2)^2 - X(1)^2 + Y(2)^2 - Y(1)^2; ...
    X(3)^2 - X(1)^2 + Y(3)^2 - Y(1)^2];

center = (A\rhs)';
radius = norm(center - [X(1),Y(1)]);
end