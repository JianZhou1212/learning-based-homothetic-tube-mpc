function [V,theta,imax] = generate_polytope(n,m)

%V = generate_polytope(n,m) 
% Generate a hyperplane description {x : V*x <= 1} of a polytopic 
% approximation of the unit hypersphere in R^n, where the minimum
% angle between any two face normals is theta = pi/2^(m+1). 
%
% theta and the number of faces N are determined from (n,m) by:
% m\n | 2   . 3    . 4     . 5    . 6    | theta
% -----------------------------------------------------
%  0  | 4   . 6    . 8     . 10   . 12   | pi/2   90
%  1  | 8   . 18   . 32    . 50   . 72   | pi/4   45
%  2  | 16  . 42   . 80    . 130  . 192  | pi/8   22.5
%  3  | 32  . 90   . 176   . 290  . 432  | pi/16  11.3
%  4  | 64  . 295  . 781   . 1657 . 2777 | pi/32  5.6
%  5  | 128 . 2199 . 12734               | pi/64  2.8
%  6  | 256 .                  etc       | pi/128 1.4
%

tol = 1e-4;
if m > 0
  V = [eye(n);-eye(n)];
else
  V = zeros(0,n);
end
if n > 1
  for k = 1:m % theta = pi/2^k
    r = size(V,1);
    Vnew = [V;zeros((r-1)*r,n)];
    row_count = r;
    disc_count = 0;
    norm_count = 0;
    for i = 1:r-1
      U = V + circshift(V,i,1);
      for j = 1:r
        u = norm(U(j,:));
        if u > tol % remove null rows
          U(j,:) = U(j,:)/u;
          d = max(sum((ones(row_count,1)*U(j,:)).*Vnew(1:row_count,:),2));
          if abs(d - cos(pi/2^(k+1))) < tol % discard new rows unless the angle subtended with the nearest row of V is equal to pi/2^(k+1)
            row_count = row_count + 1;
            Vnew(row_count,:) = U(j,:);
          else
            disc_count = disc_count + 1;
          end
        else
          norm_count = norm_count + 1;
        end
      end
    end
    % fprintf(1,'k:%d r:%d r*(r-1):%d rows added:%d discarded:%d zero:%d\n',k,r,r*(r-1),row_count-r,disc_count,norm_count);
    V = Vnew(1:row_count,:);
  end
end
% fprintf(1,'total rows:%d\n',size(V,1));

if nargout > 1
  theta = zeros(size(V,1),1); 
  imax = zeros(size(V,1),1);
  for i = 1:size(V,1) 
    [theta(i),imax(i)] = max(sum((ones(size(V,1)-1,1)*V(i,:)).*[V(1:i-1,:);V(i+1:end,:)],2));
    theta(i) = acos(theta(i))*180/pi; % minimum angle in degrees between ith normal and any other normal
    if imax(i) >= i
      imax(i) = imax(i) + 1;
    end
  end
end
