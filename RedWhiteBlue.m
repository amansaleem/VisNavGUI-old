function cm = RedWhiteBlue(n, gamma)

if nargin<1; n = 100; end
if nargin<2; gamma = 1; end

% for cyan
cm = ([n*ones(1,n), n:-1:0 ; ...
    0:n, ones(1,n)*n; ...
    0:n, ones(1,n)*n]' / n).^gamma;

% for Blue
cm = ([n*ones(1,n), n:-1:0 ; ...
    0:n, n-1:-1:0; ...
    0:n, ones(1,n)*n]' / n).^gamma;

% cm = ([n*ones(1,n), n:-1:0 ; ...
%     0:n, ones(1,n)*n; ...
%     n*ones(1,n+1), n-1:-1:0]' / n).^gamma;

cm(:,1) = cm(end:-1:1,1);
cm(:,2) = cm(end:-1:1,2);
cm(:,3) = cm(end:-1:1,3);
colormap(gcf,cm);
