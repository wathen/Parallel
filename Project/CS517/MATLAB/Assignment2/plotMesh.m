function [hout,hnodes]=plotMesh(h,varargin)
%%=========================================================================
% function [hout,hnodes]=plotMesh(h,varargin)
%
% DESCRIPTION
%
%
% Input:
%  h        - description
%  varargin - description
%
% Output:
%  hout     - description
%  hnodes   - description
%
% see also
%%=========================================================================



%{
     @author Rowan Cockett

     Last Modified: 13-Mar-2013 22:28:42
%}


O.color = 'b';
O.linewidth = 1;
O.linestyle = '-';
O.showNodes = false;
O.spacing = 1;
O = myOptions(O,varargin);

mkvc = @(x)x(:);

dim = numel(h);

[h,curvy] = meshType(h);

if ~curvy
    h = cellfun(@(hj)[0;cumsum(hj(:))],h,'uniformoutput',false);
    [h{:}] = ndgrid(h{:});
end

if ~ishold
    cla;
end
washold = ishold;
hold on

switch dim
    case 2
        vH = cellfun(mkvc,h,'uniformoutput',false);
        [vX,vY] = deal(vH{:});
        m = size(h{1});

        %  node(i,j)          node(i,j+1)
        %       A -------------- B
        %       |                |
        %       |    cell(i,j)   |
        %       |        I       |
        %       |                |
        %      D -------------- C
        %  node(i+1,j)        node(i+1,j+1)

        % compute the normals
        [A,B] = indexCube('AB',m,[m(1),m(2)-1],[O.spacing,1]);
        vertsX = mkvc([vX(A)';vX(B)';nan(1,numel(A))]);
        vertsY = mkvc([vY(A)';vY(B)';nan(1,numel(A))]);
        [A,D] = indexCube('AD',m,[m(1)-1,m(2)],[1,O.spacing]);
        vertsX = [vertsX;mkvc([vX(A)';vX(D)';nan(1,numel(A))])];
        vertsY = [vertsY;mkvc([vY(A)';vY(D)';nan(1,numel(A))])];

        %plot the boundary
        vertsX = [vertsX;mkvc(h{1}(:,end));flipud(mkvc(h{1}(end,:)))];
        vertsY = [vertsY;mkvc(h{2}(:,end));flipud(mkvc(h{2}(end,:)))];

        h = plot(vertsY,vertsX,'color',O.color);
        if O.showNodes
            uverts = unique([vertsX,vertsY],'rows');
            hnodes = plot(uverts(:,2),uverts(:,1),'r.');
        else
            hnodes = [];
        end;

        set(gca,'ydir','reverse')
        xlabel('Y')
        ylabel('X')
        view(2);
        rotate3d off;

    case 3
        vH = cellfun(mkvc,h,'uniformoutput',false);
        [vX,vY,vZ] = deal(vH{:});
        m = size(h{1});


        %        node(i,j,k+1)       node(i,j+1,k+1)
        %            E --------------- F
        %           /|               / |
        %          / |              /  |
        %         /  |             /   |
        %  node(i,j,k)         node(i,j+1,k)
        %       A -------------- B     |
        %       |    H ----------|---- G
        %       |   /cell(i,j)   |   /
        %       |  /     I       |  /
        %       | /              | /
        %       D -------------- C
        %  node(i+1,j,k)      node(i+1,j+1,k)

        vertsX = [];
        vertsY = [];
        vertsZ = [];
        [A,D] = indexCube('AD',m,[m(1)-1,m(2),m(3)],[1,O.spacing,O.spacing]);
        vertsX = [vertsX;mkvc([vX(A)';vX(D)';nan(1,numel(A))])];
        vertsY = [vertsY;mkvc([vY(A)';vY(D)';nan(1,numel(A))])];
        vertsZ = [vertsZ;mkvc([vZ(A)';vZ(D)';nan(1,numel(A))])];
        [A,B] = indexCube('AB',m,[m(1),m(2)-1,m(3)],[O.spacing,1,O.spacing]);
        vertsX = [vertsX;mkvc([vX(A)';vX(B)';nan(1,numel(A))])];
        vertsY = [vertsY;mkvc([vY(A)';vY(B)';nan(1,numel(A))])];
        vertsZ = [vertsZ;mkvc([vZ(A)';vZ(B)';nan(1,numel(A))])];
        [A,E] = indexCube('AE',m,[m(1),m(2),m(3)-1],[O.spacing,O.spacing,1]);
        vertsX = [vertsX;mkvc([vX(A)';vX(E)';nan(1,numel(A))])];
        vertsY = [vertsY;mkvc([vY(A)';vY(E)';nan(1,numel(A))])];
        vertsZ = [vertsZ;mkvc([vZ(A)';vZ(E)';nan(1,numel(A))])];

        uverts = unique([vertsX,vertsY,vertsZ],'rows');



        h = plot3(vertsY,vertsZ,vertsX);
        if O.showNodes
            hnodes = plot3(uverts(:,2),uverts(:,3),uverts(:,1),'r.');
        else
            hnodes = [];
        end;

        set(gca,'zdir','reverse')
        xlabel('Y')
        ylabel('Z')
        zlabel('X')
        view(3);
        rotate3d on;

    otherwise
        error('grid must be 2D or 3D')
end


axis equal
if ~washold, hold off; end;

if nargout > 0;hout = h;end

end