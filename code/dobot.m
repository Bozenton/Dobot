classdef dobot
    %UNTITLED2 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        link1
        link2
        link3
        linkT
        len
        mass
        center
        inertia
    end
    
    methods
        function obj = dobot(initStruct)
            % dobot 构造此类的实例
            %   此处显示详细说明
            if ~isfield(initStruct, 'mass') || ...
                    ~isfield(initStruct, 'length')
                error("The initStruct does not contain field of mass or length");
            end
            if ~isfield(initStruct, 'center') || ...
                    ~isfield(initStruct, 'inertia')
                error("The initStruct does not contain field of center or inertia");
            end
            obj.len = initStruct.length;
            obj.mass = initStruct.mass;
            l1 = obj.len(1); l2 = obj.len(2); l3 = obj.len(3);
            m1 = obj.mass(1); m2 = obj.mass(2); m3 = obj.mass(3);
            
            obj.link1 = dobotLink(0, 0, l1, 0, m1);
            obj.link2 = dobotLink(pi/2, 0, 0, 0, m2);
            obj.link3 = dobotLink(0, l2, 0, 0, m3);
            obj.linkT = dobotLink(0, l3, 0, 0, 0);
            
            obj.center = initStruct.center;
            obj.inertia = initStruct.inertia;
        end
        
        % ------------------------------------------------------------
        % 正运动学
        % ------------------------------------------------------------
        function fT = fkine(self,theta)
            % fkine 正运动学
            %   此处显示详细说明
            if length(theta)~=3
                error("The input theta should have 3 values");
            end
            fT = self.link1.transMatrix(theta(1)) * ...
            self.link2.transMatrix(theta(2)) * ...
            self.link3.transMatrix(theta(3)) * ...
            self.linkT.transMatrix(0);
        end
        
        % ------------------------------------------------------------
        % 逆运动学
        % ------------------------------------------------------------
        function tt = invkine(self, xyz)
            x = xyz(1);
            y = xyz(2);
            z = xyz(3);
            l1 = self.len(1); l2 = self.len(2); l3 = self.len(3);
            
            ll = sqrt(x^2+y^2+z^2);
            if ll > l2+l3
                error("It has exceed the working space");
            end
            
            % 根据几何关系，可以首先算出theta3
            c3 = (x*x+y*y+z*z-l2*l2-l3*l3)/(2*l2*l3);
            tt3 = -acos(c3);
            % theta1也很容易确定
            tt1 = atan2(y, x);
            % theta2
            beta = atan2(z, sqrt(x*x+y*y) );
            phi = acos( (x^2+y^2+z^2+l2^2-l3^2) / (2*l2*sqrt(x^2+y^2+z^2)) );
            tt2 = beta + phi;
            tt = [tt1; tt2; tt3];
        end 
        
        
        % ------------------------------------------------------------
        % 计算雅可比矩阵
        % ------------------------------------------------------------
        function J = jacobi(self, theta)
            q1 = theta(1); q2 = theta(2); q3 = theta(3);
            l1 = self.len(1); l2 = self.len(2); l3 = self.len(3);
            J1x = -sin(q1)*(l2*cos(q2)+l3*cos(q2+q3));
            J1y = cos(q1)*(l2*cos(q2)+l3*cos(q2+q3));
            J1z = 0;
            J2x = cos(q1)*(-l2*sin(q2)-l3*sin(q2+q3));
            J2y = sin(q1)*(-l2*sin(q2)-l3*sin(q2+q3));
            J2z = l2*cos(q2) + l3*cos(q2+q3);
            J3x = cos(q1)*(-l3*sin(q2+q3));
            J3y = sin(q1)*(-l3*sin(q2+q3));
            J3z = l3*cos(q2+q3);
            J1w = [0;0;1];
            J2w = [sin(q1); -cos(q1); 0];
            J3w = [sin(q1); -cos(q1); 0];
            J = [
                J1x J2x J3x
                J1y J2y J3y
                J1z J2z J3z
            ];
            J = [J; [J1w, J2w, J3w]];
        end
        
        function [H1, H2, H3] = Hessian(self, theta)
        l1 = self.len(1); l2 = self.len(2); l3 = self.len(3);
        th1 = theta(1);
        th2 = theta(2);
        th3 = theta(3);

        H1 = reshape([-cos(th1).*(l3.*cos(th2+th3)+l2.*cos(th2)),sin(th1).*(l3.*sin(th2+th3)+l2.*sin(th2)),l3.*sin(th2+th3).*sin(th1),sin(th1).*(l3.*sin(th2+th3)+l2.*sin(th2)),-cos(th1).*(l3.*cos(th2+th3)+l2.*cos(th2)),-l3.*cos(th2+th3).*cos(th1),l3.*sin(th2+th3).*sin(th1),-l3.*cos(th2+th3).*cos(th1),-l3.*cos(th2+th3).*cos(th1)],[3,3]);
        H2 = reshape([-sin(th1).*(l3.*cos(th2+th3)+l2.*cos(th2)),-cos(th1).*(l3.*sin(th2+th3)+l2.*sin(th2)),-l3.*sin(th2+th3).*cos(th1),-cos(th1).*(l3.*sin(th2+th3)+l2.*sin(th2)),-sin(th1).*(l3.*cos(th2+th3)+l2.*cos(th2)),-l3.*cos(th2+th3).*sin(th1),-l3.*sin(th2+th3).*cos(th1),-l3.*cos(th2+th3).*sin(th1),-l3.*cos(th2+th3).*sin(th1)],[3,3]);
        H3 = reshape([0.0,0.0,0.0,0.0,-l3.*sin(th2+th3)-l2.*sin(th2),-l3.*sin(th2+th3),0.0,-l3.*sin(th2+th3),-l3.*sin(th2+th3)],[3,3]);

        end
        
        function fplot(self, theta)
            % fplot 绘制机械臂图像
            %   此处显示详细说明
            if length(theta)~=3
                error("The input theta should have 3 values");
            end
            T1 = self.link1.transMatrix(theta(1));
            T2 = T1 * self.link2.transMatrix(theta(2));
            T3 = T2 * self.link3.transMatrix(theta(3));
            TT = T3 * self.linkT.transMatrix(0);
            Tmatrix(:, :, 1) = T1;
            Tmatrix(:, :, 2) = T2;
            Tmatrix(:, :, 3) = T3;
            Tmatrix(:, :, 4) = TT;
            num = size(Tmatrix, 3);     % 坐标点的个数（不含原点）
            figureName = string(datetime('now'));
            figure('Name',figureName);
            % 绘制旋转的坐标系
            for i = 1:num
                org = Tmatrix(1:3, 4, i);
                xyz = Tmatrix(1:3, 1:3, i)/3.0;
                clr = rand(1,3);
                txt = ["x", "y", "z"];
                for j = 1:3
                    quiver3(org(1), org(2), org(3), ...
                        xyz(1,j),xyz(2,j),xyz(3,j), ...
                        1, 'filled','Color',clr);
                    text(org(1)+xyz(1,j), org(2)+xyz(2,j), org(3)+xyz(3,j), ...
                        txt(j));
                    hold on
                end
            end
            % 绘制连杆
            orglast = [0 0 0];
            for i=1:num
                org = Tmatrix(1:3, 4, i);
                quiver3(orglast(1), orglast(2), orglast(3), ...
                        org(1)-orglast(1), org(2)-orglast(2), org(3)-orglast(3), ...
                        1, 'filled','Color',[0 0 0]);
                orglast = org;
                hold on
            end
            grid on;%绘网格
            xlabel('x');ylabel('y');zlabel('z');
            axis equal;
            title("The position and orientation of the links")
        end
        
        % ----------------------------------------------------
        % dynamic: 动力学
        %       通过牛顿欧拉法计算关节力矩
        % @param th: 关节角度
        % @param dth: 关节角速度
        % @param ddth: 关节角加速度
        % @return: 关节力矩（只需取z轴分量即可）
        % ----------------------------------------------------
        function n = dynamic(self, th, dth, ddth)
            g = 9.8; % 重力加速度
            
            num = length(self.center); % 连杆的数量
            m = self.mass;      % 质量
            I = self.inertia;   % 转动惯量
            
            % 使用牛顿欧拉法时需要用到的参数
            omega = zeros(3,num);
            domega = zeros(3,num);
            v = zeros(3, num);
            dv = zeros(3,num);
            vc = zeros(3, num);
            dvc = zeros(3,num);
            F = zeros(3, num);
            N = zeros(3, num);
            f = zeros(3, 4);
            n = zeros(3, 4);
            
            % 旋转矩阵
            R = zeros([3,3,num+1]);
            Rinv = zeros([3,3,num+1]);
            temp1 = self.link1.transMatrix(th(1));
            temp2 = self.link2.transMatrix(th(2));
            temp3 = self.link3.transMatrix(th(3));
            tempT = self.linkT.transMatrix(0);
            R(:, :, 1) = temp1(1:3, 1:3);
            R(:, :, 2) = temp2(1:3, 1:3);
            R(:, :, 3) = temp3(1:3, 1:3);
            R(:, :, 4) = eye(3);
            for i = 1:size(Rinv, 3)
                Rinv(:,:,i) = R(:,:,i)';
            end
            
            % 位置信息
            % 第i个坐标系中，第i+1个坐标系的原点的位置
            P = [temp2(1:3,4), temp3(1:3,4), tempT(1:3,4)];
            % 第i个坐标系中，第i个连杆的位置
            Pc = self.center;

            % 外推
            for i = 0:num-1
                if(i == 0)
                    omega(:, i+1) = dth(i+1)*[0;0;1];
                    domega(:,i+1) = ddth(i+1)*[0;0;1];
                    dv(:, i+1) = Rinv(:,:,i+1) * [0;0;-g];
                else
                    omega(:, i+1) = Rinv(:,:,i+1)*omega(:,i) + dth(i+1)*[0;0;1];
                    domega(:,i+1) = Rinv(:,:,i+1)*domega(:,i) + cross(Rinv(:,:,i+1)*omega(:,i), dth(i+1)*[0;0;1]) + ddth(i+1)*[0;0;1];
                    dv(:, i+1) = Rinv(:,:,i+1)*( cross(domega(:,i),P(:,i))+cross(omega(:,i),cross(omega(:,i),P(:,i))) + dv(:,i) );
                end

                dvc(:, i+1) = cross(domega(:,i+1), Pc(:,i+1)) + cross(omega(:,i+1), cross(omega(:,i+1),Pc(:,i+1))) + dv(:,i+1);
                F(:, i+1) = m(i+1)*dvc(:, i+1);
                N(:, i+1) = I(:,:,i+1)*domega(:,i+1) + cross(omega(:,i+1), I(:,:,i+1)*omega(:,i+1));
            end

            % 内推
            for i = num:-1:1
                f(:,i) = R(:,:,i+1)*f(:,i+1) + F(:,i);
                n(:,i) = N(:,i) + R(:,:,i+1)*n(:,i+1) + cross(Pc(:,i), F(:,i))+ cross(P(:,i), R(:,:,i+1)*f(:,i+1));
            end
        end % function dynamic end
    end
end

