classdef dobotLink
    %dobotLink 机器人连杆类
    %   DH参数，连杆的参数
    
    % 成员变量 ==============================================
    properties
        % DH参数
        alpha
        a
        d
        theta
        % m 质量
        m
    end
    
    % 类方法 ================================================
    methods
        function obj = dobotLink(alpha_,a_, d_, theta_, m_)
            % 构造函数
            %   利用此函数初始化类对象
            obj.alpha = alpha_;
            obj.a = a_;
            obj.d = d_;
            obj.theta = theta_;
            obj.m = m_;
        end
        
        function T = transMatrix(self,value)
            % transMatrix 计算此连杆的齐次变换矩阵
            % @param value: 旋转角度
            sa = sin(self.alpha);
            ca = cos(self.alpha);
            st = sin(value); ct = cos(value);
            dd = self.d;
            % 计算齐次变换矩阵
            T = [ ct    -st    0    self.a
                 st*ca  ct*ca  -sa  -sa*dd
                 st*sa  ct*sa  ca   ca*dd
                  0      0     0      1  ];
        end % function transMatrix
    end
end

