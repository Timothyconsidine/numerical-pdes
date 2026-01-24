% FVEMSystem Class, MATH 4340
% Timothy Considine
classdef FVEMSystem

    properties (SetAccess = public)
        partitions;
        bvp;
        A;
        B;
        Q;
        Qprime;
    end

    methods

        function obj = FVEMSystem(partition, twopointbvp) %Constructor

            obj.partitions = partition;
            obj.bvp = twopointbvp;
            alpha = sparse(obj.partitions.numintervals+1,1);
            
            %Constructing A if kappa does not depend on u
            if obj.bvp.second_order_component_is_linear == true
                obj.A = AssembleA(obj,alpha);
            end
            
            %Constructing Qprime if f exists and is linear
            if obj.bvp.f_exists == true && obj.bvp.f_is_linear == true
                obj.Qprime = AssembleQprime(obj,alpha);
            end

            %Intializing boundary vector B
            obj.B = sparse([obj.bvp.boundary_values(1); 
                         zeros(obj.partitions.numintervals-1,1); 
                         obj.bvp.boundary_values(2)]);
           
            %Checking if left boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(1) == true
                obj.B(1) = -obj.bvp.boundary_values(1);
            end
            
            %Checking if right boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(2) == true
                obj.B(obj.partitions.numintervals+1) = -obj.bvp.boundary_values(2);
            end
            
        end

        function [AA,BB,QQ,QQprime] = Get(obj,alpha) %Get function
            
            %Constructing A if kappa depends on u
            if obj.bvp.second_order_component_is_linear == false
                obj.A = AssembleA(obj,alpha);
            end
            
            %Constructing Q if f exists
            if obj.bvp.f_exists == true
                obj.Q = AssembleQ(obj,alpha);

            %Constructing Qprime if f is non-linear
                 if obj.bvp.f_is_linear == false
                    obj.Qprime = AssembleQprime(obj,alpha);
                 end

            end
            
            AA = obj.A;
            QQprime = obj.Qprime;
            QQ = obj.Q;
            BB = obj.B;

        end

    end

     methods (Access = private)
        
         function AA = AssembleA(obj,alpha) %Constructing A matrix
            
            %Intermediate variables
            alphabar = [0;(alpha(1:(length(obj.partitions.x)-1))+...
                        alpha(2:length(obj.partitions.x)))./2];
            z = obj.partitions.z;
            a = length(alpha);
            h = obj.partitions.h;
          
            kappa = [obj.bvp.kappa(z(1:end-1),alphabar(1:end))];

            %Main diagonal
            A1 = [(kappa(2))/h(1);
                  (kappa((2:a-1))./h(1:a-2))+...
                  ((kappa((3:a))./h(2:a-1)));
                  (kappa(a)/h(a-1))];
            
            %Lower diagonal
            A2 = -[kappa(2:a)./h(1:a-1) ; 0];

            %Upper diagonal
            A3 =  -[0 ; ((kappa(2:a)./h(1:a-1)))];

            %Creating tri-diagonal matrix
            AA = sparse(a,a);
            AA = spdiags([A2 A1 A3],-1:1,AA);
            
            %Checking if left boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(1) == true
                AA(1,1) = 1;
                AA(1,2) = 0;
            else
               AA(1,1) = (kappa(2)*alphabar(2))/(h(1))-obj.bvp.gamma_values(1);
               %AA(1,2) = (kappa(2)*alphabar(2))/(h(1));
            end
            
            %Checking if right boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(2) == true
                AA(a,a) = 1;
                AA(a,a-1) = 0;
            else 
                AA(a,a) = (kappa(a)*alphabar(a-1))/(h(a-1))-obj.bvp.gamma_values(2);
                %AA(a,a-1) = -(kappa(a)*alphabar(a-1))/(h(a-1));
            end
           
         end

         function QQ = AssembleQ(obj,alpha) %Constructing Q vector
            
            %Intermediate Variables
            x = obj.partitions.x;
            a = length(alpha);
            h = obj.partitions.h;
            f = [obj.bvp.f(x,alpha)];

            %Constructing Q 
            QQ = [f(1).*(h(1)/2);
                 f(2:a-1).*((h(1:a-2)+h(2:a-1))./2);
                 f(a).*(h(a-1)/2)];
            
            %Checking if left boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(1)
                QQ(1) = 0;
            end

            %Checking if right boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(2)
                QQ(a) = 0;
            end
           
         end

         function QQprime = AssembleQprime(obj,alpha) %Constructing Qprime Matrix
            
            %Intermediate Variables
            x = obj.partitions.x;
            a = length(alpha);
            h = obj.partitions.h;
            d2f = [obj.bvp.D2f(x,alpha)];

            %Constructing Qprime
            QQp = [d2f(1).*(h(1)/2);
                  d2f(2:a-1).*((h(1:a-2)+h(2:a-1))./2);
                  d2f(a).*(h(a-1)/2)];

            %Checking if left boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(1)
                QQp(1)=0;
            end
            
            %Checking if right boundary is dirichlet and setting accordingly
            if obj.bvp.type_is_Dirichlet(2)
                QQp(a)=0;
            end

            QQprime = sparse(a,a);
            QQprime = spdiags(QQp,0,QQprime);
            
         end

     end

     
end



