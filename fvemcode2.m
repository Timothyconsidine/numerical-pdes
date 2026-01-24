function [AA,BB,QQ,QQprime] = Get(obj,alpha) %Get function
            
            %Constructing A if kappa depends on u
            if obj.bvp.second_order_component_is_linear == false
                obj.A = AssembleA(obj,alpha);
            end
            
            %Constructing Q if f exists
            if obj.bvp.f_exists == true
                obj.Q = AssembleQ(obj,alpha);
            end
            
            %Constructing Qprime if f is linear
            if obj.bvp.f_is_linear == false
                obj.Qprime = AssembleQprime(obj,alpha);
            end

            AA = obj.A;
            BB = obj.B;
            QQ = obj.Q;
            QQprime = obj.Qprime;

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
            kappa = obj.bvp.kappa;

            %Main diagonal
            A1 = [(kappa(z(2), alphabar(2)))/h(1);
                  (kappa(z(2:a-1), alphabar(2:a-1))./h(1:a-2))+...
                  ((kappa(z(3:a), alphabar(3:a))./h(2:a-1)));
                  (kappa(z(a), alphabar(a))/h(a-1))];

            %Lower diagonal
            A2 = -(kappa(z(2:a), alphabar(2:a))./h(1:a-1));
            
            %Upper diagonal
            A3 = -((kappa(z(2:a), alphabar(2:a))./h(1:a-1)));