function x_new = newton_raphson(f_handle, fdx_handle, x0, max_iter, tol)

x_old = x0;
x_new = x0;
diff_rel = 1;

f = f_handle(x0);
fdx = fdx_handle(x0); 

iter = 1;
while (abs(f) > tol && diff_rel > tol && iter <= max_iter)
    
   x_new = x_old - f/fdx; 
   
   iter = iter + 1;    
   diff_rel = abs((x_new - x_old)/x_old);  
   
   x_old = x_new;  
   f = f_handle(x_new); 
   fdx = fdx_handle(x_new); 
end