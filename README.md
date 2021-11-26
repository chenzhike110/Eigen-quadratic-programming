# Eigen-quadratic-programming

quadratic-programming with Eigen

### Standard QP problem

​				$\min \frac{1}{2}x^TPx+q^Tx $

subject to         $l \leq A_cX \leq u$

### MPC problem

​		$\arg\min_{x_k,u_k}\sum_{k=0}^N(x_k-x_r)^TQ(x_k-x_r)+u_k^TRu_k$

subject to  		$x_{k+1}=Ax_k+Bu_k$

​						    $x_{min}\leq x_k \leq u_{max}$

​							$u_{min}\leq u_k \leq u_{max}$

​							$x_0=x_0$

### OSQP

