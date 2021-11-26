# Eigen-quadratic-programming

quadratic-programming with Eigen

### Standard QP problem

​				$\min \frac{1}{2}x^TPx+q^Tx $

subject to         $l \leq A_cX \leq u$

### MPC problem

​		$\arg\min_{x_k,u_k}\sum_{k=0}^N(x_k-x_{ref_{k}})^TQ(x_k-x_{ref_{k}})+u_k^TRu_k$

subject to  		$x_{k+1}=Ax_k+Bu_k$

​						    $x_{min}\leq x_k \leq u_{max}$

​							$u_{min}\leq u_k \leq u_{max}$

​							$x_0=x_0$

### MPC to QP

Hessian matrix P is equal to

$P = diag(Q, Q,\ ...,Q,R,\ ...,R)$

Gradient vector q is equal to

$q=[-Qx_{ref_{0}}, -Qx_{ref_{1}}, \ ..., -Qx_{ref_{N}}, 0, \ ...,0]^T$

Linear constraint matrix A is equal to

![linear_matrix](./images/linear_matrix.png)

Upper and lower bound are equal to

![constaint](/home/czk119/Desktop/ZJUDancer/Eigen-quadratic-programming/images/constaint.png)

### Test Example

$ \min \sum_{k=0}^N \frac{1}{2}(X-X_{Ref})Q(X-X_{Ref})^T+\frac{1}{2}URU^T  \\st. X_{Ref}=[0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1]  \\ X_{k+1} = \left[ \begin{matrix} 1 & \Delta t & \Delta t^2 \\ 0 & 1 & \Delta t  \\ 0 & 0 & 1\end{matrix} \right] X_k + \left[ \begin{matrix} 0 \\ 0 \\ 1 \end{matrix} \right]U_k \\ Q=diag(1, 1, 1) \\ R = 10\\X_0=0 \\ -1.0 \leq U \leq 1.0 \\ N=8, \Delta t=0.1$

### OSQP

