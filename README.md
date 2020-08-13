# Mathieu Functions of the first kind

Mathieu functions of first kind are those which are simply periodic, with period either ![\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpi) or ![2\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+2%5Cpi), that satisfy the equation

![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{equation}
\dfrac{d^2f}{dy^2} + \left[a-2q\cos{(2y)} \right]f = 0 \nonumber
\end{equation}](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Bequation%7D%0A%5Cdfrac%7Bd%5E2f%7D%7Bdy%5E2%7D+%2B+%5Cleft%5Ba-2q%5Ccos%7B%282y%29%7D+%5Cright%5Df+%3D+0+%5Cnonumber%0A%5Cend%7Bequation%7D)


characterized by the coefficient being periodic over the domain. There are two
free parameters, namely q and a, and in the case of simply periodic the *characteristic number* ![a = a(q)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a+%3D+a%28q%29), resulting in a necessary condition for the existance of simply periodic solutions. In other words, ![f(q, y)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+f%28q%2C+y%29) is a solution of (1) whenever ![a = a(q)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a+%3D+a%28q%29) has finite, unique values.


According to Floquet theory, Mathieu's equation can only one periodic solution, which can be either even or odd, ![\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpi) or ![2\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+2%5Cpi) periodic. With this in mind, there are 4 classes of identifyble (simply periodic) solutions to Mathieu equation, in accordance to these properties. These are

1. Even, **cosine elliptic** function, with period ![\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpi).

![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\; ce_{2n}(z,q) = \sum_{r=0}^{\infty}A_{2r}^{(2n)}(q)\cos{(2rz)}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+ce_%7B2n%7D%28z%2Cq%29+%3D+%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DA_%7B2r%7D%5E%7B%282n%29%7D%28q%29%5Ccos%7B%282rz%29%7D)


2. Odd, **Cosine elliptic** function, with period ![2\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+2%5Cpi).


![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;ce_{2n+1}(z,q) =  \sum_{r=0}^{\infty}A_{2r+1}^{(2n+1)}(q) \cos{((2r+1)z)}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3Bce_%7B2n%2B1%7D%28z%2Cq%29+%3D++%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DA_%7B2r%2B1%7D%5E%7B%282n%2B1%29%7D%28q%29+%5Ccos%7B%28%282r%2B1%29z%29%7D%0A)

3. Odd, **sine elliptic** function, with period ![2\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+2%5Cpi).

![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;se_{2n+1}(z,q) =  \sum_{r=0}^{\infty}B_{2r+1}^{(2n+1)}(q) \sin{((2r+1)z)}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3Bse_%7B2n%2B1%7D%28z%2Cq%29+%3D++%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DB_%7B2r%2B1%7D%5E%7B%282n%2B1%29%7D%28q%29+%5Csin%7B%28%282r%2B1%29z%29%7D%0A)


4. Even, **sine elliptic** function, with period ![\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpi).

![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;se_{2n+2}(z,q) =  \sum_{r=0}^{\infty}B_{2r+2}^{(2n+1)}(q) \sin{((2r+2)z)}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3Bse_%7B2n%2B2%7D%28z%2Cq%29+%3D++%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DB_%7B2r%2B2%7D%5E%7B%282n%2B1%29%7D%28q%29+%5Csin%7B%28%282r%2B2%29z%29%7D%0A)


**cosine** eliptic functions cannot coexist with **sine elliptic** functions, since only one periodic solution is allows, but a second, non-periodic solution can be obtained from each of them. 


The pair of **characteristic numbers** ![a_{n}(q), b_{n}(q)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a_%7Bn%7D%28q%29%2C+b_%7Bn%7D%28q%29) are each associated with the mode n, and correspond to **ce** and **se** functions respectively. Note that ![a_{n}(q) \neq b_{n}(q)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a_%7Bn%7D%28q%29+%5Cneq+b_%7Bn%7D%28q%29). 









