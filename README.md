# Mathieu Functions of the first kind

Mathieu's equation is the simplest case of a Hill's equation, a general 2nd order ODE with periodic coefficients, and is given by

<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Bequation%7D%0A%5Cdfrac%7Bd%5E2f%7D%7Bdy%5E2%7D+%2B+%5Cleft%5Ba-2q%5Ccos%7B%282y%29%7D+%5Cright%5Df+%3D+0%2C+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+0%5Cleq+y%5Cleq+2%5Cpi+%5Cnonumber%0A%5Cend%7Bequation%7D%0A" 
alt="\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{equation}
\dfrac{d^2f}{dy^2} + \left[a-2q\cos{(2y)} \right]f = 0, \;\;\;\;\;\;\;\;\;\; 0\leq y\leq 2\pi \nonumber
\end{equation}
">


We focus on **basically periodic** solutions to Mathieu equation that are <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpi" 
alt="\pi"> or <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+2%5Cpi" 
alt="2\pi"> periodic, and that depend on all of a, q and y, with the extra condition that two parameters are related, *i.e.* ![a = a(q)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a+%3D+a%28q%29)[[1]](#1). These solutions are also referred to as **Mathieu's functions of the first kind**, or simply **Mathieu functions**. In the limit of <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+q%5Crightarrow0" 
alt="q\rightarrow0">, each Mathieu function converges to one of <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Ccos%7B%28ny%29%7D" 
alt="\cos{(ny)}"> or <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Csin%7B%28ny%29%7D" 
alt="\sin{(ny)}">, with the exception of <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+ce_%7B0%7D%28q%2C+y%29+%5Crightarrow+1%2F%5Csqrt%7B2%7D" 
alt="ce_{0}(q, y) \rightarrow 1/\sqrt{2}">. In addition, in the same limit  <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a%3Dn%5E2" 
alt="a=n^2">, n an integer, as Mathieu's equation becomes the string equation. Following now standard notation, Mathieu's functions are written as <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+ce_%7Bn%7D%28q%2C+y%29" 
alt="ce_{n}(q, y)"> or <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+se_%7Bn%7D%28q%2C+y%29" 
alt="se_{n}(q, y)">, and usually referred to as **cosine elliptic** or **sine elliptic** functions respectively [[2]](#2).


There are 4 foundamental classes of simply-periodic solutions to Mathieu equation, each associated with a periodicity and symmetry relation, and consequenly each one satisfies a different combination of boundary conditions [[3]](#3). These are:

1. Even, **cosine elliptic** function, with period ![\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpi).

![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\; ce_{2n}(z,q) = \sum_{r=0}^{\infty}A_{2r}^{(2n)}(q)\cos{(2rz)}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+ce_%7B2n%7D%28z%2Cq%29+%3D+%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DA_%7B2r%7D%5E%7B%282n%29%7D%28q%29%5Ccos%7B%282rz%29%7D)

and

<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Bequation%7D%0A%5Cdfrac%7Bdce_%7B2n%7D%28q%2C+y%29%7D%7Bdy%7D%3D0%2C+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+%5Ctext%7Bat%7D%5C%3B%5C%3B+y%3D%5C%7B0%2C+%5Cpi%5C%7D%2C+%5Cnonumber%0A%5Cend%7Bequation%7D" 
alt="\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{equation}
\dfrac{dce_{2n}(q, y)}{dy}=0, \;\;\;\;\;\;\; \text{at}\;\; y=\{0, \pi\}, \nonumber
\end{equation}">


2. Odd, **Cosine elliptic** function, with period ![2\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+2%5Cpi).


![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;ce_{2n+1}(z,q) =  \sum_{r=0}^{\infty}A_{2r+1}^{(2n+1)}(q) \cos{((2r+1)z)}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3Bce_%7B2n%2B1%7D%28z%2Cq%29+%3D++%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DA_%7B2r%2B1%7D%5E%7B%282n%2B1%29%7D%28q%29+%5Ccos%7B%28%282r%2B1%29z%29%7D%0A)


and

<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Bequation%7D%0A%5Cdfrac%7Bdce_%7B2n%2B1%7D%7D%7Bdy%7D%3D0%2C+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+at+%5C%3B%5C%3B+y%3D0%2C++%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+ce_%7B2n%2B1%7D%3D0%2C%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+at%5C%3B%5C%3B%5C%3B+y%3D%5Cpi%5Cnonumber%0A%5Cend%7Bequation%7D" 
alt="\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{equation}
\dfrac{dce_{2n+1}}{dy}=0, \;\;\;\;\;\; at \;\; y=0,  \;\;\;\;\;\; \;\;\;\;\;\; ce_{2n+1}=0,\;\;\;\;\;\; at\;\;\; y=\pi\nonumber
\end{equation}">


3. Odd, **sine elliptic** function, with period ![2\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+2%5Cpi).

![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;se_{2n+1}(z,q) =  \sum_{r=0}^{\infty}B_{2r+1}^{(2n+1)}(q) \sin{((2r+1)z)}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3Bse_%7B2n%2B1%7D%28z%2Cq%29+%3D++%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DB_%7B2r%2B1%7D%5E%7B%282n%2B1%29%7D%28q%29+%5Csin%7B%28%282r%2B1%29z%29%7D%0A)

and

<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Bequation%7D%0Ase_%7B2n%2B1%7D%3D0+%5C%3B%5C%3B%5C%3B%5C%3B+%5Ctext%7Bat%7D+%5C%3B%5C%3B+y%3D0%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5Ctext%7Band%7D%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+%5Cdfrac%7Bdse_%7B2n%2B1%7D%7D%7Bdy%7D+%3D0%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+%5Ctext%7Bat%7D+%5C%3B%5C%3B+y%3D%5Cpi%5Cnonumber%0A%5Cend%7Bequation%7D" 
alt="\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{equation}
se_{2n+1}=0 \;\;\;\; \text{at} \;\; y=0\;\;\;\;\;\;\;\;\;\;\text{and}\;\;\;\;\;\;\;\; \dfrac{dse_{2n+1}}{dy} =0\;\;\;\;\;\; \text{at} \;\; y=\pi\nonumber
\end{equation}">

4. Even, **sine elliptic** function, with period ![\pi](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5Cpi).

![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;se_{2n+2}(z,q) =  \sum_{r=0}^{\infty}B_{2r+2}^{(2n+1)}(q) \sin{((2r+2)z)}
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3Bse_%7B2n%2B2%7D%28z%2Cq%29+%3D++%5Csum_%7Br%3D0%7D%5E%7B%5Cinfty%7DB_%7B2r%2B2%7D%5E%7B%282n%2B1%29%7D%28q%29+%5Csin%7B%28%282r%2B2%29z%29%7D%0A)


and

<img src=
"https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Bequation%7D%0Ase_%7B2n%2B2%7D%3D0%2C+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B+%5Ctext%7Bat%7D+%5C%3B%5C%3B+y%3D%5C%7B0%2C+%5Cpi%5C%7D%5Cnonumber%0A%5Cend%7Bequation%7D" 
alt="\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{equation}
se_{2n+2}=0, \;\;\;\;\;\; \text{at} \;\; y=\{0, \pi\}\nonumber
\end{equation}">

**cosine** eliptic functions cannot coexist with **sine elliptic** functions, since only one periodic solution is allows, but a second, non-periodic solution can be obtained from each of them. 


The pair of **characteristic numbers** ![a_{n}(q), b_{n}(q)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a_%7Bn%7D%28q%29%2C+b_%7Bn%7D%28q%29) are each associated with the mode n, and correspond to ![ce_{n}(q, z), se_{n}(q, z)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+ce_%7Bn%7D%28q%2C+z%29%2C+se_%7Bn%7D%28q%2C+z%29) respectively. Note that ![a_{n}(q) \neq b_{n}(q)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a_%7Bn%7D%28q%29+%5Cneq+b_%7Bn%7D%28q%29).



## Method of solution

Given the appropriate boundary conditions, there are many methods of solution ranging from regular (![|q| \ll 1](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%7Cq%7C+%5Cll+1)) and singular WKB (<img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%7Cq%7C+%5Cgg+1" 
alt="|q| \gg 1">) asymptotic expansions than help render approximations to the **cosine** and **sine** elliptic functions, as a function of varying parameter q. Mathieu's equation, a 2nd order ODE of Sturm-Louiville type, is an eigenvalue problem, and therefore here we rely on matrix methods to estimate the coefficients of the trigonometric series above that determine each of the periodic solutions to Mathieu Equation. The trigonometric series are indeed the Fourier series approximation of each of the Mathieu functions, and moreover, the series are absolutely and uniformly convergent since Mathieu functions also form a complete set in the same space as Fourier series do (Hilbert space). This also implies that integrable, periodic functions can be expanded in terms of Mathieu functions, which is particularly useful when Mathieu equation has a forcing term on the right hand side.


We illustrate the (matrix) solution method to estimate ![ce_{2n}(q, z)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+ce_%7B2n%7D%28q%2C+z%29)



![\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{eqnarray}
a_{2n}A_{0}^{(2n)} = &  qA_{2}^{(2n)}\nonumber\\
\left[a_{2n}-4 \right]A_{2}^{(2n)} = & q\left[ 2A_{0}^{(2n)} + A_{4}^{(2n)}\right], \;\;\;\;\;\;\text{and}\nonumber\\
\left[a_{2n} - 4r^2 \right]A^{(2n)}_{2r} = & q\left[A_{2r-2}^{(2n)} + A_{2r+2}^{(2n)} \right], \;\;\;\; \text{for}\;\; r\geq2\nonumber
\end{eqnarray}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Beqnarray%7D%0Aa_%7B2n%7DA_%7B0%7D%5E%7B%282n%29%7D+%3D+%26++qA_%7B2%7D%5E%7B%282n%29%7D%5Cnonumber%5C%5C%0A%5Cleft%5Ba_%7B2n%7D-4+%5Cright%5DA_%7B2%7D%5E%7B%282n%29%7D+%3D+%26+q%5Cleft%5B+2A_%7B0%7D%5E%7B%282n%29%7D+%2B+A_%7B4%7D%5E%7B%282n%29%7D%5Cright%5D%2C+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5Ctext%7Band%7D%5Cnonumber%5C%5C%0A%5Cleft%5Ba_%7B2n%7D+-+4r%5E2+%5Cright%5DA%5E%7B%282n%29%7D_%7B2r%7D+%3D+%26+q%5Cleft%5BA_%7B2r-2%7D%5E%7B%282n%29%7D+%2B+A_%7B2r%2B2%7D%5E%7B%282n%29%7D+%5Cright%5D%2C+%5C%3B%5C%3B%5C%3B%5C%3B+%5Ctext%7Bfor%7D%5C%3B%5C%3B+r%5Cgeq2%5Cnonumber%0A%5Cend%7Beqnarray%7D)


The above recurrence formula can be written as an eigenvalue problem, with an infinite tridiagonal matrix given by

<img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%5C%3B%0A%5Cbegin%7Bequation%7D%0A%5Cbegin%7Bpmatrix%7D%0A++++0+%26+%5Csqrt%7B2%7Dq+%26++%26++%26+%26+%26+%5C%5C%0A++++%5Csqrt%7B2%7Dq+%26+4+%26+q+%26++%26+%26+%26%5C%5C%0A++++%26+q+%26+16+%26+q+%26+%26+%26%5C%5C%0A++++%26+%26+q+%26+36+%26+q+%26+%26%5C%5C%0A++++%26+%26+%26+%5Cddots+%26+%5Cddots+%26+%5Cddots+%26+%5C%5C%0A++++%26+%26+%26+%26+q+%26+4r%5E2+%26+q%5C%5C%0A++++%26+%26+%26+%26+%26+%26+%5Cddots%0A%5Cend%7Bpmatrix%7D%5C%3B%0A%5Cbegin%7Bpmatrix%7D%0A%5Csqrt%7B2%7DA_%7B0%7D%5E%7B2n%7D%5C%5C%0AA_%7B2%7D%5E%7B%282n%29%7D%5C%5C%0AA_%7B4%7D%5E%7B%282n%29%7D%5C%5C%0AA_%7B6%7D%5E%7B%282n%29%7D%5C%5C%0A%5Cvdots%5C%5C%0AA_%7B2r%7D%5E%7B%282n%29%7D%5C%5C%0A%5Cvdots%0A%5Cend%7Bpmatrix%7D%0A%3D+a_%7B2n%7D%5C%3B%0A%5Cbegin%7Bpmatrix%7D%0A%5Csqrt%7B2%7DA_%7B0%7D%5E%7B2n%7D%5C%5C%0AA_%7B2%7D%5E%7B%282n%29%7D%5C%5C%0AA_%7B4%7D%5E%7B%282n%29%7D%5C%5C%0AA_%7B6%7D%5E%7B%282n%29%7D%5C%5C%0A%5Cvdots%5C%5C%0AA_%7B2r%7D%5E%7B%282n%29%7D%5C%5C%0A%5Cvdots%0A%5Cend%7Bpmatrix%7D%0A%5Cnonumber%0A%5Cend%7Bequation%7D" 
alt="\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;
\begin{equation}
\begin{pmatrix}
    0 & \sqrt{2}q &  &  & & & \\
    \sqrt{2}q & 4 & q &  & & &\\
    & q & 16 & q & & &\\
    & & q & 36 & q & &\\
    & & & \ddots & \ddots & \ddots & \\
    & & & & q & 4r^2 & q\\
    & & & & & & \ddots
\end{pmatrix}\;
\begin{pmatrix}
\sqrt{2}A_{0}^{2n}\\
A_{2}^{(2n)}\\
A_{4}^{(2n)}\\
A_{6}^{(2n)}\\
\vdots\\
A_{2r}^{(2n)}\\
\vdots
\end{pmatrix}
= a_{2n}\;
\begin{pmatrix}
\sqrt{2}A_{0}^{2n}\\
A_{2}^{(2n)}\\
A_{4}^{(2n)}\\
A_{6}^{(2n)}\\
\vdots\\
A_{2r}^{(2n)}\\
\vdots
\end{pmatrix}
\nonumber
\end{equation}">


When <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+q" 
alt="q"> is real, the matrix system is symmetric and thus all eigenvalues (**caracteristic values**) <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+a_%7B2n%7D" 
alt="a_{2n}"> are real. The **eigenvectors** of the matrix are the coefficients of the Fourier series and thus the method simultaneously determines the Fourier coefficients that approximate ![ce_{n}(q, z), se_{n}(q, z)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+ce_%7Bn%7D%28q%2C+z%29%2C+se_%7Bn%7D%28q%2C+z%29). 
When <img src=
"https://render.githubusercontent.com/render/math?math=%5Ctextstyle+q%3Dis%2C+%5C%3B+%5Ctext%7Bwith%7D%5C%3Bs%3E0+" 
alt="q=is, \; \text{with}\;s>0 "> real, the system is not Hermitian and branching occurs[[3]](#3)[[4]](#4). In such case, the eigenvalues are no longer **ordered** and additional steps are required to provide the correct **characteristic values** and **Fourier coefficients**.

The advantage of the matrix system over the perturbation approach, is the well documented python libraries for determining eigenvalues and eigenvectors, that allow the fast and efficient approximation to Mathieu functions of the first kind ![ce_{n}(q, z), se_{n}(q, z)](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+ce_%7Bn%7D%28q%2C+z%29%2C+se_%7Bn%7D%28q%2C+z%29). The greater the order of the matrix, the better the approximation to these functions.


## References
<a id="1">[1]</a> 
Arscott, F. M. (1964). 
Periodic differential equations: an introduction to Mathieu, Lame, and allied functions.
International series of monographs of pure and applied mathiematics.
Volume 66.
Pergamon Press.

<a id="2">[2]</a>
Whittaker, E. T. (1903).
A course in modern analysis.
Cambridge University Press.

<a id="3">[3]</a>
Hunter, C. and Guerrieri, B. (1981).
The eigenvalues of Mathieu's equation and their branch points.
Studies in Applied Mathematics 64(2), pages 113-142
Wiley Online Library.

<a id="4">[4]</a>
Zeiner, C. H. et al. (2012)
Mathieu functions of purely imaginary parameters.
Journal of Computational and Applied Mathematics.
Elsevier.








