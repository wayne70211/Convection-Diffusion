# Convection-Diffusion

Consider the convection-diffusion equation <br>

<p align="center"> 
<img src="https://latex.codecogs.com/gif.latex?\large&space;\frac{\partial&space;T}{\partial&space;t}&space;&plus;&space;u&space;\frac{\partial&space;T}{\partial&space;x}&space;=&space;\alpha&space;\frac{\partial^2&space;T}{\partial&space;x^2}" title="\large \frac{\partial T}{\partial t} + u \frac{\partial T}{\partial x} = \alpha \frac{\partial^2 T}{\partial x^2}" />
</p>

With the boundary conditions<br>
<p align="center"> 
<img src="https://latex.codecogs.com/gif.latex?\large&space;T(0,t)=T(1,t)=0" title="\large T(0,t)=T(1,t)=0" />
</p>

### Part 1. Pure convection (α = 0) 

Consider the following initial profile <br>
<p align="center"> 
<img src="https://latex.codecogs.com/gif.latex?\large&space;T(x,0)=\left\{\begin{matrix}&space;1-(10x-1)^2&space;&&space;at&space;\;&space;\:&space;0<x\leq&space;0.2&space;\\&space;0&space;&&space;at&space;\;&space;\:&space;0.2<x\leq&space;1&space;\end{matrix}\right." title="\large T(x,0)=\left\{\begin{matrix} 1-(10x-1)^2 & at \; \: 0<x\leq 0.2 \\ 0 & at \; \: 0.2<x\leq 1 \end{matrix}\right." />
</p>

The exact solution is <br>

<p align="center"> 
<img src="https://latex.codecogs.com/gif.latex?\large&space;T(x,t)=\left\{\begin{matrix}&space;1-(10(x-ut)-1)^2&space;&&space;at&space;\;&space;\:&space;0<x-ut\leq&space;0.2&space;\\&space;0&space;&&space;at&space;\;&space;\:&space;0.2<x-ut\leq&space;1&space;\end{matrix}\right." title="\large T(x,t)=\left\{\begin{matrix} 1-(10(x-ut)-1)^2 & at \; \: 0<x-ut\leq 0.2 \\ 0 & at \; \: 0.2<x-ut\leq 1 \end{matrix}\right." />
</p>

#### (1) Explicit Euler time advancement and second-order central difference for the spatial derivative.
Discretization:<br>

<p align="center"> 
<img src="https://latex.codecogs.com/gif.latex?\large&space;\frac{T^{n&plus;1}_{j}-T^{n}_{j}}{\Delta&space;t}&space;&plus;&space;u&space;\frac{T^{n}_{j&plus;1}-T^{n}_{j-1}}{2\Delta&space;x}&space;=&space;0" title="\large \frac{T^{n+1}_{j}-T^{n}_{j}}{\Delta t} + u \frac{T^{n}_{j+1}-T^{n}_{j-1}}{2\Delta x} = 0" />
</p>

From modified wave number analysis **Explicit Euler is unstable**.

#### (2)	Leapfrog time advancement and the second-order central difference for the spatial derivative. 
Discretization:<br>

<p align="center"> 
<img src="https://latex.codecogs.com/gif.latex?\large&space;\frac{T^{n&plus;1}_{j}-T^{n-1}_{j}}{2\Delta&space;t}&space;&plus;&space;u&space;\frac{T^{n}_{j&plus;1}-T^{n}_{j-1}}{2\Delta&space;x}&space;=&space;0" title="\large \frac{T^{n+1}_{j}-T^{n-1}_{j}}{2\Delta t} + u \frac{T^{n}_{j+1}-T^{n}_{j-1}}{2\Delta x} = 0" /></p>
 
From modified wave number analysis **Leapfrog is stable**.

## Simulation Result
Let <img src="https://latex.codecogs.com/gif.latex?\large&space;u=0.08" title="\large u=0.08" />. 
<p align="center"> 
<img src="https://github.com/wayne70211/Convection-Diffusion/blob/master/Simulation.gif">
</p>