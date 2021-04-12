### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ d54d301e-060e-4c14-b5c8-1fca26dd84f9
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="Test"), 
			Pkg.PackageSpec(name="LinearAlgebra"), 
			Pkg.PackageSpec(name="PlutoUI"), 
			Pkg.PackageSpec(name="Plots")
			])
	using Test
	using PlutoUI
	using LinearAlgebra
	using Plots
end

# ╔═╡ 606540b0-9846-11eb-2a09-8112c5081854
md"""
# SIM 203
"""

# ╔═╡ 57160c99-cccc-4aca-99fd-df9356da76a1
PlutoUI.TableOfContents(;depth=2)

# ╔═╡ 28280ba3-789f-40ec-b731-cbc43334b839
md"""!!! danger "Submission instructions"
	**Due date:** April 26, 2021
	
	Send your solutions to **luiz.maltez-faria@inria.fr**
"""

# ╔═╡ 6ea0197b-5af3-4116-b214-a27f28508c33
md"""## Introduction
This homework is intended to give you a *hands-on* experience on coding some of the fundamental numerical linear algebra algorithms that we have covered in the lectures so far. You are encouraged to consult the course material available [here](https://perso.ensta-paris.fr/~mbonnet/aln.pdf). 

Most of the questions are of practical nature; that is, I will give you guidelines and you are expected to **write code** which implements a given algorithm. I will also provide some basic tests that your code should pass, but you may want to write your own tests as you go along. Also, *please comment your code* so that I can understand the reasoning behind what you wrote! 

Some questions require a written explanation (beyond code that is), and you should create a cell with your answer. Remember that you can use `latex` syntax as your normally would by e.g. enclosing your text between the dollar sign `$`. Check some of the cells on this notebook for examples. 

Your code will be implemented in [Julia](https://docs.julialang.org/en/v1/), and the page you are viewing was created using [Pluto.jl](https://github.com/fonsp/Pluto.jl). Although you do not necessarily need to install/use Pluto.jl in your homework, doing so is highly recommended since it will allow you to navigate this notebook interactively (instead of viewing it as a static html page). [This video](https://www.youtube.com/embed/OOjKEgbt8AI) gives some basic instructions on how to install Julia and Pluto.jl (ignore the part at the end regarding the submission). If you run into a problem let me know. 

Once you are done installing *Julia* (and *Pluto.jl*), you should probably spend some time playing around with *Julia* if you have never used it before. Its syntax is very similar to *Matlab*, so you should be up and running in little time. Two important differences are:
* Indexing arrays is done through square brackets. That is, if `A` is `Matrix`, then use `A[i,j]` to index the `(i,j)` entry of `A` instead of `A(i,k)` as in *Matlab*.
* *Julia* uses a *passing by reference* convention for arguments in a function, and arrays are not copied by default when using the `=` operator. This means e.g. that if `A` is a `Matrix`, then `B=A` does not create a new copy of `A` and the variable `B` references the same memory; in particular mutating `B` will modify `A`! If you want a copy, use instead `B=copy(A)`.

Here is a [cheatsheet](https://cheatsheets.quantecon.org) summarizing some basic syntactic differences between *Julia*, *Matlab* and *Python*. Once you are comfortable with *Julia*, move on to the next section. 
"""

# ╔═╡ 8c26975d-ec0c-423c-9644-daf031c0aabb
md""" ## Part 1: QR factorization and least squares

Throughout this section we will assume that $A \in \mathbb{C}^{m\times n}$ with $m \geq n$, and $\text{rank}(A) = n$ (i.e. $A$ has full rank). This is a simplyfying assumption which is not necessarily needed in many of the algorithms (such as QR). 

The main goal of this section is to write a program that can solve the following least squares problem:
```math
	\text{least squares} \rightarrow \underset{x \in \mathbb{C}^n}{\mathrm{argmin}}|| Ax - b ||_2
```
Since we assume that $A$ has full rank, this problem has a unique solution. We will proceed to solve this problem by performing a QR factorization of $A$ through Householder reflectors, as described next. """

# ╔═╡ c707bd38-32f2-4adb-bdaf-8c8325d40ab9
md"""### Householder reflectors
Recall the definition of the Householder reflector:
```math
F(u) = I - 2 \frac{uu^*}{u^* u},
```
where $I$ is the identity matrix an $u$ is a complex vector. The notation $(\cdot)^*$ is used for the adjoint. Your first task is to code precisely this:
"""

# ╔═╡ 06ddb55b-e7a7-4ee8-91dd-171d5665df5e
function householder_matrix(v)
	return missing
end

# ╔═╡ d48b5435-ce1d-4789-a33a-91f405439138
md"""If your implementation is correct, it should pass the tests I have provided. Once they do, the cell below should change."""

# ╔═╡ 5a6887a8-e651-4e4b-96c6-d1d7058cd26d
md"""
Now, at each step of the *QR* factorization, we need to compute a Householder vector $u$ given a vector $x$ so that 
```math 
F(u)x = -\textrm{sign}(x_1)||x||e_1
```
where $x_1$ is the first entry of $x$ and $e_1$ is the cartesian basis vector. This is your second task:
"""

# ╔═╡ 1186d2ae-f51e-4237-a14f-e4275b5b304c
function householder_vector(x)
	return missing
end

# ╔═╡ a32a7c6e-2824-42ee-98dc-05781f5e7ffb
# write some code to test your function

# ╔═╡ ef865aba-8f80-4f38-9560-feecdb9013a7
md"""### *Naive* QR algorithm
We can now code a version of the Householder QR factorization similar to the one on your [lecture notes](https://perso.ensta-paris.fr/~mbonnet/aln.pdf). One difference will be that instead of simply reducing the matrix $A$ to an upper triangular form using orthogonal triangularization, we will actually form the matrix $Q$ in the process. Yes, that is inneficient, but we will do that anyways. 

Recall that the main idea is to gradually convert $A$ into an upper triangular matrix $R$ as follows:
```math
\underbrace{F(u_n)F(u_{n-1}) \ldots F(u_2)F(u_1)}_{Q^*} A = R.
```
Remember: the matrices $F(u_k)$ act only on the $k$ through $n$ columns of the matrix to its right. On your lectures notes that was achieved by some *zero padding* of the vectors $u_k$. Here we will only compute the non-zero part of these matrices. 
"""

# ╔═╡ 65523ce4-ebb8-4fb2-9869-61890e0b7e2d
function householder_qr(A)
	R   = copy(A)
    m,n = size(A)
	kmax = m>n ? n : m-1
    Qt   = diagm(0=>ones(ComplexF64,m)) # initilize Q^* to m×m identity matrix
    for k in 1:kmax
		# write your code here. 
		# Don't forget to update Qt as well as R
    end
    Q = adjoint(Qt)
    return Q,R
end

# ╔═╡ 665919fd-592f-494a-9d17-adde3434a7e0
md"""
If you got this far you have succeeded in implementing a $QR$ factorization in *Julia*! 
"""

# ╔═╡ 7c40b5d3-8afd-4648-91ef-ec6235ec376e
md"""### A compact QR algorithm
Forming the matrix $Q$ requires some extra computation (and extra memory), and having $Q$ explicitly formed is often not necessary. In what follows we will dive a little into a more efficient *compact* representation of the *QR* factorization. The idea is to use the upper triangular part of the input matrix $A$ to store $R$, and the lower triangular part (diagonal excluded) to store the Householder vectors $u_k$. We will normalize $u_k$ so that its first value is always one, meaning only $m-k$ values are needed to define it instead of $m-k+1$ values. This *trick* means $u_k$ will fit in the lower triangular part (diagonal excluded) of the $k$-th column of $A$, and therefore we should be able to perform our QR factorization using *only the memory* allocated in $A$. This will require some significant rewrite. In particular

* The `householder_vector` function should normalize the vector so that the first entry is `1`.
* We need to write a function to *apply* $H(u)$ to a matrix given the vector $u$ that defines it.  
* Some care has to be taken regarding the order of the operations so that we do not overwrite memory that is still needed for other computations. 

The function `householder_vector_normalized!(x)` shown below modifies `x[2:end]` to store the `2:end` entries of the Householder vector. Note that the first entry is not modified since it is (implicitly) one. This memory space will be reserved to the diagonal part of $R$.
"""

# ╔═╡ e24e5718-1a1a-47bc-8f25-3268ce9a5826
function householder_vector_normalized!(x)
    # norm of x
    d2 = x' * x
    d  = sqrt(d2)
    # normalized reflection vector
    scale = 1/(x[1] + sign(x[1])*d)
	m     = length(x)
	for k in 2:m
		x[k] *= scale
	end
    return x
end

# ╔═╡ 13c5b2d6-7f30-4afd-8a08-3dfe4f2db184
md"""
The `apply_hh_matrix` function below uses the special structure of $H(u)$ to efficiently computes $H(u)A$ where $u \in \mathbb{C}^m$ and $A \in \mathbb{C}^{m\times n}$ **without allocating any memory**: 
"""

# ╔═╡ 355bd86a-a278-4ade-96a7-2a18e0c8f7db
function apply_hh_matrix!(A,u)
	m,n = size(A)
	@assert length(u) == m
	beta = 2/(u'*u)
	for j in 1:n
		col   = @views A[:,j]
		coeff = u'*col*beta
		for i in 1:m
			col[i] -= coeff*u[i]
		end
	end
	return A
end

# ╔═╡ e24068de-07a2-41df-a686-01b3e21ca6b0
md"""!!! tip
	In `julia`, *slicing* an array (e.g. `x[5:10]`) creates a new array with the 	sliced entries. This means memory is allocated to store the new object. The `@views` macro used in the `apply_hh_matrix!` function above avoids such allocations by creating a `view` of the array which uses the same memory. 
"""

# ╔═╡ 4bbfd29d-832b-4e17-8f59-25463466fa37
md"Your turn now:"

# ╔═╡ b7dc76de-d6df-414e-bc76-fa9c5e2af727
function apply_hh_matrix_normalized!(A,u)
	# your code here. 
	# A good starting point is the code from apply_hh_matrix!
	return A
end

# ╔═╡ ae805230-574e-4818-bdbc-1dda19fa33b0
function householder_qr_compact!(A)
    m,n = size(A)
    for k=1:n
		# your code here
    end
    return A
end

# ╔═╡ a0af5567-02c6-4c3a-a1f2-649c622d1319
let
	m,n = 200,50
	A = rand(ComplexF64,m,n)
	t1 = @elapsed householder_qr(A)
	t2 = @elapsed householder_qr_compact!(A)
	m1 = @allocated householder_qr(A)
	m2 = @allocated householder_qr_compact!(A)
	md"""You case you got this far, the table below will give you a rough idea of execution time of the *naive* and *compact* implementations of your algorithm for a matrix of $A \in \mathbb{C}^{m\times n}$ with $ m=$m $ and $ n = $n $. The dynamically allocated memory is also shown. **These numbers only make sense if your algorithm is working**. 
	
Version | Time (s) | Memory (bytes)
:------------ | :-------------: | :----------:
naive   | $t1 | $m1
compact | $t2 | $m2
	"""
end

# ╔═╡ a1cfd851-c006-4506-be18-7571927fb5ef
md"""
You are now (almost) ready to solve the original least squares problem!
"""

# ╔═╡ f89404e4-9287-4014-bd2d-93f6833a1985
md"""### Least squares problem
"""

# ╔═╡ 8efa2b75-f85a-4be9-85cd-6cccd193bded
md"""Now that everthing is in place, solve the least squares problem:
"""

# ╔═╡ 2eb53b12-d20d-411c-9da1-96a0b6ecac39
# solve Rx = b
function solve_uppertriangular(R,b)
    m,n = size(R)
    @assert m==n
    @assert n==length(b)
    x = copy(b)
    for i in n:-1:1
        for j in i+1:n
            x[i] -= R[i,j]*x[j]
        end
        x[i] /= R[i,i]
    end
    return x
end

# ╔═╡ 7c2e9206-8403-4eab-a39e-6625945dd0c1
function leastsquares_qr(A,b)
	return missing
end

# ╔═╡ 1ab23f09-85e4-48fb-be05-2afe81717fbc
md"""## Part 2: Hessenberg decomposition and Arnoldi iterations

In this section we will focus on an important matrix decomposition: *the Hessenberg decomposition*. In what follows we assume that $A$ is a square and invertible matrix of size $\mathbb{R}^{m\times m}$, and seek a decomposition of the form 
```math
A = Q H Q^*
```
with $Q$ orthogonal and $H$ Hessenberg (i.e `H[i,j]=0` for `i>j+1`). Many eigenvalue algorithms in linear algebra proceed by first reducing $A$ to its Hessenberg form, and then computing the eigenvalues of $H$ efficiently; this makes the reduction to Hessenberg form an important tool in numerical linear algebra. 

One idea to produce a Hessenberg decomposition is to proceed in a way similar to the QR algorithm using Householder reflections; that is, you *gradually* transform the matrix $A$ into a Hessenberg matrix through left and right multiplication by *simple* unitary matrices. This is the goal of your next question:
"""

# ╔═╡ f7b8eb66-51a6-49be-b755-c6697e72bbaa
md"""
It can be shown that recursively applying the idea from the previous question leads to a Hessenberg decomposition of $A$. That is
```math
    \underbrace{F(u_n) F(u_{n-1}) \ldots F(u_2) F(u_1)}_{Q^*} \ A \ \underbrace{F(u_1)^* F(u_{2})^* \ldots F(u_{n-1})^* F(u_n)^*}_{Q}= H
```
where each $F(u_k)$ puts zeros on the `k+2:m` entries of the `k`-th column. 

This algorithm has $\mathcal{O}(n^3)$ complexity ---- for very large matrices this is often too expensive. The idea of Arnoldi iterations is to generate a *partial* Hessenberg decomposition of $A$. That is, letting $Q_n$ denote the $m \times n$ matrix  obtained from the first $n$ columns of $Q$, and $H_n$ denote the $(n+1) \times n$ matrix obtained from the first $n+1$ rows of $H$ and its first $n$ columns, we seek a decomposition of the form 
```math
A Q_n = Q_{n+1} H_n
```
The Hessenberg nature of $H$ makes for a recursive interpretation of the above equality:
"""

# ╔═╡ 26b12556-a9b8-4de7-8905-b33395ed7738
function arnoldi(A,n,b=A[:,1])
    T   = eltype(A)
    m   = size(A,1)
    Q   = zeros(T,m,n+1) # Q_{n+1}
    H   = zeros(T,n+1,n) # H_n
    Q[:,1] = b / norm(b)
    for k in 1:n
		# write here your code 
    end
    return H,Q
end

# ╔═╡ 391ba967-0b54-4827-aac1-0382b095bfe6
md"""## Part 3: GMRES

We are now in a position to better understand the inner working of GMRES. Recall that at step $n$ of the algorithm, GMRES finds the following approximate solution $x_n$:
```math
	x_n = \underset{x \in \mathcal{K}^n}{\mathrm{argmin}} || Ax - b||
```
where the Krylov subspace $\mathcal{K}_n$ is:
```math
\mathcal{K}_n = \left\{b,Ab,\ldots,A^{n-1}b \right\}
```

Since $x_n \in \mathcal{K}_n$, it makes sense to consider a basis for $\mathcal{K}_n$ and seek instead for the *coefficients* of $x_n$ in that basis. Letting $Q_n$ be an $m \times n$ matrix whose columns form an orthogonal basis for $\mathcal{K}_n$, we may therefore set $x_n = Q_ny_n$, where $y_n$ is to be found. The equivalent problem is then
```math
\begin{align}
	y_n = \underset{y \in \mathbb{C}^n}{\mathrm{argmin}} || AQ_ny - b||
\end{align}
```
And here is where the magic happens! Because of the particular choice of the space $\mathcal{K}_n$, we can use the Arnoldi iteration algorithm to compute an orthogonal basis for $\mathcal{K}_n$, so that the equation above simplifies even further:
"""

# ╔═╡ f48c878b-3b00-4778-8b82-e120acd5e188
function gmres(A,b;verbose=false,maxiter=100,tol=1e-8)
    T  = eltype(A) # element type. For you examples this is ComplexF64
    m  = size(A,1)
	β  = norm(b)
    q1 = b / β
    xn = zeros(T,m) # initialize an output vector 
    hist = []       # for keeping a long of the residue
    for n in 1:maxiter
		# fill in this part to update xn, then
		# store the residue at iteration n in the res variable
		res = norm(A*xn - b)
        verbose && println("iteration $n \t residual $res") 
        push!(hist,res)
        res < tol && break
    end
    return xn,hist
end

# ╔═╡ 06c05919-b335-4bda-8b14-8bfdb66d7975
md"""### A note on further optimizations
If the cell above has turned blue, you have managed to implement the GMRES algorithm! There are some further algorithmic optimizations which we have left on the table. In particular, solving the least-squares problem for the Hessenberg matrix $H_n$ can be done more efficiently using [*Gives rotations*](https://en.wikipedia.org/wiki/Givens_rotation), and the $QR$ factorization of $H_n$ can be used to cheaply compute the $QR$ factorization of $H_{n+1}$ at $\mathcal{O}(n)$ cost. In practice these may be important optimizations when the number of iterations needed for convergence is "large" (where "large" here is intentionally vague).  These optimizations go beyond the scope of this homework. 
"""

# ╔═╡ 40102829-a65b-4201-847e-f2c10abb27a5
md"""### Convergence of GMRES

As seen in class, the convergence of GMRES is a subtle question which depends both on the condition number of its matrix of eigenvectors and on the distribution of the eigenvalues.

Below we construct two matrices of similar condition number but very different spectrum distribution, and show that GMRES behaves very differently. This example was extracted from lecture 35 of [this book](https://people.maths.ox.ac.uk/trefethen/text.html). The first matrix $A$ is given by 
```julia
A = 2*I + 0.5*randn(m,m)/sqrt(m)
```
Its eigenvalues are (approximately) uniformly distributed on a disk of radius $0.5$ centered at $z_0 = 2 + 0i$, as can be observed below:
"""

# ╔═╡ 4a1f0c7d-9fc3-4e08-97be-4eb3b8e33705
begin
	n = m = 500
	A = 2*I + 0.5*randn(m,n)/sqrt(m);
	F = eigen(A)
	scatter(F.values,label="",aspect_ratio=1,title="spectrum of A")
end

# ╔═╡ 4790f8a5-3cb3-4ccf-a8b4-85de126f8d31
md"""The matrix $B$ is defined using:
```julia
d = [-2 + 2 * sin(k*π/(m-1)) + im*cos(k*π/(m-1)) for k in 0:m-1]
B = A + 2*diagm(d)
```
and its eigenvalues wrap around the origin:
"""

# ╔═╡ 48d0de6e-ac83-4b55-94a4-1361e6f2329e
begin
	d = [-2 + 2 * sin(k*π/(m-1)) + im*cos(k*π/(m-1)) for k in 0:m-1]
	B = A + 2*diagm(d)
	FB = eigen(B)
	scatter(FB.values,label="",aspect_ratio=1,title="spectrum of B")
end

# ╔═╡ 254bfbc4-117b-4518-a188-c19a7c123035
md"""## Utility functions
"""

# ╔═╡ 363968d9-23c9-45e7-b47d-5d7e8ad5648f
question(text,number="") = Markdown.MD(Markdown.Admonition("tip", "Question $number", [text]))

# ╔═╡ 2cded3aa-8661-4477-bcc8-ac375e01622c
question(md"Complete the `householder_matrix` function below to compute the Householder reflector.",1)

# ╔═╡ 4720b6c3-d9bf-4a46-9d38-2cee67e35142
let
	text = md"""Complete the `householder_vector` function below to implement the  	  *Householder vector* $u$ given an input vector $x$. Check that the property 
	```math 
		F(u)x = -\textrm{sign}(x_1)||x||e_1
	```
	is satisfied.
	"""
	question(text,2)
end

# ╔═╡ 9cfa6c75-cb75-470d-b037-2f01c2f22ae9
question(md"Complete the code below to implement a modified version of Algorithm 5 of the lecture notes where the matrices Q and R are explicitly constructed. You should return the *full QR*, so that $Q \in \mathbb{C}^{m\times m}$ and $R \in \mathbb{C}^{m \times n}$",3)

# ╔═╡ ee32ef8b-3b5f-4bd8-bafc-588166899782
question(md"Describe in words and formulas the steps in the code for `apply_hh_matrix!` above.",4)

# ╔═╡ 0f634e30-4f59-48d4-9843-8bbb1e244db5
question(md"""Modify `apply_hh_matrix!` so that it implicitly assumes `u[1]==1`. We will call the new function `apply_hh_matrix_normalized!`. Note that this function should **not** use the value of `u[1]`.
""",5)

# ╔═╡ 5ed2c935-22d6-4e3c-ab2a-a52dd8ee462b
question(md"""
Complete the code below to perform a compact QR factorization in place. There should be no (dynamic) memory allocated.
""",6)

# ╔═╡ dbf408ba-bf6e-440c-8855-3f7b7f254d6e
question(md"""Reformulate the least squares problem 
```math
\begin{align}
\underset{x \in \mathbb{R}^n}{\mathrm{argmin}}||Ax - b||_2
\end{align}
```
using the QR factorization of A, and complete the function `least_squares_qr` below to calculate $x$ given $A$ and $b$. Use the provided function `solve_uppertriangular(R,b)` below
to find the solution of $Rx=b$, where $R$ is a (square) upper-triangular matrix, and your version of  `householder_qr` for the QR decomposition. 
	""",7)

# ╔═╡ 6217aa0d-46e1-4dd6-a1e5-5eb2696757e4
question(md"""Suggest a way to solve the least squares problem using the `householder_qr_compact!` factorization. What method(s) do you need to implement?
""",8)

# ╔═╡ 1ab7ca93-8c5b-428c-9344-a12d087cc489
question(md"""Let $A \in \mathbb{R}^{n\times n}$. Find a unitary matrix $P$ such that $P A P^*$ has zero below the second row on the first column.""",9) 

# ╔═╡ 1e4cda7d-68c7-498a-8bb8-d8642f67f344
question(md"""Write down a recursive formula for computing the entries of $H_n$ and $Q_{n+1}$. You may initialize $\mathbf{q}_1$ to the first column of $A$ (normalized). """,10)

# ╔═╡ a8a76616-3afb-4370-b213-2ab0f950ecf6
question(md"""Complete the algorithm below to compute $H_n$ and $Q_{n+1}$ given $A$ and $n$.""",11)

# ╔═╡ 2f62830f-5b5c-4478-bc33-3b559bd5a2a3
question(md"""
Using the Arnoldi iteration algorithm of Section 2, rewrite the least squares problem for $y_n$ in terms of the Hessenberg matrices $H_n$. Then complete the code below for solving a linear system $Ax=b$ using the GMRES algorithm. Do not worry about making your code efficient: the goal here is understand how GMRES works. 
"""
,12)

# ╔═╡ 732154d4-bf4e-40db-aec1-713208032fda
question(md"""
Using your GMRES algorithm, explore the convergence of the residue when solving both $Ax = b$ and $Bx=b$. How many iterations are needed to converge to a tolerance of $10^{-8}$ in both cases? Can you explain the favorable behavior of GMRES in the first example? Can you estimate the convergence rate in that case?
""",14)

# ╔═╡ c65c438b-d97f-4857-9897-20180cc53d25
answer(text=md"Write your answer here!
") = Markdown.MD(Markdown.Admonition("warning", "Answer", [text]))

# ╔═╡ 21f916fe-21b9-43c1-b9d1-a8e541d24e43
answer(md"""
Write your answer here!
""")

# ╔═╡ f31d5edd-722e-4c55-9ef3-190b4996fd29
answer(md"""
Write your explanation here (and don't forget to complete the code).
""")

# ╔═╡ 165cfb6b-63f6-4857-849b-559ba9221768
answer(md"""
Write your answer here!
""")

# ╔═╡ ec546b10-2ab0-4d79-b81d-80532742abc2
answer(md"""
Write your answer here!	
""")

# ╔═╡ a1387fac-d1b6-46b3-86aa-8bd0ad868508
answer()

# ╔═╡ be19ddeb-cfbf-4904-8aa7-62669055f07c
answer(md"""
Answer here, and don't forget to complete the code. 
""")

# ╔═╡ 2149e010-ccc6-46a4-a3eb-d967fb99fbf1
answer()

# ╔═╡ 78e2208d-f9f8-44f7-806f-4590db85fe1a
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]))

# ╔═╡ 989434b2-0f36-4998-b213-5ba6d2ce7166
hint(md"You can use `I` for representing an identity matrix in *Julia*. It is similar to *Matlab*'s `eye`, but it does not form a matrix. Check also `diagm` for a more general method. If you are running this notebook interactively through *Pluto.jl*, check the live docs!")

# ╔═╡ 5fc13f1f-5818-4ebd-8531-fa5900c47bc9
hint(md"Think about the operation $Q^*b$ and how to implement that when only the normalized Householder vectors have been stored.")

# ╔═╡ 06f60aa7-2348-4313-a0c7-b748539659d8
hint(md"""Consider the $P$ to be the Householder reflector constructed using `A[2:m,1]`. By construction $PA$ will have zeroes on all entries below the second row on the first column. What happens when you multiply on the right by $P^*$?""")

# ╔═╡ a7aa6fe4-3ee2-4dff-a792-4289d4fca83c
keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]))

# ╔═╡ e25e883e-e0f9-4848-adfd-24f155ecf902
correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("note", "Got it!", [text]))

# ╔═╡ 21f882ba-fa33-4695-9324-b058ba64d258
note(text) = Markdown.MD(Markdown.Admonition("note", "Note", [text]))

# ╔═╡ 8c0961b4-4982-4355-9f7a-6aa37e342f3c
note(md"The cell below will add the dependencies that you need for this notebook. This may take some time (say a few minutes) the first time you run it since some packages will be downloaded and precompiled.")

# ╔═╡ 92791a5f-02f8-461e-9bb1-4e1077f24be9
note(md"""In *Julia*, functions whose name end with `!` typically mutate their arguments. This is merely a convention, but you stick with it. """)

# ╔═╡ 74097a18-b4fe-4ec8-b47c-b279fce28598
note(md"""
In case your `householder_qr` function does not work, you may use `LinearAlgebra.qr` for perfoming the $QR$ factorization. Check its documentation for more details on how to use it.
""")

# ╔═╡ 167a446c-24b7-414c-9bdf-bae7c90df2ab
warning(text) = Markdown.MD(Markdown.Admonition("warning", "Warning", [text]))

# ╔═╡ 95fd022d-3797-4f93-94cf-c7d5a928261c
warning(md"To keep things simple, we will implement suboptimal versions of the algorithms studied, focusing mostly on algorithmic aspects.")

# ╔═╡ 7e94c5d7-ba40-4850-9602-18fe63381ee8
function check_answer(f)
	try
		t = f()
		correct()
	catch E
		keep_working()
	end
end

# ╔═╡ bb7b2a43-4e93-4bfb-a41f-65ff6376babd
check_answer() do
	@testset "Householder matrix" begin
		v = rand(10)
		F = householder_matrix(v)
		@test F*v ≈ -v # reflection
		@test F*F' ≈ F' * F ≈ I # unitary
	end
end

# ╔═╡ 15410cfb-db47-48e0-a7e0-49a061aaaec5
check_answer() do 
	A    = rand(5,5) .- 0.5 # a random matrix for testing
	a1   = A[:,1]
	u1   = householder_vector(a1)
	F    = householder_matrix(u1)
	A1   = F*A
	@testset "Householder vector" begin
		@test A1[1,1] ≈ -sign(a1[1])*norm(a1)
		@test norm(A1[2:end,1],Inf) < 1e-10
	end
end


# ╔═╡ 9fc30c4d-e5c2-4359-9bec-e93ca5876763
check_answer() do 
	@testset "Householder qr" begin
		for (m,n) in [(10,5),(10,10)]
			A   = rand(ComplexF64,m,n)
			Q,R = householder_qr(A)
			@test size(Q) == (m,m)
			@test size(R) == (m,n)
			# test ortogonality
			@test Q'*Q ≈ I
			# test factorization
			@test Q*R  ≈ A
			# test that R is lower triangular
			@test norm(tril(R,-1),Inf) < 1e-15
		end
	end
end

# ╔═╡ e15fdd6e-ed3e-49bf-a3a7-7530cda0d479
check_answer() do 
	A    = rand(10,5)
	x    = rand(10)
	x[1] = 1
	ex = (I - 2*x*x'/(x'*x))*A
	x[1] = rand()
	apply_hh_matrix_normalized!(A,x)
	@testset "Appling householder matrix normalized" begin
		@test ex ≈ A
		#mem = @allocated apply_hh_matrix_normalized!(A,x)
		#@test mem == 0
	end
end

# ╔═╡ 2e865221-74b1-4023-962c-77203686cf03
check_answer() do 
	A = rand(10,5)
	Q,R = householder_qr(A)
	householder_qr_compact!(A)
	@testset "Compact QR" begin
		@test triu(A) ≈ R
		mem = @allocated householder_qr_compact!(A)
		@test mem == 0
	end
end

# ╔═╡ c727e963-6afb-41f5-9de8-38026ab126e0
check_answer() do
	@testset "Least squares QR" begin
		m,n = 10,5
		A = rand(ComplexF64,10,5)
		b = rand(10)
		x = leastsquares_qr(A,b)
		xe = A\b
		@test norm(x-xe,Inf) < 1e-10
	end
end

# ╔═╡ c4153dcf-6ff6-419d-bb5b-efe91656ea28
check_answer() do
	@testset "Arnoldi iteration" begin
		m = 20
		n = 5
		A = rand(ComplexF64,m,m)
		H,Q = arnoldi(A,n,rand(m))
		Qn  = Q[:,1:n]
		@test Q' * Q ≈ I
		@test norm(tril(H,-2),Inf) < 1e-10
		@test A*Qn ≈ Q*H
	end
end

# ╔═╡ e5b9f7d1-9af4-4534-80c2-6019430da0c9
check_answer() do
	@testset "GMRES test" begin
		n      = 100
		A      = 2*I + 0.5*randn(n,n)/sqrt(n) + im*I
		b      = rand(ComplexF64,n)
		xe     = A\b
		x,hist = gmres(A,b)
		@test norm(xe - x,Inf)/norm(xe,Inf) < 1e-5
	end
end

# ╔═╡ 8b8a87c1-871d-4fab-bf1a-6b3d8bc9a0fd
function show_tests(f)
	with_terminal() do
		try
			f()
		catch
		end
	end
end

# ╔═╡ 297f8790-3973-462e-b14a-bff39f399e8e
show_tests() do 
	A = rand(10,5)
	x = rand(10)
	exact = (I - 2*x*x'/(x'*x))*A
	apply_hh_matrix!(A,x)
	@testset "Appling householder matrix" begin
		@test exact ≈ A
		mem = @allocated apply_hh_matrix!(A,x)
		@test mem == 0
	end
end

# ╔═╡ 9ad83c33-f843-4d63-b0f3-1a0401068c80
show_tests() do 
	try
		@testset "Upper triangular solver" begin
			A = rand(10,10)
			b = rand(10)
			R = triu(A)
			xe = R\b
			x = solve_uppertriangular(R,b)
			@test norm(x - xe,Inf) < 1e-10
		end
	catch
		@warn "Exception encoutered"
	end
end

# ╔═╡ Cell order:
# ╟─606540b0-9846-11eb-2a09-8112c5081854
# ╟─57160c99-cccc-4aca-99fd-df9356da76a1
# ╟─8c0961b4-4982-4355-9f7a-6aa37e342f3c
# ╠═d54d301e-060e-4c14-b5c8-1fca26dd84f9
# ╟─28280ba3-789f-40ec-b731-cbc43334b839
# ╟─6ea0197b-5af3-4116-b214-a27f28508c33
# ╟─95fd022d-3797-4f93-94cf-c7d5a928261c
# ╟─8c26975d-ec0c-423c-9644-daf031c0aabb
# ╟─c707bd38-32f2-4adb-bdaf-8c8325d40ab9
# ╟─2cded3aa-8661-4477-bcc8-ac375e01622c
# ╟─989434b2-0f36-4998-b213-5ba6d2ce7166
# ╠═06ddb55b-e7a7-4ee8-91dd-171d5665df5e
# ╟─d48b5435-ce1d-4789-a33a-91f405439138
# ╟─bb7b2a43-4e93-4bfb-a41f-65ff6376babd
# ╟─5a6887a8-e651-4e4b-96c6-d1d7058cd26d
# ╟─4720b6c3-d9bf-4a46-9d38-2cee67e35142
# ╠═1186d2ae-f51e-4237-a14f-e4275b5b304c
# ╠═a32a7c6e-2824-42ee-98dc-05781f5e7ffb
# ╟─15410cfb-db47-48e0-a7e0-49a061aaaec5
# ╟─ef865aba-8f80-4f38-9560-feecdb9013a7
# ╟─9cfa6c75-cb75-470d-b037-2f01c2f22ae9
# ╠═65523ce4-ebb8-4fb2-9869-61890e0b7e2d
# ╟─9fc30c4d-e5c2-4359-9bec-e93ca5876763
# ╟─665919fd-592f-494a-9d17-adde3434a7e0
# ╟─7c40b5d3-8afd-4648-91ef-ec6235ec376e
# ╟─92791a5f-02f8-461e-9bb1-4e1077f24be9
# ╠═e24e5718-1a1a-47bc-8f25-3268ce9a5826
# ╟─13c5b2d6-7f30-4afd-8a08-3dfe4f2db184
# ╠═355bd86a-a278-4ade-96a7-2a18e0c8f7db
# ╟─297f8790-3973-462e-b14a-bff39f399e8e
# ╟─e24068de-07a2-41df-a686-01b3e21ca6b0
# ╟─4bbfd29d-832b-4e17-8f59-25463466fa37
# ╟─ee32ef8b-3b5f-4bd8-bafc-588166899782
# ╟─21f916fe-21b9-43c1-b9d1-a8e541d24e43
# ╟─0f634e30-4f59-48d4-9843-8bbb1e244db5
# ╠═b7dc76de-d6df-414e-bc76-fa9c5e2af727
# ╟─e15fdd6e-ed3e-49bf-a3a7-7530cda0d479
# ╟─5ed2c935-22d6-4e3c-ab2a-a52dd8ee462b
# ╠═ae805230-574e-4818-bdbc-1dda19fa33b0
# ╟─2e865221-74b1-4023-962c-77203686cf03
# ╟─a0af5567-02c6-4c3a-a1f2-649c622d1319
# ╟─a1cfd851-c006-4506-be18-7571927fb5ef
# ╟─f89404e4-9287-4014-bd2d-93f6833a1985
# ╟─8efa2b75-f85a-4be9-85cd-6cccd193bded
# ╟─dbf408ba-bf6e-440c-8855-3f7b7f254d6e
# ╟─74097a18-b4fe-4ec8-b47c-b279fce28598
# ╠═2eb53b12-d20d-411c-9da1-96a0b6ecac39
# ╟─9ad83c33-f843-4d63-b0f3-1a0401068c80
# ╠═f31d5edd-722e-4c55-9ef3-190b4996fd29
# ╠═7c2e9206-8403-4eab-a39e-6625945dd0c1
# ╟─c727e963-6afb-41f5-9de8-38026ab126e0
# ╟─6217aa0d-46e1-4dd6-a1e5-5eb2696757e4
# ╟─165cfb6b-63f6-4857-849b-559ba9221768
# ╟─5fc13f1f-5818-4ebd-8531-fa5900c47bc9
# ╟─1ab23f09-85e4-48fb-be05-2afe81717fbc
# ╟─1ab7ca93-8c5b-428c-9344-a12d087cc489
# ╟─ec546b10-2ab0-4d79-b81d-80532742abc2
# ╟─06f60aa7-2348-4313-a0c7-b748539659d8
# ╟─f7b8eb66-51a6-49be-b755-c6697e72bbaa
# ╟─1e4cda7d-68c7-498a-8bb8-d8642f67f344
# ╟─a1387fac-d1b6-46b3-86aa-8bd0ad868508
# ╟─a8a76616-3afb-4370-b213-2ab0f950ecf6
# ╠═26b12556-a9b8-4de7-8905-b33395ed7738
# ╟─c4153dcf-6ff6-419d-bb5b-efe91656ea28
# ╟─391ba967-0b54-4827-aac1-0382b095bfe6
# ╟─2f62830f-5b5c-4478-bc33-3b559bd5a2a3
# ╟─be19ddeb-cfbf-4904-8aa7-62669055f07c
# ╠═f48c878b-3b00-4778-8b82-e120acd5e188
# ╟─e5b9f7d1-9af4-4534-80c2-6019430da0c9
# ╟─06c05919-b335-4bda-8b14-8bfdb66d7975
# ╟─40102829-a65b-4201-847e-f2c10abb27a5
# ╟─4a1f0c7d-9fc3-4e08-97be-4eb3b8e33705
# ╟─4790f8a5-3cb3-4ccf-a8b4-85de126f8d31
# ╟─48d0de6e-ac83-4b55-94a4-1361e6f2329e
# ╟─732154d4-bf4e-40db-aec1-713208032fda
# ╟─2149e010-ccc6-46a4-a3eb-d967fb99fbf1
# ╟─254bfbc4-117b-4518-a188-c19a7c123035
# ╟─363968d9-23c9-45e7-b47d-5d7e8ad5648f
# ╟─c65c438b-d97f-4857-9897-20180cc53d25
# ╟─78e2208d-f9f8-44f7-806f-4590db85fe1a
# ╟─a7aa6fe4-3ee2-4dff-a792-4289d4fca83c
# ╟─e25e883e-e0f9-4848-adfd-24f155ecf902
# ╟─21f882ba-fa33-4695-9324-b058ba64d258
# ╟─167a446c-24b7-414c-9bdf-bae7c90df2ab
# ╟─7e94c5d7-ba40-4850-9602-18fe63381ee8
# ╟─8b8a87c1-871d-4fab-bf1a-6b3d8bc9a0fd
