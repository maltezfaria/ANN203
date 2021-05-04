### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d54d301e-060e-4c14-b5c8-1fca26dd84f9
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="LinearAlgebra"), 
			Pkg.PackageSpec(name="FileIO"), 
			Pkg.PackageSpec(name="PlutoUI"), 
			Pkg.PackageSpec(name="Images"),
			Pkg.PackageSpec(name="FixedPointNumbers")
			])
	using LinearAlgebra
	using FileIO
	using PlutoUI
	using Images
	using FixedPointNumbers
end

# ╔═╡ 606540b0-9846-11eb-2a09-8112c5081854
md"""
# SIM 203
"""

# ╔═╡ 57160c99-cccc-4aca-99fd-df9356da76a1
PlutoUI.TableOfContents(;depth=2)

# ╔═╡ 28280ba3-789f-40ec-b731-cbc43334b839
md"""!!! danger "Submission instructions"
	**Due date:** May 6, 2021
	
	Send your completed notebook to **luiz.maltez-faria@inria.fr**
"""

# ╔═╡ 6ea0197b-5af3-4116-b214-a27f28508c33
md"""## Introduction

In this second homework we will explore the connection between images and linear algebra. In particular, we will see how images can be viewed as matrices of pixels, where the actual representation of pixel will depend, among other things, on whether we are dealing with colored or grayscale images. 

As before, you are encouraged to consult the course material available [here](https://perso.ensta-paris.fr/~mbonnet/aln.pdf) for further discussion on the theoretical aspects covered in this homework. *Please comment your code* so that I can understand the reasoning behind what you wrote.

Unlike the first assignement, there are no formal tests that your code has to pass. Since we will be dealign with images, most of things we will do will be judged by the *eye norm*; that is, it should look correct. Of course I will inspect your code to make what you wrote is in fact sensible. 

Remember that if you running this notebook using `Pluto`, everything you see is interactive so you can play with the images as you go along!

"""

# ╔═╡ 8c26975d-ec0c-423c-9644-daf031c0aabb
md""" ## Images as matrices

You have probably heard this one before, but here it goes again: images are *just* matrices. What I mean by that is that you can think of an image as a `Matrix` of pixels. For grayscale images, a pixel is typically represented by a single byte (8 bits) enconding the *grayness*. Colored images on the other hand typically store in each pixel three bytes encoding the amount of red, green, and blue (`RGB`). To make our life easier (and because this is not a course on image processing) we will rely on the types `Gray` and `RGB` defined on the package `Images` to represent such pixels. Since `Pluto` knows how to display those, rendering a Matrix of `RGB` or `Gray` pixels should be effortless. Here is an example of how to create a vector of `Gray` objects:
"""

# ╔═╡ 532fc383-556d-4a50-a198-6d28d15c28a7
[Gray(i/255) for i in 0:1:255]

# ╔═╡ e1f3a12c-8667-4f38-84af-10ba265bb302
md"""Notice that that this exactly how you would create a vector of e.g. integers. Similarly, you can see first hand where `RGB` gets its name from:"""

# ╔═╡ 4b58ae89-2581-46a5-8a08-8e82c96ef441
[RGB(1,0,0) RGB(0,1,0) RGB(0,0,1)]

# ╔═╡ 31e0090a-a286-4221-96dc-32e6e7e5afa5
md"""
As you may have noticed a `Gray` pixel takes values ranging from $0$ (black) to $1$ (white). The underlying representation is variable, but for most images a single grayscale pixel takes a byte of memory (like a `char` in C). Likewise, the `RGB` type uses three numbers between $0$ and $1$ to determine the "amount" of reed/green/blue in the pixel, and each of these is usually represented using a byte so that an `RGB` pixel occupies 24 bits. 

With this in mind, we can generate some basic images now by manipulating the underying `Matrix`:
"""

# ╔═╡ fddfc105-2ffc-4a48-b483-edb6fc5b60f9
function french_flag(h,w)
	imag = zeros(RGB{N0f8},h,w) # initalize a zero "matrix" of appropriate type
	for i in 1:h
		for j in 1:w
			# fill in pixels
		end
	end
	return imag
end

# ╔═╡ bbb1c6fe-1858-426d-b298-ec11dbcd7981
imag = french_flag(200,400)

# ╔═╡ 9fc1aacf-91ef-4070-b5eb-6721d3e6763a
md"""
We will later work with other images, both created artificially and downloaded from the internet. For now, let us focus on the `imag` of the French flag you just created and see what type of operations we can do with it. 
"""

# ╔═╡ d0c2213a-050d-4c50-850b-c72072ad4341
md"""## Convolutions and filters
Now that we know that images are (roughly) just a matrix of pixels, and that a pixel is just a number, we can start modifying images by modifying the corresponding matrix entries. Over the next few questions, we will implement various functions that iterate over an image modifying the pixels to achieve a desired functionality. 

### Conversion to grayscale

The simplest conversion is probably to transform an `RGB` image to grayscale using an  average of the colors. For reasons related to how light is perceived, a weighted average between the red, green, and red intensities is more appropriate. Here is how the gray color will be calculated:
```math
c = 0.298r + 0.588g + 0.114b
```
where $r,g,b$ are the red, green, and blue values. 
"""

# ╔═╡ 5d05d2c2-5f05-4c66-bfb8-096615b690cd
function convert_to_grayscale_pix(pix)
	# the weird letters at the end of the constants specify the number type
	c = (0.298N0f8*pix.r)+(0.588N0f8*pix.g)+(0.114N0f8*pix.b)
	Gray(c)
end

# ╔═╡ 5fab3cba-0a5a-4bf6-8913-d0ad4bc5d28b
# how to use it
convert_to_grayscale_pix(RGB(1,0,1)) 

# ╔═╡ 39353f86-49e3-49d1-8a20-e0daad972470
function convert_to_grayscale(imag)
	out = zeros(Gray{N0f8},size(imag))
	# fill in the pixels
	return out
end

# ╔═╡ 6b74e301-fd5b-4cc3-9b59-fd2db1f14667
imag_gray = convert_to_grayscale(imag)  

# ╔═╡ c7ba64fa-0ab8-4d1e-9e30-35e6b2831d07
md"""### Blurring 
Next we will implement a *blurring* filter where we will replace each pixel in the image by the average of its $9$ neighbors (itself included). That is, `imag[i,j]` is replaced by the average of the pixels in `imag[i-1:i+1,j-1:j+1]`. The border cases require some special attention. 
"""

# ╔═╡ b86fc0bd-ee79-4020-9e50-4b849b14d665
"""
	blur(imag)

Blur a grayscale image by replacing each pixel with its average over the `3×3` block surrounding it. 
"""
function blur(imag::Matrix)
	# make sure this is a grayscale image (i.e. elements are `Gray` pixels)
	@assert eltype(imag) <: Gray
	# comment the rest of this function!
	out = zero(imag)
	m,n = size(imag)
	for i in 1:m
		for j in 1:n
			acc = 0.0
			for k in -1:1
				for l in -1:1
					row = i-k
					col = j-l
					row = min(max(row,1),m)
					col = min(max(col,1),n)
					acc += imag[row,col]/9
				end
			end
			out[i,j] = acc
		end
	end
	return out
end

# ╔═╡ e73031d2-907d-4ab1-9746-26725e9c0312
blur(imag_gray)

# ╔═╡ 9d74220f-2265-439f-a96c-b933e1bb5b2c
md"### Convolutions
Many interesting filters are based on the idea of *convolving* an image with a given 
kernel. Recall the definition of a two-dimensional (continuous) convolution: 
```math
(g \star f)(x_1,x_2) = \int_{-\infty}^\infty\int_{-\infty}^\infty k(y_1,y_2) f(x_1-y_1,x_2-y_2) dy_1 dy_2,
```
where $g : \mathbb{R}^2 \to \mathbb{R}$ and $F : \mathbb{R}^2 \to \mathbb{R}$. 

Now the discrete version of this can be thought of as 

```math
(G \star F)[i,j] = \sum_{k=-\infty}^\infty\sum_{l=-\infty}^\infty G[k,l] F[i-k,j-l] 
```
where $G : \mathbb{Z}^2 \to \mathbb{R}$ and $F : \mathbb{Z}^2 \to \mathbb{R}$.

Look at the core of the `blur` function and compare it with the formulae above. Do you see how to relate discrete convolutions to the `blur` filter? 
"

# ╔═╡ 16c0bd4d-4734-4894-b18b-4999b5d80001
md"""
The discussion above suggests writting a somewhat *generic* function to apply a kernel to an image. This is your next task:
"""

# ╔═╡ 22789eba-d3be-4287-89f5-9c1681741946
"""
	convolve(G,imag,w)

Compute `G ⋆ imag`, i.e. the convolution of `imag` with the kernel `G`. The argument `w` is the window size, and `size(G) == (2w+1)×(2w+1)`. The pixel `imag[i,j]` is the discrete convolution of `G` with the `imag[i-w:i+w,j-w:j+w]` pixels. 

For example, the `blur` function implemented above can be viewed as convolving with 
`
G = [	1/9 1/9 1/9;
		1/9 1/9 1/9;
		1/9 1/9 1/9
	]	
`
so that `w=1` in that case.
"""
function convolve(G::Matrix,imag::Matrix,w::Int)
	@assert eltype(imag) <: Gray
	out = zero(imag)
	m,n = size(imag)
	# fill in here
	return out
end

# ╔═╡ 91449998-8036-4ffa-8193-9421ccff4720
md"If your convolve function is both *generic* and *correct*, moving the slider below should increase the blur of the `imag_gray` (i.e. the French flag):" 

# ╔═╡ 5d5c9419-2e01-4968-b3be-f2ed53cbfb0c
@bind w Slider(1:50)

# ╔═╡ ad06fea7-3b95-41df-b019-79380085ddd8
md"smoothing width --> $w"

# ╔═╡ 5b050842-779e-4b8d-b6b1-5909b87afa2c
let
	sz = 2w+1
	G = 1/(sz^2)*ones(sz,sz)
	convolve(G,imag_gray,w)
end

# ╔═╡ 9de16c5e-bb86-4a31-8013-815a85bdee13
md""" ### Gaussian blur
The nice thing about the generic version using `convolve` is that it now becomes rather easy to apply other *filters*. For instance, a more common type of *blurring* is the [Gaussian blur](https://en.wikipedia.org/wiki/Gaussian_blur). It consists of convoling your image with a kernel which approximates
```math
G(x,y) = \frac{1}{2\pi \sigma^2} e^{-\frac{x^2 + y^2}{2\sigma^2}}
```
Since this function integrates to $1$ over $\mathbb{R}^2$, it acts as a weighted average of the pixels in our image, giving more weight to pixels which are close to the target pixel. In many applications, this is a better filter than the simpler average implemented in `blur`. 

In practice a truncated version of this kernels is employed over a window of width `w`, as done in the following function:
"""

# ╔═╡ c1fd08c3-6325-4459-a0cf-064bc966b300
"""
	gaussian_kernel(σ,w)

A matrix of size `(2w+1)×(2w+1)` implementing a truncated Gaussian kernel with starndard deviation σ. 
"""
function gaussian_kernel(σ,w) 
	out = [exp(-(i^2 + j^2)/(2σ^2)) for i in -w:w, j in -w:w]
	out /= sum(out) # rescale
end

# ╔═╡ 7b8970bd-29ee-4a00-8e6a-ce8382b3186e
# how to use it
gaussian_kernel(1,2)

# ╔═╡ 32190a4b-8b41-465b-89a3-0110942f94c1
md"""The averaging nature of the [Gaussian blur](https://en.wikipedia.org/wiki/Gaussian_blur) filter has the effect of smoothening features of an image, and this can be used to reduce image noise in certain circumstances. The next question is an exploratory one: you are to play with the parameters of the Gaussian blur filter to try and produce a clearer image of the following noised flag:"""

# ╔═╡ 904faf2e-7bc3-4f66-afa4-2cd607bceec9
flag_gray_noise = imag_gray + 0.1*randn(size(imag_gray)...)

# ╔═╡ 44429ea9-86df-4e01-93a1-69677dc6bf69
@bind wi Slider(0:20)

# ╔═╡ 18e7cdc3-aed8-4aa7-9ab4-5d5b9bdff19b
@bind σ Slider(0.1:0.1:10)

# ╔═╡ 95a2a42d-8aba-438b-8b5f-ef1140b915ce
md"""
**Guassian blur parameters**:

| w   | σ |           
| --- |---| 
| $wi  | $σ| 
"""

# ╔═╡ b25bac36-df56-4f4a-8a94-f8603683404b
let
	convolve(gaussian_kernel(σ,wi),flag_gray_noise,wi)
end

# ╔═╡ 18b1d7da-7503-4eec-9232-cdec03dc619a
md"""## Compressing an image

In this last part we will look at how to compress images. The resolution of an image is the number of pixels in it; that is, the `size` of the underlying matrix. How much resolution is needed to properly represent an image depends both on the characteristics of the image, as well as on the size that you want to display it on. In what follows we will look at some basic techniques to compress an image. 

To change things a bit, we will work with an image downloaded from the internet:
"""

# ╔═╡ c5c584d2-5e31-4fc4-bd59-6ee14e60225e
function get_image(url)
	fname = download(url)
	load(fname)
end

# ╔═╡ 4653ba8a-2ece-4e35-b26f-1ccc6fa35a23
bedroom = get_image("https://www.artic.edu/iiif/2/25c31d8d-21a4-9ea1-1d73-6a2eca4dda7e/full/843,/0/default.jpg")

# ╔═╡ b294dc8a-95c8-4448-830e-c1837f3ea7ac
md"""### Downsampling
The most basic technique is downsampling. That is, extracting every $k$ pixel from the image, where $k \in \mathbb{N}$ is a downsampling factor. For the sake of completeness, implement it below.
"""

# ╔═╡ e4393765-62f3-47f7-b384-219bd65338c9
size(bedroom)

# ╔═╡ 9eda5160-1fce-4ba7-86ac-6834ddb0ae8c
function down_sample(im,k)
	nothing
end

# ╔═╡ 41ee8b95-9687-4467-b75c-8a1f1289d988
bedroom_downsampled = down_sample(bedroom,4)

# ╔═╡ 40d177a6-36e8-47bc-b4d4-af7e21518c25
md"### SVD compresion"

# ╔═╡ 4fed707c-129d-43c9-854e-ae465c72c3ef
md"""
A more interesting compression technique is based on computing the SVD of the underlying matrix, and keeping only the first $k$ singular values and singular vectors. This yields a compressed representation of the matrix $M \in \mathbb{R}^{m\times n}$ provided $k \ll min(m,n)$. More precisely, if 
```math
M = U \Sigma V^*
```
then we will denote by *truncated SVD* the matrix
```math
M_k = U_{:,1:k} \Sigma_{1:k,1:k} V^*_{1:k,:}
```
Over the next few questions you will analyze this representation and implement it for image compresion. 
"""

# ╔═╡ 6a656d99-0cde-47d1-b4ab-5716d7d53da4
md"""
We will now implement a compression technique based on the truncated SVD. To make things simpler, I have provided a `struct` to represent the compressed image, as well as the code to display it:
"""

# ╔═╡ 19ab5bcc-b09b-4464-9d2c-2e27a722c4b3
struct CompressedGrayImage
	U::Matrix{Float64}
	Σ::Vector{Float64}
	V::Matrix{Float64}
end

# ╔═╡ 3ebb52fa-ba57-4f51-aaa0-1726c411b0b0
# this will assemble the Mₖ matrix and prepare it for display
decompress(im::CompressedGrayImage) = Gray.(im.U*Diagonal(im.Σ)*adjoint(im.V))

# ╔═╡ 90643cc3-f8aa-45a1-99bf-6e954ce42877
function compress_gray_image(img,k)
end

# ╔═╡ bd28703b-e871-4723-a85f-1fed859187d7
function circle_image(h,w)
	im = Matrix{Gray{N0f8}}(undef,h,w)
	ic,jc = h/2, w/2
	r = min(h,w)/4
	for i in 1:h
		for j in 1:w
			if ((i-ic)^2 + (j-jc)^2 < r^2)
				im[i,j] = Gray(0)
			else
				im[i,j] = Gray(1)
			end
		end
	end
	return im
end

# ╔═╡ dabb1c2d-ab10-4d33-b0be-52fd2df3ac70
circ_imag = circle_image(200,200) 

# ╔═╡ db7968c4-8a32-4bb0-9499-69adba3c940a
md"## Playground
*A place to play with what you coded!*
"

# ╔═╡ 254bfbc4-117b-4518-a188-c19a7c123035
md"""## Utility functions
"""

# ╔═╡ 363968d9-23c9-45e7-b47d-5d7e8ad5648f
question(text,number="") = Markdown.MD(Markdown.Admonition("tip", "Question $number", [text]))

# ╔═╡ 18903e3c-d037-4272-b225-715d2027621d
question(md"""
	Now that you understand how images are just matrices of pixels, complete the function below to generate an image of the French flag.""",1)

# ╔═╡ 315dc0ea-cc1e-4608-bec7-6983b7794f8e
question(md"Use the provided `convert_to_grayscale_pix(pix::RGB)` method to complete the function `convert_to_grayscale` below.",2)

# ╔═╡ 9d25d28f-45d3-415a-9b20-e437d779cfa5
question(md"Comment the `blur` function below to explain how it works. Make sure to explain where/how the border cases are handled. Feel free to suggest improvements or modify it if you can justify your changes.",3)

# ╔═╡ fdfdded0-3e5d-4e40-a689-c78760ad8f36
question(md"Based on the previous discussion on discrete convolutions, re-express the `blur` function as a convolution of your image with a given kernel $G$. What is the complexity of `blur` (number of flops) in terms of size of the image?",4)

# ╔═╡ 72c7379c-5397-4f6d-99a6-3b7812d19c65
question(md"Complete the function below to convolve your image with a generic kernel $G$. A good staring point is the `blur` function provided.",5)

# ╔═╡ 2bdc912e-d5a9-433f-b26d-ca922f69424c
question(md"Play with the window and standard deviation of the Gaussian blur to reduce the noise in the image above. Can you explain why such a process can reduce the perceived noise?",6)

# ╔═╡ 4de21fb4-3c01-4f0d-875c-9b0e35c11f37
question(md"Complete the function below to downsample an image. Apply it to the `bedroom` high resolution image below with a factor of $k=4$. How much smaller is the downsampled image?",7)

# ╔═╡ ff260677-06bf-49d4-a95e-ceb667869a7c
question(md"What is the error $||M - M_k||_2$ of the approximation?",8)

# ╔═╡ 39dee65e-8fac-4000-a4aa-c515098707cb
question(md"""Estimate the storage cost of the truncated SVD representation. What is the critical value of $k$ beyond which storing the truncated SVD becomes more expensive than storing the full matrix. You may assume that the representation of the individual entries of both the SVD and the full matrix require the same amount of bits.""",9)

# ╔═╡ 104fcfd8-b43a-4c56-a5cc-6394a2b8fec7
question(md"Complete the code below perform a grayscale image `img`. You may use the `svd` function from `LinearAlgebra`. How many singular values do you need to represent the grayscale version of the French flag? Hint: look at the singular values.",10)

# ╔═╡ 40748a9e-5af7-48ea-8988-3362dce5b6b5
question(md"For the `circ_imag` below, how large must `k` be so as to obtain a visually acceptable approximation of the circle? Why is this larger than the French flag example?",11)

# ╔═╡ 9a16103d-643b-42f7-9447-d68c42a5373a
question(md"""Finally, apply your compression technique to a grayscale version of the `bedroom` image. How large does `k` have to be obtain a good approximation under the *eye norm*?
""",12)

# ╔═╡ c65c438b-d97f-4857-9897-20180cc53d25
answer(text=md"Write your answer here!
") = Markdown.MD(Markdown.Admonition("warning", "Answer", [text]))

# ╔═╡ 5e5fe53d-af80-48b1-bcc2-b74cdb8662ec
answer(md"""
Write here.
""")

# ╔═╡ 1dc98212-7494-4617-99ec-0be0c172864c
answer(md"""
Answer here.
""")

# ╔═╡ 93b6d6d7-7b76-4df2-8fdb-0e6dc969fa0d
answer()

# ╔═╡ b0a2201f-7afd-41b4-8f36-e334283960d0
answer()

# ╔═╡ bcfbdb37-c763-45ca-9ac6-c5b4872e9536
answer()

# ╔═╡ d5ebb24d-6d40-4669-b423-adadbb6530c0
answer()

# ╔═╡ 78e2208d-f9f8-44f7-806f-4590db85fe1a
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]))

# ╔═╡ 4f880d3c-11e9-4968-9e9d-7cbed480e1f6
hint(md"""
Let $F$ be an $m \times n$ matrix with indices $1 \leq i \leq m$, $1 \leq j \leq n$. Consider the following operation:
	```math
(G \star F)[i,j] = \sum_{k=-1}^1\sum_{l=-1}^1 G[k,l] F[i-k,j-l], 
```
where $\quad i = 2, \ldots, m-1$, $j = 2,\ldots,n-1$.
Find the appropriate $3\times 3$ matrix $G$ so that this operations reduces to the `blur` function. Don't forget to mentiond how the border cases $i=1,m$ and $j=1,n$ should be handled.
""")

# ╔═╡ e25e883e-e0f9-4848-adfd-24f155ecf902
correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("note", "Got it!", [text]))

# ╔═╡ 21f882ba-fa33-4695-9324-b058ba64d258
note(text) = Markdown.MD(Markdown.Admonition("note", "Note", [text]))

# ╔═╡ 8c0961b4-4982-4355-9f7a-6aa37e342f3c
note(md"The cell below will add the dependencies that you need for this notebook. This may take some time the first time you run it since some packages will be downloade and precompiled.")

# ╔═╡ ecc3904b-a8f9-4b5f-a7ff-a2ecca9ba70d
note(md"""In reality things are a bit more complicated. The `Gray` type from `Colors` is actually a parametric type, which is similar to a templated `C++` struct. The templated parameter is the underlying representation of the number. So for instance `Gray(0f0)` is of type `Gray{Float32}` while `Gray(0.0)` is of type `Gray{Float64}`. The "normal" representation uses the `N0f8` type and occupies $1$ byte -- this is the reason why you many see `N0f8` floating around. You can ignore these details.""")

# ╔═╡ 99142b26-928a-4e4a-acee-9370f478738d
note(md"""
From here one, almost everything we will do will be on grayscale images. This is for sake of simplicity, and most operations/filters in the next few sections can be applied to `RGB` images by independently appling the operation/filter on each of color channels.  	
""")

# ╔═╡ c7cadb29-addf-4c6a-9587-d7b4f30c3291
note(md"""
There are many other interesting filters that you can apply to images to e.g. sharpen edges, change colors, increase brightness, etc. See for example the [Sobel operator](https://en.wikipedia.org/wiki/Sobel_operator). Feel free to implement and play around with your custom filters at the end of this notebook.
""") 

# ╔═╡ f13a3759-8c03-46f8-b39c-ad459c281e15
note(md"""Assuming the kernel $G$ is known, an interesting question (with practical consequences)  is whether or not one may recover the original image `imag` given observation of `G ⋆ imag`. This is an an inverse problem, and linear algebra can be quite useful to help answer such questions. 
""")

# ╔═╡ 167a446c-24b7-414c-9bdf-bae7c90df2ab
warning(text) = Markdown.MD(Markdown.Admonition("warning", "Warning", [text]))

# ╔═╡ Cell order:
# ╟─606540b0-9846-11eb-2a09-8112c5081854
# ╟─57160c99-cccc-4aca-99fd-df9356da76a1
# ╟─8c0961b4-4982-4355-9f7a-6aa37e342f3c
# ╠═d54d301e-060e-4c14-b5c8-1fca26dd84f9
# ╟─28280ba3-789f-40ec-b731-cbc43334b839
# ╟─6ea0197b-5af3-4116-b214-a27f28508c33
# ╟─8c26975d-ec0c-423c-9644-daf031c0aabb
# ╠═532fc383-556d-4a50-a198-6d28d15c28a7
# ╟─e1f3a12c-8667-4f38-84af-10ba265bb302
# ╠═4b58ae89-2581-46a5-8a08-8e82c96ef441
# ╟─31e0090a-a286-4221-96dc-32e6e7e5afa5
# ╟─ecc3904b-a8f9-4b5f-a7ff-a2ecca9ba70d
# ╟─18903e3c-d037-4272-b225-715d2027621d
# ╠═fddfc105-2ffc-4a48-b483-edb6fc5b60f9
# ╠═bbb1c6fe-1858-426d-b298-ec11dbcd7981
# ╟─9fc1aacf-91ef-4070-b5eb-6721d3e6763a
# ╟─d0c2213a-050d-4c50-850b-c72072ad4341
# ╠═315dc0ea-cc1e-4608-bec7-6983b7794f8e
# ╠═5d05d2c2-5f05-4c66-bfb8-096615b690cd
# ╠═5fab3cba-0a5a-4bf6-8913-d0ad4bc5d28b
# ╠═39353f86-49e3-49d1-8a20-e0daad972470
# ╠═6b74e301-fd5b-4cc3-9b59-fd2db1f14667
# ╟─99142b26-928a-4e4a-acee-9370f478738d
# ╟─c7ba64fa-0ab8-4d1e-9e30-35e6b2831d07
# ╠═9d25d28f-45d3-415a-9b20-e437d779cfa5
# ╠═b86fc0bd-ee79-4020-9e50-4b849b14d665
# ╠═e73031d2-907d-4ab1-9746-26725e9c0312
# ╟─9d74220f-2265-439f-a96c-b933e1bb5b2c
# ╟─fdfdded0-3e5d-4e40-a689-c78760ad8f36
# ╟─4f880d3c-11e9-4968-9e9d-7cbed480e1f6
# ╠═5e5fe53d-af80-48b1-bcc2-b74cdb8662ec
# ╟─16c0bd4d-4734-4894-b18b-4999b5d80001
# ╟─72c7379c-5397-4f6d-99a6-3b7812d19c65
# ╠═22789eba-d3be-4287-89f5-9c1681741946
# ╟─91449998-8036-4ffa-8193-9421ccff4720
# ╟─5d5c9419-2e01-4968-b3be-f2ed53cbfb0c
# ╟─ad06fea7-3b95-41df-b019-79380085ddd8
# ╠═5b050842-779e-4b8d-b6b1-5909b87afa2c
# ╟─9de16c5e-bb86-4a31-8013-815a85bdee13
# ╠═c1fd08c3-6325-4459-a0cf-064bc966b300
# ╠═7b8970bd-29ee-4a00-8e6a-ce8382b3186e
# ╟─c7cadb29-addf-4c6a-9587-d7b4f30c3291
# ╟─32190a4b-8b41-465b-89a3-0110942f94c1
# ╟─904faf2e-7bc3-4f66-afa4-2cd607bceec9
# ╟─2bdc912e-d5a9-433f-b26d-ca922f69424c
# ╟─44429ea9-86df-4e01-93a1-69677dc6bf69
# ╟─18e7cdc3-aed8-4aa7-9ab4-5d5b9bdff19b
# ╟─95a2a42d-8aba-438b-8b5f-ef1140b915ce
# ╠═b25bac36-df56-4f4a-8a94-f8603683404b
# ╠═1dc98212-7494-4617-99ec-0be0c172864c
# ╟─f13a3759-8c03-46f8-b39c-ad459c281e15
# ╟─18b1d7da-7503-4eec-9232-cdec03dc619a
# ╠═c5c584d2-5e31-4fc4-bd59-6ee14e60225e
# ╠═4653ba8a-2ece-4e35-b26f-1ccc6fa35a23
# ╟─b294dc8a-95c8-4448-830e-c1837f3ea7ac
# ╟─4de21fb4-3c01-4f0d-875c-9b0e35c11f37
# ╠═e4393765-62f3-47f7-b384-219bd65338c9
# ╠═9eda5160-1fce-4ba7-86ac-6834ddb0ae8c
# ╠═41ee8b95-9687-4467-b75c-8a1f1289d988
# ╟─93b6d6d7-7b76-4df2-8fdb-0e6dc969fa0d
# ╟─40d177a6-36e8-47bc-b4d4-af7e21518c25
# ╟─4fed707c-129d-43c9-854e-ae465c72c3ef
# ╟─ff260677-06bf-49d4-a95e-ceb667869a7c
# ╠═39dee65e-8fac-4000-a4aa-c515098707cb
# ╟─6a656d99-0cde-47d1-b4ab-5716d7d53da4
# ╠═19ab5bcc-b09b-4464-9d2c-2e27a722c4b3
# ╠═3ebb52fa-ba57-4f51-aaa0-1726c411b0b0
# ╟─104fcfd8-b43a-4c56-a5cc-6394a2b8fec7
# ╠═90643cc3-f8aa-45a1-99bf-6e954ce42877
# ╠═b0a2201f-7afd-41b4-8f36-e334283960d0
# ╠═40748a9e-5af7-48ea-8988-3362dce5b6b5
# ╠═bcfbdb37-c763-45ca-9ac6-c5b4872e9536
# ╠═bd28703b-e871-4723-a85f-1fed859187d7
# ╠═dabb1c2d-ab10-4d33-b0be-52fd2df3ac70
# ╠═9a16103d-643b-42f7-9447-d68c42a5373a
# ╠═d5ebb24d-6d40-4669-b423-adadbb6530c0
# ╟─db7968c4-8a32-4bb0-9499-69adba3c940a
# ╟─254bfbc4-117b-4518-a188-c19a7c123035
# ╟─363968d9-23c9-45e7-b47d-5d7e8ad5648f
# ╟─c65c438b-d97f-4857-9897-20180cc53d25
# ╟─78e2208d-f9f8-44f7-806f-4590db85fe1a
# ╟─e25e883e-e0f9-4848-adfd-24f155ecf902
# ╟─21f882ba-fa33-4695-9324-b058ba64d258
# ╟─167a446c-24b7-414c-9bdf-bae7c90df2ab
