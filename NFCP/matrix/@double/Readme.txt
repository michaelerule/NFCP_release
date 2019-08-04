
				Optimized MATLAB toolbox
		________________________________________________________

				  Elementary functions
		________________________________________________________


** Contents

	1. Introduction
	2. Requirements
	3. Installation
	4. Comments
	5. Copyright
	6. Warranty
	7. History
	8. Download
	9. Trademarks

** Publisher

	Marcel Leutenegger		marcel.leutenegger@epfl.ch
	EPFL - STI - LOB
	Bâtiment BM 5.143		Tel:	+41 21 693 77 37
	CH-1015 Lausanne



1. Introduction

	I regularily deal with a huge number of elementary functions as rounding,
	transcendentals, exponential and logarithm and was unsatisfied with the
	average performance. So, I started to provide my own assembler subroutines
	executing with full floating-point performance on any Intel Pentium II+ or
	compatible system. On recent Pentium 4 computers, I was able to improve the
	average performance by two to five times.

	This archive is dedicated to any number cruncher using MATLAB R12 build 6.0.88
	or newer.


2. Requirements

	• An Intel-based computer architecture with an Pentium II or better processor.
	• MATLAB R12 build 6.0.88 or newer running.
	• Need for speed :-)


3. Installation

	Unpack the archive in a folder that is part of the MATLAB path. There are two
	ways to improve the performance.

	      •	The subfolder '@double' contains libraries that directly overload the
		built-in functions. This has the advantage, that every existing MATLAB
		code automatically benefits from the increased performance.
		Just say "builtin(func,arguments)" to call the original version of
		'func'. Please note that this works only for built-in functions - this
		means that "which func" reports "func is a built-in function".

	      •	You may rename the functions 'func' into 'ffunc'. Existing MATLAB code
		will work with the original functions.
		In this case, you can always selectively decide whether to use the
		built-in or the external functions. Just say "ffunc" to call the fast
		version of 'func'.


4. Comments

	The libraries should always reside in a subfolder called '@double' to make sure
	they are not called for any data except of type 'double'. They do not check the
	data type of input arguments.

	Due to the overhead for calling external functions, the built-in functions work
	faster if called for matrices with less than about 4 to 8 elements. Therefore,
	if you know the size of your matrix, you can choose the appropriate function.

	The functions 'angle' and 'mod' are not built-in but MATLAB scripts.
	Therefore, the functions within this package work much faster in any case.

   Accuracy

	The results of the external functions slightly differ from the built-in functions.
	The accuracy of the external functions benefits from the floating-point registers
	with a 64bit mantissa on Intel processors. Intermediate values are kept in floa-
	ting-point registers such that rounding takes place mostly once - when writing
	the result into the output matrix with a 53bit mantissa.

	In general, a statement of "inverseFunction(function(value))" produces "value"
	with a relative error of less than 1e-13. The roundoff error of about 1.1e-16
	leads to relatively important deviations in exponentiation. Note also that in
	particular addition/subtraction in the argument of a logarithm are critical
	operations due to a relative amplification of the rounding error. The logarithm
	itself works accurate over the full complex plane R x iR. For the sake of per-
	formance, the inverse transcendental functions are currently not implemented in
	that way. See also the summary about complex functions (available online).

   Transcendental functions: pi

	An exception is made for every periodic transcendental function, where the constant
	2*pi is truncated to a 53bit mantissa. This guarantees that a statement of the
	form "sin(x)" with real 'x' always produces the expected value "sin(rem(x,2*pi))"
	at the same accuracy. The reminder is computed internaly and prescaled to a 64bit
	mantissa.

	Example 1:

		x=pow2(pi,0:256);
		plot(x,builtin('sin',x),'r',x,sin(x),'b');
		set(gca,'XScale','log');


5. Copyright

	These routines are published in terms of freeware. The author reserves the right
	to modify any of the routines listed below.

	You are allowed to distribute these routines as long as you deliver for free the
	entire package.

		Path		File		Description

		/		Linux.txt	Linux specific information
				Readme.txt	This manual
		@double/	abs.mexglx	Absolute value
				acos.mexglx	Inverse cosine
				acosh.mexglx	Inverse hyperbolic cosine
				angle.mexglx	Phase angle of complex value
				asin.mexglx	Inverse sine
				asinh.mexglx	Inverse hyperbolic sine
				atan2.mexglx	Four-quadrant inverse tangens
				atan.mexglx	Inverse tangens
				atanh.mexglx	Inverse hyperbolic tangens
				ceil.mexglx	Rounding towards +infinity
				cis.mexglx	Sine and cosine
				cis.m		Sine and cosine help
				cos.mexglx	Cosine
				cosh.mexglx	Hyperbolic cosine
				cot.m		Cotangens
				coth.m		Hypberbolic cotangens
				exp.mexglx	Exponential
				fix.mexglx	Rounding towards zero
				floor.mexglx	Rounding towards -infinity
				log2.mexglx	Binary logarithm
				log10.mexglx	Decimal logarithm
				log.mexglx	Natural logarithm
				mod.mexglx	Signed modulus
				pow2.mexglx	Scaling
				rem.mexglx	Signed remainder
				round.mexglx	Rounding to nearest
				sign.mexglx	Signum
				sin.mexglx	Sine
				sinh.mexglx	Hyperbolic sine
				sqrt.mexglx	Square root
				tan.mexglx	Tangens
				tanh.mexglx	Hyperbolic tangens


6. Warranty

	Any warranty is strictly refused and you cannot anticipate any financial or
	technical support in case of malfunction or damage.

	Feedback and comments are welcome. I will try to track reported problems and
	fix bugs.


7. History

   • January 18, 2004
	Initial release

   • February 28, 2004
	Bug fixed in rem(matrix,value): the result was stored back to the value causing
	an assertion failure.

   • June 20, 2004
	Service release thanks to a bug reported by Tom Minka in exp(matrix): the output
	was NaN for infinite input.

	This bug fix made me think about affine inputs. They are now all handled as par-
	ticular values for two reasons:
	     1.	The output is well defined. In cases with more than one possible solution,
		the function limit towards that value has been used.
	     2.	The performance does not degreade but increases considerably (table look-
		up instead of calculation). Any floating-point operation producing an
		affine result tries to throw an exception. Even if the exception is masked
		as within MATLAB, the processor calls up an internal assist slowing down
		the computation to about 10%-20% of normal performance.

8. Download

	Optimized MATLAB routines are available online at:

		   http://dmtwww.epfl.ch/~mleutene/MatLabToolbox/


	A summary is also published at MATLAB central:

			http://www.mathworks.com/matlabcentral/


9. Trademarks

	MATLAB is a registered trademark of The MathWorks, Inc. Pentium II is a
	registered trademark of Intel Corporation. Other product or brand names
	are trademarks or registered trademarks of their respective holders.

		________________________________________________________

			    Site map • EPFL © 2004, Lausanne
				Webmaster • 20/6/2004
