/*****************************************************************************
 * square & rectangular -> non-zero, binary
 * square:
 *	  Hermitian xor skew-Hermitian xor not
 *		invertible xor singular
 *		invertible:
 *					orthogonal
 *					orthonormal(unitary)
 *		Hermitian:
 *					positive definite xor semi-positive definite xor negative definite xor semi-negative definite xor indefinite
 *****************************************************************************/
module numeric.sekitk.elmcheck;

import std.bitop: bitfields;
import numeric.tribool;

/*************************************************************
 * Matrix Traits
 *************************************************************/
enum ListOfTraits: ubyte{Zero, Binary, Hermitian, SkewHermitian, Invertible, Unitary, Definite}

struct MatTraits{
	mixin(bitfields!(Tribool, is_zero, 2,
									 Tribool, is_binary, 2,
									 TriBool, is_hermitian, 2,
									 TriBool, is_skew_hermitian, 2,
									 TriBool, is_invertible, 2,
									 TriBool, is_unitary, 2,
									 TriBool, is_definite, 2));
}

/+
 +	true: 					0b11= 0x03u
 +	false: 					0b00= 0x00u
 +	indeterminate:	0b10= 0x02u
 +
 +	0b|definite|unitary|invertible|skew-hermitan|hermitian|binary|zero
 +	0b????_?***_****_****
2*4+2+1= 8+3 bytes
+/

mixin template ElementCheckCommon(CT, RT, RT threshold){
	@safe pure nothrow const{
		bool checkZeroMat(CT)(in CT[] num){
			import std.algorithm.searching: any;
			import numeric.approx: approxEqualZero;

			return num.any!( (in CT a) => !approxEqualZero!(CT, threshold)(a) );
		}
/+
		bool checkBinary(CT)(CT[] num)
		in{assert(num.length > 1u);}
		body{
			import std.algorithm.searching: any;
			CT elm0= num[0];
			num.popFront;
		  num.findAdjacent!( (in CT a) => !approxEqual);
			return num.any!();
		}+/
	}
}

mixin template ElmentCheckSq(CT, RT, RT threshold, MajorOrder matodr, ubyte size, MatrixType shape){
	@safe pure nothrow const{
		bool checkHermitian() @nogc{
			typeof(return) is_hermitian= true;

			final switch(shape){
/+
			MatrixType.Zero:
				break;
+/
			case MatrixType.Diag:
				static if(isComplex!CT){
					foreach(num; v){
						if(!approxEqualZero(num.im)){
							is_hermitian= false;
							break;
						}
					}
				}
				break;
			case MatrixType.Dense:
				break;
			case MatrixType.TriDiag:
				break;
			case MatrixType.UpperTri:
				foreach(idx; idxSetUpperTri){
					if(!approxEqualZero(v[idx])){
						is_hermitian= false;
						break;
					}
				}
				break;
			case MatrixType.LowerTri:
				foreach(idx; idxSetLowerTri){
					if(!approxEqualZero(v[idx])){
						is_hermitian= false;
						break;
					}
				}
			}
			return is_hermitian;
		}
	}

	bool checkSkewHermitian() @nogc{
		typeof(return) is_skew_hermitian= true;

		final switch(shape){
/+
		case MatrixType.Zero:
			break;
+/
		case MatrixType.Diag:
			foreach(num; v){
				if(!approxEqualZero(num)){
					is_skew_hermitian= false;
					break;
				}
			}
			break;
		case MatrixType.Dense:
		}
	}
}

/+
					bool hermitian= true;
					bool skew_hermitian= true;
					bool invertible;
					bool unitary;

					final switch(shape){
					case MatrixType.Dense:
						break;

					case MatrixType.Diag:
						flags= SqMatConstructor!(T, size, threshold).constructDiag(num);
						flags |= (CHECK_SINGULARITY | CHECK_UNITARY);
						ev.eigenvalue[]= v[];
						ev.isObtained= true;
						break;

					case MatrixType.TriDiag:
						flags= constructTriDiag(num);
						break;

					case MatrixType.UpperTri:
						// singularity
						size_t j;
						foreach(idx; idxSetDiag(shape)){
							ev.eigenvalue[j++]= v[idx];
							if(!approxEqualZero(v[idx])) skew_hermitian= false;
							else{
								invertible= false;
								static if(isComplex!T){
									if(!approxEqualZero(eval.im)) hermitian= false;
								}
							}
						}
						ev.isObtained= true;
						flags |= CHECK_SINGULARITY;

						auto upelm= idxUptri(shape);
						foreach(idx; upelm){
							if(!approxEqualZero(v[idx])){
								skew_hermitian= false;
								hermitian= false;
							}
						}
						break;

					case MatrixType.LowerTri:
						// singularity
						size_t j;
						foreach(idx; idxDiag(shape)){
							ev.eigenvalue[j++]= v[idx];
							if(!approxEqualZero(v[idx])) skew_hermitian= false;
							else{
								invertible= false;
								static if(isComplex!T){
									if(!approxEqualZero(eval.im)) hermitian= false;
								}
							}
						}
						ev.isObtained= true;
						flags |= CHECK_SINGULARITY;

						auto loelm= idxLotri(shape);
						foreach(idx; loelm){
							if(!approxEqualZero(v[idx])){
								skew_hermitian= false;
								hermitian= false;
							}
						}
						break;

					case MatrixType.Dense:
						// Hermitian or symmetric
						size_t[size] diag_idxs= idxDiag(shape);
						auto upelm= idxUptri(shape);
						auto loelm= idxLotri(shape);

						static if(isComplex!T){	// 複素數行列
							foreach(idx; diag_idxs){
								if(!approxEqualZero(v[idx])) skew_hermitian= false;
								else if(!approxEqualZero(v[idx].imag)) hermitian= false;

								if(!skew_hermitian && !hermitian) goto Sym_complex;
							}
							for(size_t i=0u ; i<upelm.length ; ++i){
								if( hermitian && (v[upelm[i]] != conj(v[loelm[i]])) ) hermitian= false;
								else if( skew_hermitian && (v[upelm[i]] != -conj(v[loelm[i]])) ) skew_hermitian= false;

								if(!skew_hermitian && !hermitian) break;
							}
						Sym_complex:
						}
						else{	// 實數行列
							foreach(idx; diag_idxs){
								if(!approxEqualZero(v[idx])){
									skew_hermitian= false;
									break;
								}
							}
							for(size_t i=0u ; i<upelm.length ; ++i){
								if(v[upelm[i]] != v[loelm[i]]) hermitian= false;
								else if(v[upelm[i]] != -v[loelm[i]]) skew_hermitian= false;

								if(!skew_hermitian && !hermitian) break;
							}
						}

						// unitary or orthonomal
						unitary= true;
						SimpleVector!row[column] colvec;
						SimpleVector!row lhs;
						T temp= T(0.0);
						foreach(size_t i; 0..column) colvec[i]= SimpleVector!row(sliceColumn(i));

					Loop: foreach(size_t pivot; 0..column){
							lhs= colvec[pivot];
							static if(isComplex!T){
								foreach(elm; lhs.a)	temp += sqAbs!double(elm);
							}
							else{foreach(elm; lhs.a) temp += elm^^2;}

							if(approxEqual(temp, 1.0, 1.0e-5, threshold)){
								unitary= false;
								break;
							}
							foreach(size_t rhs; pivot+1u..column){
								if(!approxEqualZero(lhs*colvec[rhs])){
									unitary= false;
									break Loop;
								}
							}
						}
						flags |= CHECK_UNITARY;
					}	// final switch
					flags |= (CHECK_HERMITIAN | CHECK_SKEW_HERMITIAN);
					if(hermitian) flags |= TYP_HERMITIAN;
					if(skew_hermitian) flags |= TYP_SKEW_HERMITIAN;
					if(unitary) flags |= (TYP_INVERTIBLE | TYP_UNITARY);
					else if(invertible) flags |= TYP_INVERTIBLE;
+/
