/*****************************************************************************
 * LU decomposition
 *
 * This file contains some algorithm using LU decomposition.
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 0.01
 *****************************************************************************/
module sekitk.lu;

import std.math: abs;
import std.complex;
import std.traits: isFloatingPoint;
import std.typecons;
import std.algorithm;
import std.array;
import std.range;
import mathematic.basic: Identity;
import numeric.approx: approxEqualZero;
import numeric.pseudocmplx: isComplex;

import sekitk.base;
import sekitk.indexset: IndexSetDiag, IndexSetSubDiag,
	IndexSetStrictTriR, IndexSetStrictTriL;

/*********************************************
 * Algorithm types of LU decomposition
 *********************************************/
enum LUdAlgorithm: ubyte{doolittle, crout}

/**************************************************************
 * template LUdecomposition
 **************************************************************/
template LUdecomposition(T, real Threshold,
												 TypeOfSize Size: 1u,
												 MatrixType Shape,
												 LUdAlgorithm Scheme,
												 MajorOrder MatOdr)
if(isFloatingPoint!T || isComplex!T){
	/********************************************
	 * class
	 ********************************************/
	class ArrayLUP{
		@safe pure nothrow @nogc{
			this(){}

			this(in typeof(this) other){
				this.lu= other.lu;
			}
		}

		/**************************
		 * Methods
		 **************************/
		@safe pure nothrow @nogc const{
			bool isInvertible() @property{
				return approxEqualZero!(T, Threshold)(lu)? false: true;
			}

		  // determinant
			alias det= lu;

			// inverse
			T inverseArray() @property{
				return Identity!(T, "*")/lu;
			}
		}

	private:
		T lu;
	}
}

template LUdecomposition(T, real Threshold,
												 TypeOfSize Size, MatrixType Shape,
												 LUdAlgorithm Scheme, MajorOrder MatOdr)
if((isFloatingPoint!T || isComplex!T) && Size > 1){
	/*********************************************
	 * Structure
	 *********************************************/
	enum T VALUE_ZERO= T(0.0L);
	/+++++++++++++++++++++++++++
	 + struct FwElimArray
	 +++++++++++++++++++++++++++/
	struct ResultOfFwElim{
		this(in T[Size][Size] matUpTri, in ubyte permCounter, in ubyte step) @safe pure nothrow @nogc{
			_rank= step;
			perm= permCounter;
			size_t st= 0u, en= void;
			foreach(i; 0..Size){
				en= st+Size-i;
				buf[st..en]= matUpTri[i][i..Size];
		 		st= en;
			}
		}

		@property @safe pure nothrow @nogc const{
			auto arrayUpperTri(){return buf;}

			uint rank(){return _rank;}
			uint numberOfPermutation(){return perm;}

			T det(){
				T result;
				if(_rank < Size) result= VALUE_ZERO;	// singular matrix
				else{
					TypeOfIndex j= 0u;
					result= buf[j];
					foreach(i; 1u..Size){
						j += Size-(i-1u);
						result *= buf[j];
					}
					if(perm%2 != 0) result= -result;
				}
				return result;
			}
		}

	private:
		T[arrayLength!(Size, Size, MatrixType.upperTri)] buf;
		ubyte _rank;
		ubyte perm;
	}

	/*********************************************
	 * class ArrayLUP
	 *********************************************/
  class ArrayLUP{
		/**************************
		 * Constructors
		 **************************/
		@safe pure nothrow{
			// default constructor
			this() @nogc{
			  invertible= false;
				counterEx= 0;
			}

			// copy constructor
			this(in typeof(this) other){
				this.invertible= other.invertible;
				this.lu[]= other.lu[];
			  this.p[]= other.p[];
			  this.counterEx= other.counterEx;
			}

			// normal constructor
			this(in bool isInvertible, in T[Size][Size] mat, in ubyte[Size] perm, in ubyte ex){
				invertible= isInvertible;
			  foreach(TypeOfSize i; 0u..Size)
				  foreach(TypeOfSize j; 0u..Size){{
					  lu[indexMap(i, j)]= mat[i][j];
			  }}
			  p[]= perm[];
			  counterEx= ex;
			}
		}

		/**************************
		 * Other methods
		 **************************/
		@safe pure nothrow const{
			bool isInvertible() @nogc @property{return invertible;}

			T[arrayLength!(Size, Size, MatrixType.lowerTri)] matrix(DecomposedMat Mat: DecomposedMat.lowerTri)() @nogc{
				enum TypeOfIndex LEN= arrayLength!(Size, Size, MatrixType.lowerTri);
				typeof(return) num;

				final switch(Scheme){
				case LUdAlgorithm.doolittle:
					{	// diagonal elements
					  auto idxSetDest= IndexSetDiag!(Size, Size, MatrixType.lowerTri, MatOdr)();
						foreach(TypeOfIndex idx; idxSetDest) num[idx]= Identity!(T, "*");
					}
					{	// strict lower triangular elements
						auto idxSetSrc= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
					  auto idxSetDest= IndexSetStrictTriL!(Size, MatrixType.lowerTri, MatOdr)();
						foreach(idxSrc, idxDest; zip(idxSetSrc, idxSetDest)) num[idxDest]= lu[idxSrc];
					}
					break;
				case LUdAlgorithm.crout:
					auto idxSetSrc0= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
					auto idxSetSrc1= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
					foreach(idxSrc, idxDest;
									zip(merge(idxSetSrc0, idxSetSrc1), iota(LEN))) num[idxDest]= lu[idxSrc];
				}
				return num;
			}

			T[arrayLength!(Size, Size, MatrixType.upperTri)] matrix(DecomposedMat Mat: DecomposedMat.upperTri)() @nogc{
				enum TypeOfIndex LEN= arrayLength!(Size, Size, MatrixType.upperTri);
				typeof(return) num= -3.0;

				final switch(Scheme){
				case LUdAlgorithm.doolittle:
					auto idxSetSrc0= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
					auto idxSetSrc1= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
					foreach(idxSrc, idxDest;
									zip(merge(idxSetSrc0, idxSetSrc1), iota(LEN))
									) num[idxDest]= lu[idxSrc];
					break;
				case LUdAlgorithm.crout:
					{	// diagonal elements
						auto idxSetDest= IndexSetDiag!(Size, Size, MatrixType.upperTri, MatOdr)();
						foreach(TypeOfIndex idx; idxSetDest) num[idx]= Identity!(T, "*");
					}
					{	// strict upper triangular elements
						auto idxSetSrc= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
						auto idxSetDest= IndexSetStrictTriR!(Size, MatrixType.upperTri, MatOdr)();
						foreach(idxSrc, idxDest; zip(idxSetSrc, idxSetDest)) num[idxDest]= lu[idxSrc];
					}
				}
				return num;
			}

			T[arrayLength!(Size, Size, MatrixType.dense)] matrix(DecomposedMat Mat: DecomposedMat.permutationLeft)() @nogc{
				typeof(return) num= VALUE_ZERO;
				foreach(i; 0u..Size) num[i*Size+p[i]]= Identity!(T, "*");
				return num;
			}

			// determinant
			T det() @nogc{
				typeof(return) result= lu[0];
				auto idxSet= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
				foreach(TypeOfIndex idx; idxSet.dropOne) result *= lu[idx];
				if(counterEx%2) result= -result;
				return result;
			}

			// inverse
			T[arrayLength!(Size, Size, MatrixType.dense)] inverseArray() @nogc{
				typeof(return) result= void;
				T[Size][Size] inv2d;

				foreach(TypeOfSize j; 0u..Size){
					foreach(TypeOfSize i; 0u..Size){
						inv2d[i][j]= (p[i] == j)? Identity!(T, "*"): VALUE_ZERO;
            foreach(TypeOfSize k; 0u..i) inv2d[i][j] -= lu[indexMap(i, k)]*inv2d[k][j];
					}

					foreach_reverse(i; 0u..Size){	// it must be an down counting
						foreach(k; i+1u..Size) inv2d[i][j] -= lu[indexMap(i, k)]*inv2d[k][j];
						inv2d[i][j] /= lu[indexMap(i, i)];
					}
				}

				size_t st;
				foreach(i; 0..Size){
					st= i*Size;
					result[st..st+Size]= inv2d[i];
				}
				return result;
			}
		}

	private:
		TypeOfIndex indexMap(in size_t i, in size_t j) @safe pure nothrow @nogc const{
			return cast(TypeOfIndex)(i*Size+j);
		}

		ubyte[Size] p;
		T[arrayLength!(Size, Size, MatrixType.dense)] lu;
		ubyte counterEx;
		bool invertible;
	}

	/************************************************************
	 * Functions
	 ************************************************************/
	/********************************************
	 * Function forwardElimination
	 *
	 * forward elimination
	 *
	 * Params:
	 *  mat= matrix as 2-dimensional array
	 *
	 * Returns:
	 *  struct "ResultOfFwElim" as the result
	 ********************************************/
	ResultOfFwElim forwardElimination(in T[Size][Size] mat) @safe pure nothrow @nogc{
		ubyte perm= 0u, rank;
	  TypeOfSize idxMaxVal;
		T[Size][Size] a= mat;
		T[Size] temp;
		T absMaxVal, abs_temp;

		final switch(Scheme){
		case LUdAlgorithm.doolittle:
			TypeOfSize j;
			for(j= 0u; j < Size; ++j){	// column
				absMaxVal= a[j][j];
				foreach(TypeOfSize i; j..Size){	// row
					abs_temp= abs(a[i][j]);
					if(abs_temp > absMaxVal){
						absMaxVal= abs_temp;
						idxMaxVal= i;
					}
				}
				if(absMaxVal < Threshold) break;
				else if(idxMaxVal > j){	// permutation
					swapRanges(a[j][], a[idxMaxVal][]);
					++perm;
				}
				foreach(i; j+1u..Size) a[i][j+1u..Size] -= (a[i][j]/a[j][j])*a[j][j+1u..Size];
			}
			rank= cast(typeof(rank))j;
			break;

		case LUdAlgorithm.crout:
			TypeOfSize i;
			for(i= 0u; i < Size; ++i){	// row
				absMaxVal= a[i][i];
				foreach(TypeOfSize j; i..Size){	// column
					abs_temp= abs(a[i][j]);
					if(abs_temp > absMaxVal){
						absMaxVal= abs_temp;
						idxMaxVal= j;
					}
				}
				if(absMaxVal < Threshold) break;
				else if(idxMaxVal > i){	// permutation
					swapRanges(a[i][], a[idxMaxVal][]);
					++perm;
				}
				//foreach(size_t j; i+1u..size) a[i+1u..size][j] -= (a[j][i]/a[i][i])*a[i+1u..size][i];
				/*
				foreach(size_t j; i+1u..size){
					foreach(){
						a[i+1u..size][j] -= (a[j][i]/a[i][i])*a[i+1u..size][i];
					}
					}*/
			}
			rank= cast(ubyte)i;
		}

		return ResultOfFwElim(a, perm, rank);
	}

	/*********************************************
	 * Function partialPivLU
	 *********************************************/
	ArrayLUP partialPivLU(in T[arrayLength!(Size, Size, Shape)] mat) @safe pure nothrow{
		bool invertible= true;
		T[Size][Size] a= (){
			T[Size][Size] num= void;
			foreach(i; 0u..Size) num[i]= mat[i*Size..(i+1)*Size];
			return num;
		}();
	  ubyte[Size] p= iota(Size).staticArray!Size;
		ubyte counterEx= 0;

		TypeOfSize idxMax;
		T elmMax;

		final switch(Scheme){
		case LUdAlgorithm.doolittle:
			foreach(TypeOfSize j; 0u..Size){	// column index
				elmMax= VALUE_ZERO;
				idxMax= j;
				foreach(TypeOfSize i; j..Size){	// row index
					if(a[i][j].abs > elmMax){
						elmMax= a[i][j].abs;
						idxMax= i;
					}
				}
				if(elmMax < Threshold){
					invertible= false;
					break;
				}

				if(idxMax != j){
					swap(p[j], p[idxMax]);
					swapRanges(a[j][], a[idxMax][]);
					++counterEx;
				}

				foreach(i; (j+1u)..Size){	// row index
					a[i][j] /= a[j][j];
					foreach(k; (j+1u)..Size) a[i][k] -= a[i][j]*a[j][k];
				}
			}
			break;

		case LUdAlgorithm.crout:
			foreach(TypeOfSize i; 0u..Size){	// row index
				elmMax= VALUE_ZERO;
				idxMax= i;
				foreach(TypeOfSize j; i..Size){
					if(a[i][j].abs > elmMax){
						elmMax= a[i][j].abs;
						idxMax= j;
					}
				}
				if(elmMax < Threshold){
					invertible= false;
					break;
				}

				if(idxMax != i){
					swap(p[i], p[idxMax]);
					swapRanges(a[i][], a[idxMax][]);
				}

				foreach(j; i+1u..Size){
					a[i][j] /= a[j][j];
					foreach(k; i+1u..Size) a[k][j] -= a[i][j]*a[k][j];
				}
			}
		}

		return new ArrayLUP(invertible, a, p, counterEx);
	}

	/// ditto
	ArrayLUP!Shape partialPivLU(MatrixType Shape: MatrixType.band3
															)(in T[arrayLength!(Size, Size, Shape)]){
		return new ArrayLUP(invertible, a, p, counterEx);
	}
}	// template `LUdecomposition'
