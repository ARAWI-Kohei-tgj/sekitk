module sekitk.algorithm.fwelim;

import std.traits: isFloatingPoint;
import sekitk.complex.pseudo: isComplex;

template ForwardEliminationImpl(T, T Threshold)
if(isFloatingPoint!T || isComplex!T){
	import sekitk.qvm.common: TypeOfSize, TypeOfIndex, MatrixType, MajorOrder, arrayLength, arrayTo2Dim;

	template ForwardElimination(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr){
		/******************************************
		 * Function forwardElimination
		 *
		 * forward elimination
		 *
		 * Params:
		 *  mat= matrix as 2-dimensional array
		 *
		 * Returns:
		 *  struct "ResultOfFwElim" as the result
		 ******************************************/
		Result transform(ArgType)(in ArgType mat) @safe pure nothrow @nogc
	  if(is(ArgType: T[arrayLength!(Size, Size, MatrixType.dense)]) ||
			 is(ArgType: T[Size][Size])){
			import std.algorithm: swapRanges;
			import std.math: abs;

			TypeOfSize perm= 0u;
			TypeOfSize idxMaxVal;
			T[Size] temp;
			T absMaxVal, absTemp;

			static if(is(ArgType: T[arrayLength!(Size, Size, MatrixType.dense)])){
				T[Size][Size] a= arrayTo2Dim!(T, Size, Shape)(mat);
			}
			else{
				T[Size][Size] a= mat;
			}

			TypeOfSize j;
			for(j= 0u; j < Size; ++j){	// column if MatOdr is MajorOrder.row
				absMaxVal= a[j][j];
				foreach(TypeOfSize i; j..Size){	// row if MatOdr is MajorOrder.row
					absTemp= abs(a[i][j]);
					if(absTemp > absMaxVal){
						absMaxVal= absTemp;
						idxMaxVal= i;
					}
				}
				if(absMaxVal < Threshold) break;

				if(idxMaxVal > j){	// permutation
					swapRanges(a[j][], a[idxMaxVal][]);
					++perm;
				}

				foreach(i; j+1u..Size) a[i][j+1u..Size] -= (a[i][j]/a[j][j])*a[j][j+1u..Size];
			}

			return Result(a, perm, j);
		}

		/******************************************
		 * This class stores the result of the forward elimination
		 ******************************************/
		struct Result{
			import std.algorithm: merge;
			import std.range: iota, zip;
			import sekitk.indexset: IndexSetDiag, IndexSetStrictTriR;

			/************************
			 * Params:
			 *	matUpTri= result of foreard elimination
			 *	permCounter= number of the permutation in forward elimination
			 *	step= number of the elimination steps
			 ************************/
			this(in T[Size][Size] matUpTri,
					 in TypeOfSize permCounter, in TypeOfSize step) @safe pure nothrow @nogc{
				_rank= step;
				_perm= permCounter;

				T[arrayLength!(Size, Size, MatrixType.dense)] temp= (in T[Size][Size] array2D){
					T[arrayLength!(Size, Size, MatrixType.dense)] result= void;
					TypeOfIndex st, en= 0u;
					foreach(TypeOfIndex i; 0u..Size){
						st= en;
						en += Size;
						result[st..en]= array2D[i][];
					}
					return result;
/+
			  import std.array: join;
				return join(matUpTri);
+/
				}(matUpTri);

				auto idxSetSrc= (){
					auto idxSet0= IndexSetDiag!(Size, Size, Shape, MatOdr)();
					auto idxSet1= IndexSetStrictTriR!(Size, Shape, MatOdr)();
					return merge(idxSet0, idxSet1);
				}();
				auto idxSetDest= iota!TypeOfIndex(0u, arrayLength!(Size, Size, MatrixType.upperTri));
				foreach(idxDest, idxSrc; zip(idxSetDest, idxSetSrc)) _values[idxDest]= temp[idxSrc];
			}

			@property @safe pure nothrow @nogc const{
				/**********************
				 * Returns:
				 *	Raw data of forward elimination
				 **********************/
				auto arrayUpperTri(){return _values;}

				/**********************
				 * Returns:
				 *	Rank of the source matrix
				 **********************/
				uint rank() @property{return _rank;}

				///
				uint numberOfPermutation(){return _perm;}

				/**********************
				 * Returns:
				 *	determinant
				 **********************/
				T det() @property{
					import std.range: dropOne;
					T result= void;
					if(_rank < Size) result= VALUE_ZERO;	// singular matrix
					else{
						auto idxSet= IndexSetDiag!(Size, Size, MatrixType.upperTri, MatOdr)();
						result= _values[0];
						foreach(TypeOfIndex idx; idxSet.dropOne) result *= _values[idx];
						if(_perm%2 != 0) result= -result;
					}
					return result;
				}
			}

		private:
			enum T VALUE_ZERO= T(0.0L);
			TypeOfSize _rank;
			TypeOfSize _perm;
			T[arrayLength!(Size, Size, MatrixType.upperTri)] _values;
		}
	}
}
