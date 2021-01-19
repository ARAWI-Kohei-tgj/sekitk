module sekitk.qvm.algorithm.decomposition; 

import std.traits: isFloatingPoint;
import sekitk.complex.pseudo: isComplex;

template MatrixDecomp(T, real Threshold)
if(isFloatingPoint!T || isComplex!T){
	import std.meta: templateOr;
	import sekitk.qvm.common: TypeOfSize, MatrixType, MajorOrder, DecompScheme;

	alias isNormalLU= templateOr!(isDoolittle, isCrout);

	enum T VALUE_ZERO= T(0.0L);//Identity!(T, "+");
	enum T VALUE_ONE= T(1.0L);//Identity!(T, "*");

	private{
		template isDoolittle(DecompScheme Method){
			static if(Method is DecompScheme.luWithNonPivDoolittle ||
								Method is DecompScheme.luWithPartialPivDoolittle ||
								Method is DecompScheme.luWithFullPivDoolittle){
				enum bool isDoolittle= true;
			}
			else{
				enum bool isDoolittle= false;
			}
		}

		template isCrout(DecompScheme Method){
			static if(Method is DecompScheme.luWithNonPivCrout ||
								Method is DecompScheme.luWithPartialPivCrout ||
								Method is DecompScheme.luWithFullPivCrout){
				enum bool isCrout= true;
			}
			else{
				enum bool isCrout= false;
			}
		}

		template isNonPivotting(DecompScheme Method){
			static if(Method is DecompScheme.luWithNonPivDoolittle ||
								Method is DecompScheme.luWithNonPivCrout){
				enum bool isNonPivotting= true;
			}
			else{
				enum bool isNonPivotting= false;
			}
		}

		template isPartialPivotting(DecompScheme Method){
			static if(Method is DecompScheme.luWithPartialPivDoolittle ||
								Method is DecompScheme.luWithPartialPivCrout){
				enum bool isPartialPivotting= true;
			}
			else{
				enum bool isPartialPivotting= false;
			}
		}

		// transversal for 2D array
		T[Major] transversal(TypeOfSize Major, TypeOfSize Sub)(T[Sub][Major] a, TypeOfSize idxSub){
			typeof(return) result= void;
			foreach(TypeOfSize idxMajor; 0u..Major) result[idxMajor]= a[idxMajor][idxSub];
			return result;
		}
	}

	template DecompLU(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr){
		import sekitk.qvm.common;
		import sekitk.qvm.indexset: IndexSetDiag, IndexSetStrictTriL, IndexSetStrictTriR;
		static if(Shape is MatrixType.dense){	// dense
			import sekitk.qvm.algorithm.decomposition.lu_dense: LU_DENSE_IMPL1, LU_DENSE_IMPL2;
			mixin((Size > 1)? LU_DENSE_IMPL2: LU_DENSE_IMPL1);
		}
		else static if(Shape is MatrixType.band1){
			import sekitk.qvm.algorithm.decomposition.lu_dummy: LU_BAND1_IMPL;
			mixin(LU_BAND1_IMPL);
		}
		else static if(Shape is MatrixType.band3){
			import sekitk.qvm.algorithm.decomposition.lu_band3: LU_BAND3_IMPL;
			mixin(LU_BAND3_IMPL);
		}

	private:
		enum string PARTIAL_PIV_IMPL= q{
			idxMax= idxD +((idxD < Size-1)?
										 maxIndex!((a, b) => abs(a) < abs(b))(a.transversal(idxD)[idxD..Size])
										 : 0u);
			if(abs(a[idxMax][idxD]) < Threshold){
				invertible= false;
				break;
			}
			else if(idxMax != idxD){
				swap(p[idxD], p[idxMax]);
				swapRanges(a[idxD][], a[idxMax][]);
				++counterEx;
			}
		};

		enum string FULL_PIV_IMPL= q{
			idxMax= maxIndex2D(idxSub, Size, idxSub, Size);
			if(idxMax[0] != idxSub){
				swap(p[idxSub], p[idxMax[0]]);
				swapRanges(a[idxSub][], a[idxMax[0]][]);
				++counterEx[0];
			}

			if(idxMax[1] != idxSub){
				swap(p[idxSub], p[idxMax[1]]);
				foreach(){}
				++counterEx[1];
			}
		};

		/**************************
		 * one dimensional array -> two dimensional array
		 **************************/
		T[Sub][Major] arrayTo2Dim(T, TypeOfSize Major, TypeOfSize Sub)(in T[Major*Sub] num) @safe pure nothrow @nogc{
			import std.array: staticArray;
			import std.range: chunks;
			typeof(return) result= void;

			auto temp= num[].chunks(Sub);
			foreach(i; 0u..Major) result[i]= temp[i].staticArray!Sub;

			return result;
		}

		/**
		 *
		 **/
		size_t[2] maxIndex2D(T[Size][Size] a,
														 TypeOfSize stMajor, TypeOfSize enMajor,
														 TypeOfSize stSub, TypeOfSize enSub) @safe pure nothrow @nogc
		in(stMajor < enMajor)
		in(stSub < enSub)
  	out(result; result[0] < TypeOfSize.max && result[1] < TypeOfSize.max){
			import std.algorithm: maxIndex;
			import std.container.array: Array;
			import std.math: abs;

		  immutable TypeOfSize LEN_MAJOR= cast(TypeOfSize)(enMajor-stMajor);
			immutable TypeOfSize LEN_SUB= cast(TypeOfSize)(enSub-stSub);
		  T[Size*Size] num= void;
			//num.reserve(LEN_SUB*LEN_MAJOR);
			//foreach(idxMajor; stMajor..enMajor) num.insertBack(a[idxMajor][stSub..enSub]);
			foreach(size_t i; 0u..LEN_MAJOR){
				num[i*LEN_SUB..(i+1)*LEN_SUB]= a[i][];
			}
			size_t idx= maxIndex!((a, b) => abs(a) < abs(b))(num[0u..LEN_SUB*LEN_MAJOR]);
			return [idx/LEN_SUB+stMajor, idx%LEN_SUB+stSub];
		}
	}
}
