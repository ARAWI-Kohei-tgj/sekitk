/**************************************************************
 * diagonal matrix
 *************************************************************/
module sekitk.band1;

import sekitk.base: TypeOfSize, MajorOrder;

mixin template Band1(T, real Thresold, TypeOfSize Size, MajorOrder MatOdr)
if(Size > 1u){
	import std.traits: Unqual;
	import sekitk.base;
	import sekitk.indexset: IndexSetDiag;

public:
	/******************************
	 * Operators
	 *****************************/
	@safe pure nothrow const{
		/++++++++++++++++++++++
		 + add & sub operator
		 +++++++++++++++++++++/
		Matrix!(Size, Size,
						ShapeR, MatOdr) opBinary(string Op,
																		 MatrixType ShapeR)(in Matrix!(Size, Size,
																																	 ShapeR, MatOdr) rhs)
		if(isPlusOrMinusSign!Op && ShapeR !is Shape){
			return mixin("rhs" ~Op ~"this;");
		}

		/++++++++++++++++++++++
		 + matrix multiplications
		 +++++++++++++++++++++/
		Matrix!(Row, ColumnR,
						ShapeR, MatOdr) opBinary(string OP: "*",
																		 TypeOfSize ColumnR,
																		 MatrixType ShapeR)(in Matrix!(Column, ColumnR,
																																	 ShapeR, MatOdr) rhs) @nogc
		if(ShapeR is MatrixType.dense && ColumnR > 0u ||
			 (ShapeR is MatrixType.band1 || ShapeR == MatrixType.band3 ||
				ShapeR is MatrixType.upperTri || ShapeR == MatrixType.lowerTri) &&
			 ColumnR == Column){
			T[arrayLength!(Row, ColumnR, ShapeR)] num= void;

			final switch(ShapeR){
			case MatrixType.zero:
				break;
			case MatrixType.dense:
				size_t st= 0u, en= void;
				foreach(i; 0u..Row){
					en= (i+1)*ColumnR;
					num[st..en]= rhs.v[st..en]*this.v[i];
					st= en;
				}
				break;
			case MatrixType.band1:
				num[]= this.v[]*rhs.v[];
				break;
			case MatrixType.band3:
				size_t st= 0u, en= void;
				foreach(i; 0u..Row){
					en= 3u*i +(i < Row-1u)? 2u : 1u;
					num[st..en]= rhs.v[st..en]*v[i];
					st= en;
				}
				break;
			case MatrixType.upperTri:
				size_t st= 0u, en= ColumnR;
				foreach(i; 0u..Row){
					num[st..en]= rhs.v[st..en]*this.v[i];
					st= en;
					en += (ColumnR-(i+1u));
				}
				break;
			case MatrixType.lowerTri:
				size_t st= 0u, en= 1u;
				foreach(i; 0u..Row){
					num[st..en]= rhs.v[st..en]*this.v[i];
					st= en;
					en += i+2u;
				}
			}
			return typeof(return)(num);
		}

		/++++++++++++++++++++++
		 + power
		 +++++++++++++++++++++/
		Matrix!(Row, Column,
						MatrixType.band1,
						MatOdr) opBinary(string Pow: "^^", N)(in N rhs) @nogc
	  if(isUnsignedInt!N){
			typeof(return) result;
			if(rhs == 0u){
				result= typeof(return)(Identity!(T, "*"));
			}
			else{
				T[arrayLength!(Row, Column, MatrixType.band1)] num;
				num[]= this.v[]^^rhs;
				result= typeof(return)(num);
			}
			return result;
		}

		/++++++++++++++++++++++
		 + cast
		 +++++++++++++++++++++/
	  Matrix!(Row, Column,
						ShapeDest,
						MatOdrDest) opCast(MatrixType ShapeDest,
															 MajorOrder MatOdrRet)() @nogc
		if(ShapeDest is MatrixType.band3 ||
			 ShapeDest is MatrixType.upperTri ||
			 ShapeDest is MatrixType.lowerTri){
			T[arrayLength!(Row, Column, ShapeDest)] num;

			static if(MatTypRet is MatrixType.band3){
				auto idxset= IndexSetSubDiag!(Size, ShapeDest, MatOdrDest)();
			}
			else static if(MatTypRet is MatrixType.upperTri){
				auto idxset= IndexSetStrictTriR!(Size, ShapeDest, MatOdrDest)();
			}
			else static if(MatTypRet is MatrixType.lowerTri){
				auto idxset= IndexSetStrictTriL!(Size, ShapeDest, MatOdrDest)();
			}
			foreach(i; idxset) num[i]= VALUE_ZERO;

			{
				auto idxset= IndexSetDiag!(Row, Column, ShapeDest, MatOdrDest)();
				size_t j= 0u;
				foreach(i; idxset) num[i]= this.v[j++];
			}
			return typeof(return)(num);
		}
	}

	/******************************
	 * Other methods
	 *****************************/
	@safe pure nothrow const{
		// transpose
		Matrix!(Column, Row, Shape, MatOdr) transposed() @nogc @property{
			T[arrayLength!(Column, Row, Shape)] num= void;
			static if(isFloatingPoint!T) import numeric.pseudocmplx: conj;

			foreach(idx; 0u..Size) num[idx]= v[idx].conj;
			return typeof(return)(num);
		}
	}

private:
	import std.traits: isIntegral, isUnsigned;

	@safe pure nothrow const{
		/++++++++++++++++++++++++++++++
		 + This function generates the internal array of Matrix!()
		 ++++++++++++++++++++++++++++++/
		T[arrayLength!(Size, Size, MatrixType.dense)] arrayDense() @nogc{
			typeof(return) temp= VALUE_ZERO;
			size_t j= 0u;
		  auto idxSet= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
			foreach(TypeOfIndex i; idxSet) temp[i]= this.v[j++];
			return temp;
		}

		/++
		 +
		 +++/
		T detImpl() @nogc{
			typeof(return) result= v[0];
			foreach(num; v[1u..$]) result *= num;
			return result;
		}

		/++
		 +
		 +++/
	  TypeOfInternalArray inverseImpl() @nogc{
			typeof(return) num= v[];
			num[]= 1.0/num[];
			return num;
		}
	}
}
