/*****************************************************************************
 * Lower triangular matrix Impl
 *****************************************************************************/
module sekitk.lowertri;

import sekitk.base: TypeOfSize, MajorOrder;
/+
mixin template LowerTriangular(T, real Threshold, TypeOfSize Size, MajorOrder MatOdr){
+/
enum string LOWER_TRI= q{
	import std.range: dropOne;
	import sekitk.base;

public:
	/********************************************
	 * Operators
	 ********************************************/
	@safe pure nothrow const{
		// LowerTri pm Zero
		Matrix!(Row, Column, Shape) opBinary(string PM,
																				 MatrixType ShapeR: MatrixType.zero
																				 )(in Matrix!(Row, Column, ShapeR) rhs)
		if(isPlusOrMinusSign!PM){
			return new Matrix!(Row, Column, Shape)(this);
		}

		// LowerTri pm Diag
	  Matrix!(Row, Column, Shape) opBinary(string PM,
																				 MatrixType ShapeR: MatrixType.band1
																				 )(in Matrix!(Row, Column, ShapeR) rhs)
		if(isPlusOrMinusSign!PM){
			TypeOfInternalArray num= this.v[];
			TypeOfIndex j= 0;
			auto idxSet= IndexSetDiag!(Row, Column, Shape, MatOdr)();
			foreach(idx; idxSet) mixin("num[idx]" ~PM ~"= rhs._values[j++];");
			return new typeof(return)(num);
		}

		// LowerTri pm Dense, LowerTri pm TriDiag, LowerTri pm UpperTri
		Matrix!(Row, Column, MatrixType.dense) opBinary(string PM,
																										MatrixType ShapeR
																										)(in Matrix!(Row, Column, ShapeR) rhs)
		if(isPlusOrMinusSign!PM &&
			 (ShapeR is MatrixType.dense || ShapeR is MatrixType.upperTri)){
			enum MatrixType ShapeOfReturn= MatrixType.dense;
			T[arrayLength!(Row, Column, ShapeOfReturn)] num= void;

			// lhs
			{
				TypeOfIndex j= 0u;
				{
					auto idxSet= IndexSetStrictTriR!(Size, ShapeOfReturn, MatOdr)();
					foreach(idx; idxSet) num[idx]= VALUE_ZERO;
				}
				{
					auto temp0= IndexSetDiag!(Row, Column, ShapeOfReturn, MatOdr)();
					auto temp1= IndexSetStrictTriL!(Size, ShapeOfReturn, MatOdr)();
					foreach(idx; merge(temp0, temp1)) num[idx]= this._values[j++];
				}
			}
			// rhs
			static if(ShapeR is MatrixType.dense){
				foreach(TypeOfIndex idx; 0u..(num.length)) mixin("num[idx] " ~PM ~"= rhs._values[idx];");
			}
			else static if(ShapeR is MatrixType.band3){
				TypeOfIndex j= 0u;
				auto temp0= IndexSetDiag!(Row, Column, MatrixType.dense, MatOdr)();
				auto temp1= IndexSetSubDiag!(Size, MatrixType.dense, MatOdr)();
				foreach(idx; merge(temp0, temp1)) mixin("num[idx] " ~PM ~"= rhs._values[j++];");
			}
			else{	// UpperTri
				TypeOfIndex j= 0u;
				auto temp= IndexSetStrictTriR!(Size, ShapeOfReturn, MatOdr)();
				foreach(idx; temp) mixin("num[idx] " ~PM ~"= rhs._values[j++];");
			}

			return new typeof(return)(num);
		}

		// matrix multiplication
		/******************************************
		 * lower triangular * dense
		 ******************************************/
		Matrix!(Row, ColumnR,
						MatOpReturnType!(Op, Shape, ShapeR),
						MatOdr) opBinary(string Op: "*",
														 TypeOfSize ColumnR,
														 MatrixType ShapeR)(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs)
		if(ShapeR !is MatrixType.zero){
			enum TypeOfIndex LEN= arrayLength!(Row, ColumnR, MatOpReturnType!(Op, Shape, ShapeR));
		  T[LEN] num= void;

			final switch(ShapeR){
			case MatrixType.zero:
				assert(false);
			case MatrixType.dense:
				break;
			case MatrixType.band1:
				break;
			case MatrixType.band3:
				break;
			case MatrixType.upperTri:
				import mathematic.progression: sumFromZero;
				size_t idxPrd= 0, idxLhs= void, idxRhs= void;
				foreach(i; 0u..Row){
					idxLhs= sumFromZero(i);
					num[idxPrd++]= this._values[idxLhs]*rhs._values[0];
					foreach(j; 1u..ColumnR){
						idxRhs= j;
						num[idxPrd]= this._values[idxLhs]*rhs._values[idxRhs];
						foreach(k; 1u..((i+1u > j)? j+1u : i+1u)){
							idxRhs += ColumnR-k;
							num[idxPrd] += this._values[idxLhs+k]*rhs._values[idxRhs];
						}
						++idxPrd;
					}
				}
				break;
			case MatrixType.lowerTri:
			}

			return new typeof(return)(num);
		}
		unittest{	// lower triangular * upper triangular
			SekiTK!double.Matrix!(5, 5, MatrixType.lowerTri) lhsLowtri;
			SekiTK!double.Matrix!(5, 5, MatrixType.upperTri) rhsUptri;
			{
				immutable double[15] tempLowtri= [5,
																					 53, 61,
																					 67, 137, 241,
																					 337, 349, 383, 389,
																					 13, 557, 643, 701, 19];
				lhsLowtri= new typeof(lhsLowtri)(tempLowtri);
				immutable double[15] tempUptri= [929, 937, 1039, 1091, 1123,
																					1153, 1223, 1277, 1327,
																					1399, 1409, 1433,
																					1481, 11,
																					17];
				rhsUptri= new typeof(rhsUptri)(tempUptri);
			}
			immutable double[25] result= [4645, 4685, 5195, 5455, 5615,
																		49237, 119994, 129670, 135720, 140466,
																		62243, 220740, 574323, 587615, 602393,
																		313073, 718166, 1312787, 1929096, 1394692,
																		12077, 654402, 1594275, 2669640, 1683191];
			auto mul= lhsLowtri*rhsUptri;
			assert(mul._values[] == result[]);
		}
	}

private:
	@safe pure nothrow @nogc const{
		/******************************************
		 *
		 ******************************************/
		T[arrayLength!(Row, Column, MatrixType.dense)] arrayDense(){
			import std.algorithm: merge;

			typeof(return) num= void;
			{
				auto idxSet= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
				foreach(TypeOfIndex idx; idxSet) num[idx]= VALUE_ZERO;
			}
			{
				TypeOfIndex j= 0u;
				auto temp0= IndexSetDiag!(Row, Column, MatrixType.dense, MatOdr)();
				auto temp1= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
				foreach(idx; merge(temp0, temp1)) num[idx]= this._values[j++];
			}

			return num;
		}

		/******************************************
		 *
		 ******************************************/
		T detImpl(){
			typeof(return) result= _values[0];
			auto idxSet= IndexSetDiag!(Row, Column, Shape, MatOdr)();
			foreach(TypeOfIndex idx; idxSet.dropOne) result *= _values[idx];
			return result;
		}

		/******************************************
		 *
		 ******************************************/
		T[arrayLength!(Size, Size, MatrixType.lowerTri)] inverseImpl(){
			typeof(return) num= void;

			static if(Size == 1){
				num= [Identity!(T, "*")/_values[0]];
			}
			else static if(Size == 2){
				num= [_values[2], -_values[1], _values[0]];
				num[] /= this.det;
			}
			else static if(Size == 3){
				num= [_values[3]*_values[5], -_values[1]*_values[5], _values[1]*_values[4]-_values[2]*_values[3],
							_values[0]*_values[5], -_values[0]*_values[4],
							_values[0]*_values[3]];
				num[] /= this.det;
			}
			else assert(false, "FIXME: inverse of n > 3 lower triangular matrix");

			return num;
		}
	}
};
