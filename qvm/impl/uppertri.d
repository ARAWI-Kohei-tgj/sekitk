/*****************************************************************************
 * upper triangular matrix
 *****************************************************************************/
module sekitk.qvm.impl.uppertri;

import sekitk.qvm.common: TypeOfSize, MajorOrder;

/+
mixin template UpperTriangular(T, real Threshold,
															 TypeOfSize Size, MajorOrder MatOdr)
if(Size > 0){
+/
enum string UPPER_TRI= q{
	import std.algorithm: map;
	import std.range: dropOne;
	import sekitk.qvm.common;

	// Operators
	@safe pure nothrow const{
		/******************************************
		 * upperTri pm band1
		 ******************************************/
		Matrix!(Row, Column,
						MatrixType.band1,
						MatOdr) opBinary(string Op,
														 MatrixType ShapeR: MatrixType.band1
														 )(in Matrix!(Row, Column, ShapeR, MatOdr) rhs)
		if(isPlusOrMinusSign!Op){
			enum MatrixType ShapeOfReturn= MatOpReturnType!(Op, Shape, ShapeR);
			T[arrayLength!(Row, Column, ShapeOfReturn)] num;

			TypeOfIndex j= 0u;
			auto idxSet= IndexSetDiag!TemplateArgs();
			num= this._values[];
			foreach(TypeOfIndex i; idxSet) mixin("num[i]" ~Op ~"= rhs._values[j++];");

			return new typeof(return)(num);
		}

	  /******************************************
		 * upperTri pm others
		 ******************************************/
		Matrix!(Row, Column,
						MatrixType.dense,
						MatOdr) opBinary(string PM,
														 MatrixType ShapeR
														 )(in Matrix!(Row, Column, ShapeR, MatOdr) rhs)
		if(isPlusOrMinusSign!Op && (ShapeR is MatrixType.dense ||
																ShapeR is MatrixType.band3 ||
																ShapeR is MatrixType.lowerTri)){
			enum MatrixType ShapeDest= MatOpReturnType!(Op, Shape, ShapeR);
			T[arrayLength!(Row, Column, ShapeDest)] num;

			// lhs
			{
				{
					auto idxSet= IndexSetStrictTriL!(Size, ShapeDest, MatOdr)();
					foreach(idx; idxSet) num[idx]= VALUE_ZERO;
				}
				{
					TypeOfIndex j= 0u;
					auto temp0= IndexSetDiag!(Row, Column, ShapeDest, MatOdr)();
					auto temp1= IndexSetStrictTriR!(Size, ShapeDest, MatOdr)();
					foreach(idx; merge(temp0, temp1)) num[idx]= this._values[j++];
				}
			}
			// rhs
			static if(ShapeR is MatrixType.dense){
				foreach(TypeOfIndex idx; 0u..(num.length)) mixin("num[idx]=" ~PM ~"rhs._values[idx];");
			}
			else static if(ShapeR is MatrixType.band3){
				TypeOfIndex j= 0u;
				auto temp0= IndexSetDiag!(Row, Column, MatrixType.dense, MatOdr)();
				auto temp1= IndexSetSubDiag!(Size, MatrixType.dense, MatOdr)();
				foreach(idx; merge(temp0, temp1)) mixin("num[idx] " ~PM ~"= rhs._values[j++];");
			}
			else{	// LowerTri
				TypeOfIndex j= 0u;
				auto temp= IndexSetStrictTriL!(Size, ShapeDest, MatOdr)();
				foreach(idx; temp) mixin("num[idx]= " ~PM ~"rhs._values[j++];");
			}

			return new typeof(return)(num);
		}

		/******************************************
		 * upperTri * band1
		 ******************************************/
		Matrix!(Row, ColumnR,
						Shape,
						MatOdr) opBinary(string Op: "*",
														 TypeOfSize ColumnR,
														 MatrixType ShapeR: MatrixType.band1
														 )(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs){
			import std.range: zip;
			T[arrayLength!(Row, ColumnR, Shape)] num= void;
			TypeOfIndex k= 0;
			T[] elmLhs, elmRhs;
		  auto idxSet= IndexSetDiag!(Size, Size, Shape, MatOdr)();

			foreach(i; 0u..Size){
				elmLhs= _values[idxSet.front .. idxSet.front+Size-i].dup;
				idxSet.popFront;
				elmRhs= _values[i..$].dup;
				foreach(numL, numR; zip(elmLhs, elmRhs)) num[k++]= numL+numR;
			}

			return new typeof(return)(num);
		}

		/******************************************
		 * uptri mul other
		 ******************************************/
		Matrix!(Row, ColumnR,
						MatOpReturnType!(Op, Shape, ShapeR),
						MatOdr) opBinary(string Op: "*",
														 TypeOfSize ColumnR,
														 MatrixType ShapeR
														 )(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs){
			enum ShapeDest= MatOpReturnType!(Op, Shape, ShapeR);
			T[arrayLength!(Row, ColumnR, ShapeDest)] num= void;
			size_t idxPrd= void, idxLhs= void, idxRhs= void;

			final switch(ShapeR){
			case MatrixType.dense:
				alias idxMapDest= indexMap!(Row, ColumnR, ShapeDest, MatOdr);
				idxLhs= Column-1u;
				foreach(i; iota!TypeOfSize(Row)){
					foreach(j; iota!TypeOfSize(ColumnR)){
						idxPrd= idxMapDest(i, j);
						idxRhs= (Column-1u)*ColumnR+j;
						num[idxPrd]= this._values[idxLhs]*rhs._values[idxRhs];
						foreach(k; 1u..Column-i) num[idxPrd] += this._values[idxLhs-k]*rhs._values[idxRhs-k*ColumnR];
					}
					idxLhs += Column-(i+1u);
				}
				break;
			case MatrixType.band3:
				break;
			case MatrixType.upperTri:
				auto idxsetLhs= IndexSetDiag!TmplArgsOfThis;
				auto idxsetRhs= IndexSetDiag!(Row, ColumnR, ShapeR, MatOdr)();
				foreach(i; 0u..Row){
					foreach(j; 0u..ColumnR-i){
						idxPrd= idxsetLhs[j]+i;
						idxLhs= idxsetLhs[j];
						num[idxPrd]= this._values[idxLhs]*rhs._values[idxLhs+i];
						foreach(k; 1u..i+1u) num[idxPrd] += this._values[idxLhs+k]*rhs._values[idxsetRhs[j+k]+(i-1u)-(k-1u)];
					}
				}
				break;
			case MatrixType.lowerTri:
/+
			 foreach(size_t i; ){
			 foreach(size_t j;){
			 idx_prd= ;
			 num[idx_prd]= v.[]*rhs.v.[];
			 foreach(size_t k;);
			 }
			 }+/

					/+
 0  1  2  3  4 |  0             |  0  1  2  3  4
    5  6  7  8 |  1  2          |  5  6  7  8  9
       9 10 11 |  3  4  5       | 10 11 12 13 14
         12 13 |  6  7  8  9    | 15 16 17 18 19
            14 | 10 11 12 13 14 | 20 21 22 23 24

c0= a0*b0 +a1*b1  +a2*b3  +a3*b6 +a4*b10
c1=        a1*b2  +a2*b4  +a3*b7 +a4*b11
c2=                a2*b5  +a3*b8 +a4*b12
c3=                        a3*b9 +a4*b13
c4=                               a4*b14

c5=        a5*b1  +a6*b3  +a7*b6 +a8*b10
c6=        a5*b2  +a6*b4  +a7*b7 +a8*b11
c7=                a6*b5  +a7*b8 +a8*b12
c8=                        a7*b9 +a8*b13
c9=                               a8*b14

c10=               a9*b3 +a10*b6 +a11*b10
c11=               a9*b4 +a10*b7 +a11*b11
c12=               a9*b5 +a10*b8 +a11*b12
c13=                      a10*b9 +a11*b13
c14=                              a11*b14

c15=                      a12*b6 +a13*b10
c16=                      a12*b7 +a13*b11
c17=                      a12*b8 +a13*b12
c18=                      a12*b9 +a13*b13
c19=                              a13*b14

c20=                              a14*b10
c21=                              a14*b11
c22=                              a14*b12
c23=                              a14*b13
c24=                              a14*b14
+/
				break;
			case MatrixType.zero, MatrixType.band1:
			  assert(false);
			}

			return new typeof(return)(num);
		}

		unittest{
			SekiTK!double.Matrix!(5, 5, MatrixType.upperTri) lhsUptri;
			SekiTK!double.Matrix!(5, 2, MatrixType.dense) rhsDense;
			SekiTK!double.Matrix!(5, 5, MatrixType.upperTri) rhsUptri;
			{
				immutable double[15] temp_uptri= [2.0, 5, 11, 17, 23,
																					31, 41, 47, 59,
																					67, 73, 83,
																					97, 103,
																					109];
				lhsUptri= new typeof(lhsUptri)(temp_uptri);
				immutable double[10] temp_dense= [113, 107,
																					101, 89,
																					79, 71,
																					61, 53,
																					43, 37];
				rhsDense= new typeof(rhsDense)(temp_dense);
				immutable double[15] temp_uptri2= [113, 107, 101, 89, 79,
																					 71, 61, 53, 43,
																					 37, 29, 19,
																					 13, 7,
																					 3];
				rhsUptri= new typeof(rhsUptri)(temp_uptri2);
			}

			{	// UpperTri * Dense
				immutable double[10] result= [3626, 3192,
																			11774, 10344,
																			13315, 11697,
																			10346, 8952,
																			4687, 4033];
				auto mul= lhsUptri*rhsDense;
				assert(mul._values[] == result[]);
			}

			// upper triangular * diagonal
			// upper triangular * tridiagonal

			{	// upper triangular * upper triangular
				immutable double[15] result= [226, 569, 914, 983, 770,
																			2201, 3408, 3443, 2618,
																			2479, 2892, 2033,
																			1261, 988,
																			327];
				auto mul= lhsUptri*rhsUptri;
				assert(mul._values[] == result[]);
			}

			// upper triangular * lower triangular
		}
	}

package:
	@safe pure nothrow @nogc const{
		/******************************************
		 * Generating the converting array of upper triangular to dense
		 ******************************************/
		T[arrayLength!(Row, Column, MatrixType.dense)] arrayDense(){
			import std.algorithm: merge;

			typeof(return) num= void;
			{
				auto idxSet= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
			  foreach(TypeOfIndex idx; idxSet) num[idx]= VALUE_ZERO;
			}
			{
				TypeOfIndex j= 0u;
				auto temp0= IndexSetDiag!(Row, Column, MatrixType.dense, MatOdr)();
				auto temp1= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
				foreach(idx; merge(temp0, temp1)) num[idx]= this._values[j++];
			}

			return num;
		}

		/****************************************
		 * extracts the specified column
		 *
		 * Params:
		 * 	idxColumn= index of column (0 ≤ idxColumn < Column)
		 ****************************************/
		T[Row] sliceColumn(IdxType)(in IdxType idxColumn)
		if(isIntegral!IdxType)
	  in(idxColumn >= 0u && idxColumn < Column){
			typeof(return) result= void;

			auto idxSet= (){
				final switch(MatOdr){
				case MajorOrder.row:
			  	return recurrence!((a, n) => a[n-1]+Column-n)(idxColumn).take(idxColumn+1u);
					break;
				case MajorOrder.column:
					enum TypeOfIndex IDX_ST= sumFromZero(idxColumn);
					return iota(IDX_ST, IDX_ST+idxColumn+1u);
				}
			}();

			foreach(idxDest, idxInternal;
							zip(iota(idxColumn+1u), idxSet)) result[idxDest]= this._values[idxInternal];
			result[idxColumn+1u..$]= VALUE_ZERO;
			return result;
	  }
		@safe pure @nogc unittest{
			double[10] num= [4.75, 3.0, -2.25, 8.675,
											 -2.825, -6.125, -5.5,
											 -7.0, 1.25,
											 5.625];
			auto foo= new SekiTK!double.Matrix!(4, 4, MatrixType.upperTri, MajorOrder.row)(num);
			assert(foo.sliceColumn(0u) == [4.75, 0.0, 0.0, 0.0]);
			assert(foo.sliceColumn(1u) == [3.0, -2.825, 0.0, 0.0]);
			assert(foo.sliceColumn(2u) == [-2.25, -6.125, -7.0, 0.0]);
			assert(foo.sliceColumn(3u) == [8.675, -5.5, 1.25, 5.625]);
			/+
 0  1  2  3
 *  4  5  6
 *  *  7  8
 *  *  *  9
+/
		}
		@safe pure @nogc unittest{
			double[10] num= [4.75,
											 3.0, -2.25,
											 8.675, -2.825, -6.125,
											 -5.5, -7.0, 1.25, 5.625];
			auto foo= new SekiTK!double.Matrix!(4, 4, MatrixType.upperTri, MajorOrder.column)(num);
			assert(foo.sliceColumn(0u) == [4.75, 0.0, 0.0, 0.0]);
			assert(foo.sliceColumn(1u) == [3.0, -2.25, 0.0, 0.0]);
			assert(foo.sliceColumn(2u) == [8.675, -2.85, -6.125, 0.0]);
			assert(foo.sliceColumn(3u) == [-5.5, -7.0, 1.25, 5.625]);
			/+
 0  1  3  6
 *  2  4  7
 *  *  5  8
 *  *  *  9
+/
		}

		/******************************************
		 * determinant
		 *
		 * Time complexity:
		 * 	O(n)= n
		 ******************************************/
		T detImpl(){
			typeof(return) result= _values[0];
			auto idxSet= IndexSetDiag!(Row, Column, Shape, MatOdr)();
			foreach(TypeOfIndex idx; idxSet.dropOne) result *= _values[idx];
			return result;
		}

		/******************************************
		 * Internal array of inverse matrix
		 *
		 * decompose U= D+U_s
		 * U= D(I+D^{-1} U_s)
		 * U^{-1}= (I+D^{-1} U_s)^{-1} D^{-1}
		 * =\right[ sum_{i=0}^{n-1} (-D^{-1} U_s)^i \left] D^{-1}
		 *
		 * Time complexity:
		 * 	O(n)= n^2+3n-6
		 ******************************************/
		TypeOfInternalArray inverseImpl(){
			typeof(return) result= void;
			{
				auto idxSetDiag= IndexSetDiag!TemplateArgs();
				foreach(idx; idxSetDiag) result[idx]= Identity!(T, "*")/this._values[idx];
			}

			static if(Size == 1){}
			else static if(Size == 2){
				result[1]= -_values[1]/this.det;
			}
			else static if(Size == 3){
				final switch(MatOdr){
				case MajorOrder.row:
					result[1]= -_values[1]/(_values[0]*_values[3]);
					result[4]= -_values[4]/(_values[3]*_values[5]);
					result[2]= (_values[1]*_values[4]-_values[2]*_values[3])/this.det;
					break;
				case MajorOrder.column:
					result[1]= -_values[1]/(_values[0]*_values[2]);
					result[4]= -_values[4]/(_values[2]*_values[5]);
					result[3]= (_values[1]*_values[4]-_values[2]*_values[3])/this.det;
				}
			}
			else static if(Size == 4){
				final switch(MatOdr){
				case MajorOrder.row:
					result[1]= -_values[1]/(_values[0]*_values[4]);
					result[2]= (_values[2]-_values[1]*_values[5]/_values[4])/(_values[0]*_values[7]);
					result[3]= ((_values[3]*_values[4]-_values[1]*_values[6])*_values[7]
											-(_values[2]*_values[4]+_values[1]*_values[5])*_values[8])/this.det;
					result[5]= _values[5]/(_values[4]*_values[7]);
					result[6]= (_values[6]-_values[5]*_values[8]/_values[7])/(_values[4]*_values[9]);
					result[8]= -_values[8]/(_values[7]*_values[9]);
					break;
				case MajorOrder.column:
					result[1]= -_values[1]/(_values[0]*_values[2]);
					result[3]= (_values[3]-_values[1]*_values[4]/_values[2])/(_values[0]*_values[5]);
					result[4]= _values[4]/(_values[2]*_values[5]);
					result[6]=((_values[2]*_values[6]-_values[1]*_values[7])*_values[5]
										 -(_values[2]*_values[3]+_values[1]*_values[4])*_values[8])/this.det;
					result[7]= (_values[7]-_values[4]*_values[8]/_values[5])/(_values[2]*_values[9]);
					result[8]= -_values[8]/(_values[5]*_values[9]);
				}
			}
			else{
				assert(false, "I shall be implemented.");	// FIXME:
			}

			return result;
		}
	}
};
