/******************************************************************************
 * dense matrix
 ******************************************************************************/
module sekitk.dense;

import sekitk.base: TypeOfSize, MajorOrder;

/**************************************************************
 * common methods of dense matrix both rectangular and square
 **************************************************************/
mixin template DenseCommon(T, TypeOfSize Row, TypeOfSize Column, MajorOrder MatOdr)
if(Row > 0 && Column > 0){
	import std.traits: Unqual, ParameterTypeTuple;
	import std.array: array;
	import std.range: repeat, iota, zip;
	import mathematic.basic: Identity;

	import sekitk.base;
	import sekitk.indexset;

private:
	/********************************************
	 * Operators
	 ********************************************/
	@safe pure nothrow const{
		// add & sub (different shape)
		TypeOfThis opBinary(string Op,
												MatrixType ShapeR: MatrixType.band1
												)(in Matrix!(Row, Column,
																		 ShapeR, MatOdr) rhs) @nogc
		if(isPlusOrMinusSign!Op){
		  TypeOfInternalArray num= this.v[];
			TypeOfIndex j= 0u;
			auto idxSet= IndexSetDiag!TmplArgsOfThis();
			foreach(i; idxSet) mixin("num[i] " ~PM ~"= rhs.v[j++];");

			return typeof(return)(num);
		}

		// other shape
		TypeOfThis opBinary(string Op,
												MatrixType ShapeR
												)(in Matrix!(Row, Column,
																		 ShapeR, MatOdr) rhs) @nogc
		if(isPlusOrMinusSign!Op && ShapeR != Shape){
			import std.algorithm: merge;
			TypeOfInternalArray num= this.v[];
			TypeOfIndex j= 0u;

			auto idxSet0= IndexSetDiag!TmplArgsOfThis();

			static if(ShapeR is MatrixType.band3){
				auto idxSet1= IndexSetSubDiag!(Size, Shape, MatOdr)();
				foreach(idx; merge(idxSet0, idxSet1)) mixin("num[idx] " ~PM ~"= rhs.v[j++];");
			}
			else static if(ShapeR is MatrixType.upperTri){
				auto idxSet1= IndexSetStrictTriR!(Size, Shape, MatOdr)();
				foreach(idx; merge(idxSet0, idxSet1)) mixin("num[idx] " ~PM ~"= rhs.v[j++];");
			}
			else static if(ShapeR is MatrixType.lowerTri){
				auto idxSet1= IndexSetStrictTriL!(Size, Shape, MatOdr)();
				foreach(idx; merge(idxSet0, idxSet1)) mixin("num[idx] " ~PM ~"= rhs.v[j++];");
			}
			else{
				assert(false);
			}

			return typeof(return)(num);
		}

		// matrix multiplications
		Matrix!(Row, ColumnR,
						MatOpReturnType!(Op, Shape, ShapeR, MatOdr)
						) opBinary(string Op: "*",
											 TypeOfSize ColumnR,
											 MatrixType ShapeR)(in Matrix!(Column, ColumnR,
																										 ShapeR, MatOdr) rhs) @nogc{
			enum TypeOfIndex LEN= arrayLength!(Row, ColumnR, MatOpReturnType!(OP, Shape, ShapeR));
			T[LEN] num= void;

			final switch(ShapeR){
			case MatrixType.dense:
				alias idxMapDest= indexMap!(templateArgsOf!(typeof(return)));
				alias idxMapLhs= indexMap!TmplArgsOfThis;
				alias idxMapRhs= indexMap!(templateArgsOf!(typeof(rhs)));
				foreach(i; 0u..Row)
					foreach(j; 0u..ColumnR){{
						num[idxMapDest(i, j)]= this.v[idxMapLhs(i, 0u)] * rhs.v[idxMapRhs(j, 0u)];
						foreach(k; 1u..Column) num[idxMapDest(i, j)] += this.v[idxMapLhs(i, k)]*rhs.v[idxMapRhs(k, j)];
					}}
				break;
			case MatrixType.band1:
				alias idxMapLhs= indexMap!TmplArgsOfThis;
				TypeOfIndex k= 0u;
				foreach(i; 0u..Row)
					foreach(j; 0u..Column){{
					  num[k++]= this.v[idxMapLhs(i, j)]*rhs.v[j];
			  }}
				break;
			case MatrixType.band3:
				size_t idxPrd= 0u, idxLhs= void, idxRhs= void, en= void;
				foreach(i; 0u..Row){
					idxLhs= i*Column;
					idxRhs= 0u;
					foreach(j; 0u..ColumnR){
						en= (j == 0u || j == ColumnR-1u)? 2u : 3u;
						num[idxPrd]= this.v[idxLhs]*rhs.v[idxRhs];
						foreach(k; 1u..en) num[idxPrd] += this.v[idxLhs+k]*rhs.v[idxRhs+2u*k];
						++idxPrd;
						if(j == 0u){
							idxRhs += 1u;
						}
						else{
							++idxLhs;
							idxRhs += 3u;
						}
					}
				}
				break;
			case MatrixType.upperTri:
				size_t temp= void;
				alias idxMapDest= indexMap!(templateArgsOf!(typeof(return)));
				alias idxMapLhs= indexMap!TmplArgsOfThis;
				foreach(i; iota!TypeOfSize(ColumnR)){
					foreach(j; iota!TypeOfSize(Row)){
						num[idxMapDest(j, i)]= this.v[idxMapLhs(j, 0u)]*rhs.v[i];
						temp= 0u;
						foreach(k; 1u..i+1u){
							num[idxMapDest(j, i)] += this.v[idxMapLhs(j, k)]*rhs.v[i-1u+ColumnR+temp];
							temp += ColumnR-(k+1);
						}
					}
				}
				break;
			case MatrixType.lowerTri:
				import mathematic.progression: sumFromZero;

				size_t indexRhs(in size_t i, in size_t k) @safe pure nothrow @nogc const{
					return i*(i+3u)/2u+sumFromZero(k)-sumFromZero(i);
				}

				size_t idxPrd= void, idxRhs= void, temp= void;
				foreach(i; 0u..Column){
					foreach(j; 0u..Row){
						idxPrd= j*ColumnR+i;
						idxRhs= indexRhs(i, ColumnR-1u);
						num[idxPrd]= this.v[Column*(j+1u)-1u]*rhs.v[idxRhs];
						foreach_reverse(k; i..ColumnR-1u) num[idxPrd] += this.v[j*Column+k]*rhs.v[indexRhs(i, k)];
					}
				}
				break;
			case MatrixType.zero:
				assert(false);
			}
			return new typeof(return)(num);
		}
		unittest{
			immutable float[21] tempDense1= [2.0, 5, 11, 17, 23, 31, 41,
																				47, 59, 67, 73, 83, 97, 103,
																				109, 127, 137, 149, 157, 167, 173];
			immutable float[35] tempDense2= [3.0, 7, 13, 19, 29,
																				37, 43, 53, 61, 71,
																				79, 89, 101, 107, 113,
																				131, 139, 151, 163, 179,
																				181, 191, 193, 197, 199,
																				211, 223, 227, 229, 233,
																				239, 241, 251, 257, 263];
			immutable float[15] result= [23790.0, 24758, 25736, 26458, 27282,
																	 77287, 81283, 85419, 88595, 92397,
																	 140369, 148049, 156117, 162397, 169983];

			auto mDense1= new SekiTK!float.Matrix!(3, 7, MatrixType.dense)(tempDense1);
			auto mDense2= new SekiTK!float.Matrix!(7, 5, MatrixType.dense)(tempDense2);
			auto mul= mDense1*mDense2;
			assert(mul.v[] == result[]);
		}
		unittest{
			SekiTK!float.Matrix!(4, 5, MatrixType.dense) lhsDense;
			SekiTK!float.Matrix!(5, 5, MatrixType.band1) rhsDiag;
			SekiTK!float.Matrix!(5, 5, MatrixType.band3) rhsBand3;
			SekiTK!float.Matrix!(5, 5, MatrixType.upperTri) rhsUptri;
			SekiTK!float.Matrix!(5, 5, MatrixType.lowerTri) rhsLotri;
			{
				immutable float[20] tempDense= [2.0, 3, 5, 7, 11,
																				 13, 17, 19, 23, 29,
																				 31, 37, 41, 43, 47,
																				 53, 59, 61, 67, 71];
				lhsDense= new typeof(lhsDense)(tempDense);
				immutable float[5] tempDiag= [73, 79, 83, 89, 97];
				rhsDiag= new typeof(rhsDiag)(tempDiag);
				immutable float[13] tempBand3= [157, 151,
																				 149, 139, 137,
																				 131, 127, 113,
																				 107,101, 89,
																				 79, 73];
				rhsBand3= new typeof(rhsBand3)(tempBand3);
				immutable float[15] tempUptri=[73.0, 79, 83, 89, 97,
																				101, 103, 107, 109,
																				113, 127, 131,
																				137, 139,
																				149];
				rhsUptri= new typeof(rhsUptri)(tempUptri);
				immutable float[15] tempLotri=[73.0,
																				79, 83,
																				89, 97, 101,
																				103, 107, 109, 113,
																				127, 131, 137, 139, 149];
				rhsLotri= new typeof(rhsLotri)(tempLotri);
			}
			// dense * diagonal
			{
				immutable float[20] result= [146, 237, 415, 623, 1067,
																		 949, 1343, 1577, 2047, 2813,
																		 2263, 2923, 3403, 3827, 4559,
																		 3869, 4661, 5063, 5963, 6887];
				auto mul= lhsDense*rhsDiag;
				assert(mul.v[] == result[]);
			}
			// dense * band3
			{
				immutable float[20] result= [761, 1374, 1795, 2141, 1426,
																		 4574, 6815, 7203, 6761, 4164,
																		 10380, 15195, 14877, 12689, 7258,
																		 17112, 24195, 22999, 19269, 11146];
				auto mul= lhsDense*rhsBand3;
				assert(mul.v[] == result[]);						
			}
			// dense * upper triangular
			{
				immutable float[20] result= [146.0, 461, 1040, 2093, 3788,
																		 949, 2744, 4977, 8540, 13121,
																		 2263, 6186, 11017, 17816, 25391,
																		 3869, 10146, 17369, 27956, 39455];
				auto mul= lhsDense*rhsUptri;
				assert(mul.v[] == result[]);
			}
			// dense * lower triangular
			{
				immutable float[20] result= [2946.0, 2924, 2775, 2320, 1639,
																		 10035, 9514, 8399, 6630, 4321,
																		 19233, 17806, 15267, 11392, 7003,
																		 29877, 27284, 23191, 17440, 10579];
				auto mul= lhsDense*rhsLotri;
				assert(mul.v[] == result[]);
			}
		}	// end of unittest
	}

  const{
		int opApply(Dg)(scope Dg dg)
		if(ParameterTypeTuple!Dg.length == 1){
			typeof(return) result= 0;

			foreach(elm; v){
				result= dg(elm);
			}
			return result;
		}

		int opApply(Dg)(scope Dg dg)
		if(ParameterTypeTuple!Dg.length == 1){
			typeof(return) result= 0;

			foreach(idx; 0..Column){
				result= dg(Vector!Row(sliceColumn(idx)));
			}
			return result;
		}
	}
}


/*************************************************************
 * Rectangular dense matrix
 *************************************************************/
mixin template RectDense(T, real Threshold,
												 TypeOfSize Row,
												 TypeOfSize Column,
												 MajorOrder MatOdr)
if(Row > 0 && Column > 0){
	mixin DenseCommon!(T, Row, Column, MatOdr);

	version(future){
		void pseudoInverse(){};
	}
}


/*************************************************************
 * Square dense matrix
 *************************************************************/
mixin template SqDense(T, real Threshold, TypeOfSize Size, MajorOrder MatOdr)
if(Size > 0){
	mixin DenseCommon!(T, Size, Size, MatOdr);

private:
	/****************************
	 * determinant
	 ****************************/
	static if(Size == 1){
		T detImpl() @safe pure nothrow @nogc const{return v[0];}
	}
	else static if(Size == 2){
		/// ditto
		/// Standars: rule of Sarrus
		T detImpl() @safe pure nothrow @nogc const{
			typeof(return) result;

			final switch(MatOdr){
			case MajorOrder.row, MajorOrder.column:
				result= v[0]*v[3]-v[1]*v[2];
				break;
			case MajorOrder.diag:
				result= v[0]*v[1]-v[2]*v[3];
			}

			return result;
		}
	}
	else static if(Size == 3){
		/// ditto
		/// Standards: rule of Sarrus
		T detImpl() @safe pure nothrow @nogc const{
			typeof(return) result;

			final switch(MatOdr){
			case MajorOrder.row, MajorOrder.column:
				result= v[0]*(v[4]*v[8]-v[5]*v[7])
					+v[1]*(v[5]*v[6]-v[3]*v[8]) +v[2]*(v[3]*v[7]-v[4]*v[6]);
				break;
			case MajorOrder.diag:
				result= v[0]*(v[1]*v[2]-v[4]*v[6])
					+v[3]*(v[4]*v[8]-v[2]*v[5]) +v[7]*(v[5]*v[6]-v[1]*v[8]);
			}

			return result;
		}
	}
	else{
		/// ditto
		T detImpl() @safe pure{
			import numeric.pseudocmplx: assumeRealNum;
			typeof(return) result;

			if(ev !is null){
				static if(isComplex!T) result= ev.det;
				else result= assumeRealNum!(CT, Threshold)(ev.det);
			}
			else if(lup !is null) result= lup.det;	// LU decomposition is obtained
			else{
				auto temp= LUdecomposition!(T, Threshold, Size, Shape, AlgoLU, MatOdr).forwardElimination(this.opCast!(T[Column][Row])());	// forward elimination
				// if DIP 1008 is enable, @nogc
				result= temp.det;
			}

			return result;
		}
	}

	/****************************
	 * inverse
	 ****************************/
	static if(Size < 3){
		TypeOfInternalArray inverseImpl() @safe pure nothrow @nogc const{
			typeof(return) result= void;

			static if(Size == 1){
			  result= [Identity!(T, "*")/v[0]];
			}
			else{
				final switch(MatOdr){
				case MajorOrder.row, MajorOrder.column:
					result= [v[3], -v[1], -v[2], v[0]];
					break;
				case MajorOrder.diag:
					result= [v[1], v[0], -v[2], -v[3]];
				}
				result[] /= this.detImpl;
			}
			
			return result;
		}
	}
	else{
		TypeOfInternalArray inverseImpl() @safe pure{
			typeof(return) num= void;
			if(lup is null) this.decompose!(DecompScheme.lup);

			if(lup !is null && lup.isInvertible) num= lup.inverseArray;

			return num;
		}
	}

	/+
a11, a12
a21, a22

1,   0| u11 u12| u11      u12
l21, 1|   0 u22| l21*u11  l21*u12+u22

u11= a11
u12= a12
l21*u11= a21
    l21= a21/a11
l21*u12+u22= a22
u22= a22-a21*a12/a11
+/
}
