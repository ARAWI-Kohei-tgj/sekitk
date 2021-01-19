module sekitk.qvm.algorithm.decomposition.lu_dummy;

enum string LU_BAND1_IMPL= q{
	enum size_t LEN= arrayLength!(Size, Size, MatrixType.band1);
	Result!Method decompose(DecompScheme Method)(in T[LEN] num){
		return typeof(return)(num,
													num[].all!(a => !approxEqualZero!(T, Threshold)(a)));
	}

	struct Result(DecompScheme Method){
		@safe pure:
		// Constructors
		nothrow @nogc{
			/// Copy constructor
			this(ref return scope inout){}

			///
			this(ArrayType)(in ArrayType num, in bool invertibility)
			if((isDynamicArray!ArrayType && is(ArrayType: T[])) ||
				 (isStaticArray!ArrayType && is(ArrayType: T[LEN]))){
				_invertible= invertibility;
				_values[]= num[];
			}
		}

		const{
			bool isInvertible() nothrow @nogc{
				return _invertible;
			}

			/****************************************
			 * Determinant
			 ****************************************/
			T det() nothrow @nogc{
				typeof(return) result= _values[0];
				auto idxSet= IndexSetDiag!(Size, Size, Shape, MatOdr)();
				foreach(idx; idxSet.dropOne) result *= _values[idx];
				return result;
			}

			/****************************************
			 * Inverse
			 *
			 * Throws:
			 *	SingularMatrix
			 ****************************************/
			T[LEN] inverseArry(){
				typeof(return) result= void;
				foreach(idx; 0u..LEN) result[idx]= 1.0/_values[idx];
				return result;				
			}

			/****************************************
			 * Lower triangular matrix $(I L)
			 ****************************************/
			T[arrayLength!(Size, Size, MatrixType.band1)]
			matrix(DecomposedMat Mat: DecomposedMat.lowerTri)() nothrow @nogc{

			}
		}

	private:
		bool _invertible;
		T[LEN] _values;
	}
};
