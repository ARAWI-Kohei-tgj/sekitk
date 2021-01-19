module sekitk.qvm.algorithm.decomposition.lu_band3;

enum string LU_BAND3_IMPL= q{
	/********************************************
	 *
	 ********************************************/
	Result!Method decompose(DecompScheme Method, ArgType)(in ArgType a)
	if(is(Argtype == T[arrayLength!(Size, Size, Shape)])){ // || is(ArgType == T[3u][Size])
		T[arrayLength!(Size, Size, MatrixType.band3)] num= void;

/+
		static if(is(ArgType == T[3u][Size])){
			// TODO:
		}
+/
		num[0]= a[0];
		TypeOfSize idxDiagBef= 0u;
		auto idxSetU= IndexSetUpperTri!(Size, MatrixType.band3, MatOdr)();
		auto idxSetDiag= IndexSetDiag!(Size, Size, MatrixType.band3, MatOdr)();
		auto idxSetL= IndexSetLowerTri!(Size, MatrixType.band3, MatOdr)();
		foreach(idxL, idxU, idxDiagCurr;
						zip(idxSetL, idxSetU, idxSetDiag.dropOne)){
			num[idxU]= a[idxU];
			num[idxL]= a[idxL]/a[idxDiagBef];
			num[idxDiagCurr]= a[idxDiagBef]-num[idxL]*num[idxU];
			idxDiagBef= idxDiagCurr;
		}

		return typeof(return)(num);
	}

	/******************************************
	 * 
	 ******************************************/
	struct Result(DecompScheme Method){
		@safe pure nothrow @nogc:
		// copy constructor
		this(ref return scope inout typeof(this) other){}

		const{
			bool isInvertible(){return _isInvertible;}

			/************************
			 * Lower triangular  matrix $(I L)
			 ************************/
			T[arrayLength!(Size, Size, MatrixType.lowerTri)]
			matrix(DecomposedMat Mat: DecomposedMat.lowerTri)(){
				typeof(return) result= void;
				{
					auto idxSet= IndexSetDiag!(Size, Size, MatrixType.lowerTri, MatOdr)();
					foreach(idx; idxSet) result[idx]= VALUE_ONE;
				}
				auto idxSetSrc= IndexSetStrictTriL!(Size, MatrixType.band3)();
			}

			/************************
			 * Upper triangular matrix $(I U)
			 ************************/
			TypeOfSize[Size] matrix(DecomposedMat Mat: DecomposedMat.permutationLeft)(){
				import std.array: staticArray;
				import std.range: iota;

				return iota!TypeOfSize(Size).staticArray!Size;
			}

			/************************
			 * Permutation matrix $(I P)
			 ************************/
			TypeOfSize[Size] matrix(DecomposedMat Mat: DecomposedMat.permutationLeft)(){
				import std.array: staticArray;
				import std.range: iota;

				return iota!TypeOfSize(Size).staticArray!Size;
			}
		}

	private:
		T[arrayLength!(Size, Size, MatrixType.band3)] _resultDecomp;
		bool _isInvertible;
	}
};
