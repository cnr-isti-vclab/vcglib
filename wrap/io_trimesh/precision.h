#ifndef __VCGLIB_PRECISION
#define __VCGLIB_PRECISION

namespace vcg
{
	namespace tri
	{
		namespace io
		{
			template<typename SCALAR>
			class Precision
			{
				static int decimalPrecision; // the desired decimal precision
			public:
				static int digits() {return 0;}
				static void setDigits(int digits) { decimalPrecision = digits; }
				static const char* typeName() {return "";}
			};

			template <typename SCALAR>
			int Precision<SCALAR>::decimalPrecision = 0;

			// Precision specializations with reasonable defaults

			// Precision<float> specializations
			template<>
			int Precision<float>::digits() { return decimalPrecision <= 0 ? 7 : decimalPrecision; }

			template<>
			const char* Precision<float>::typeName() { return "float"; }

			// Precision<double> specializations
			template<>
			int Precision<double>::digits() { return decimalPrecision <= 0 ? 16 : decimalPrecision; }

			template<>
			const char* Precision<double>::typeName() { return "double"; }

		}
	}
}


#endif
