export module lab:constants;

import :core;

#ifdef __STDCPP_FLOAT128_T__
export namespace lab
{
	template<typename T = double>
	struct constants
	{
		static_assert(std::floating_point<T>);

		static constexpr T pi = 3.1415926535897932384626433832795028841971f128;
	};
}
#else
export namespace lab
{
	template<typename T = double>
	struct constants
	{
		static_assert(std::floating_point<T>);

		static constexpr T pi = 3.1415926535897932384626433832795028841971L;
	};
}
#endif