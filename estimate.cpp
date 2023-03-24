module;

#include <cassert>

export module lab:estimate;

import :core;

export namespace lab
{
	struct from_stddev_t {} inline constexpr from_stddev{};

	template<typename T = double>
	struct estimate_t
	{
		using value_type = T;
		static_assert(std::floating_point<value_type>);
	private:
		value_type _value, _variance;
	public:
		constexpr estimate_t()
			: _value(std::numeric_limits<value_type>::quiet_NaN()), _variance(_value)
		{}

		constexpr estimate_t(value_type value, value_type variance)
			: _value(value), _variance(variance)
		{
			assert(variance >= 0);
		}
		
		constexpr estimate_t(from_stddev_t, value_type value, value_type stddev)
			: estimate_t(value, stddev* stddev)
		{
			assert(stddev >= 0);
		}

		constexpr value_type value() const { return _value; }
		constexpr value_type variance() const { return _variance; }
		constexpr value_type stddev() const { return std::sqrt(variance()); }
	};

	template<typename T>
	inline constexpr bool is_estimate = false;

	template<typename T>
	inline constexpr bool is_estimate<estimate_t<T>> = true;


	template<typename T>
	estimate_t<T> operator+(estimate_t<T> lhs, estimate_t<T> rhs)
	{
		return { lhs.value() + rhs.value(), lhs.variance() + rhs.variance() };
	}

	template<typename T>
	estimate_t<T> operator-(estimate_t<T> lhs, estimate_t<T> rhs)
	{
		return { lhs.value() - rhs.value(), lhs.variance() + rhs.variance() };
	}

	template<typename T>
	estimate_t<T> operator*(estimate_t<T> lhs, estimate_t<T> rhs)
	{
		auto lhs_mean2 = lhs.value() * lhs.value(), rhs_mean2 = rhs.value() * rhs.value();
		return { lhs.value() * rhs.value(), lhs.variance() * rhs.variance() + lhs.variance() * rhs_mean2 + lhs_mean2 * rhs.variance() };
	}

	// ratio defined below


	template<typename T, typename U>
	estimate_t<T> operator+(estimate_t<T> e, U v)
	{
		return { e.value() + v, e.variance() };
	}

	template<typename T, typename U>
	estimate_t<T> operator+(U v, estimate_t<T> e)
	{
		return e + v;
	}

	template<typename T, typename U>
	estimate_t<T> operator-(estimate_t<T> e, U v)
	{
		return { e.value() - v, e.variance() };
	}

	template<typename T, typename U>
	estimate_t<T> operator-(U v, estimate_t<T> e)
	{
		return { v - e.value(), e.variance() };
	}

	template<typename T, typename U>
	estimate_t<T> operator*(estimate_t<T> e, U v)
	{
		return { e.value() * v, e.variance() * v * v };
	}

	template<typename T, typename U>
	estimate_t<T> operator*(U v, estimate_t<T> e)
	{
		return e * v;
	}

	template<typename T, typename U>
	estimate_t<T> operator/(estimate_t<T> e, U v)
	{
		return e * (T(1) / v);
	}

	template<typename T, typename U>
	estimate_t<T> operator/(U v, estimate_t<T> e)
	{
		return estimate_t<T>(v, 0) / e;
	}
}

export namespace std
{
	template<typename T, typename CharT>
	struct formatter<lab::estimate_t<T>, CharT> : formatter<T, CharT>
	{
		template<typename FormatContext>
		constexpr auto format(lab::estimate_t<T> estimate, FormatContext& format_context) const
		{
			auto out = formatter<T, CharT>::format(estimate.value(), format_context);
			out = ranges::copy_n(" +- ", 4, out).out;
			format_context.advance_to(out);
			return formatter<T, CharT>::format(estimate.stddev(), format_context);
		}
	};
}

/* internal */ namespace lab
{
	template<typename T>
	struct _value_type
	{
		using type = T;
	};
	template<typename T> requires is_estimate<T>
	struct _value_type<T>
	{
		using type = T::value_type;
	};
}

export namespace lab
{
	template<typename T, size_t N, typename Range = decltype(stdv::repeat(0))>
	auto estimate(
		auto const& function,
		std::array<T, N> const& arguments,
		Range&& covariance_triangular_matrix = stdv::repeat(0) // stricly triangular if is_estimate<T>
	)
	{
		static_assert(std::floating_point<T> || is_estimate<T>);
		static_assert(stdr::range<Range>);
		
		using value_type = _value_type<T>::type;

		auto [value, derivative] = [&]<size_t... In>(std::index_sequence<In...>)
		{
			if constexpr(std::floating_point<T>)
				return std::pair(function.value_at(arguments[In]...), function.derivative_at(arguments[In]...));
			else
				return std::pair(function.value_at(arguments[In].value()...), function.derivative_at(arguments[In].value()...));
		}(std::make_index_sequence<N>());

		value_type variance = 0, compensation = 0;
		auto compensated_add = [&](value_type new_term)
		{
			value_type t = variance + new_term;
			if (std::abs(variance) >= std::abs(new_term))
				compensation += (variance - t) + new_term;
			else
				compensation += (new_term - t) + variance;
			variance = t;
		};
		auto cov_it = stdr::begin(covariance_triangular_matrix);
		for (size_t i = 0; i != N; ++i)
		{
			if constexpr (std::floating_point<T>)
				compensated_add(derivative[i] * derivative[i] * (*cov_it++));
			else
				compensated_add(derivative[i] * derivative[i] * arguments[i].variance());
			for (size_t j = i; j != N; ++j)
				compensated_add(2 * derivative[i] * derivative[j] * (*cov_it++));
		}
		variance += compensation;
		return estimate_t(value, variance);
	}

	
	template<typename T>
	estimate_t<T> operator/(estimate_t<T> lhs, estimate_t<T> rhs)
	{
		struct
		{
			static T value_at(T l, T r)
			{
				return l / r;
			}
			static auto derivative_at(T l, T r)
			{
				return std::array{1 / r, -l / (r * r) };
			}
		} quotient;
		return lab::estimate(quotient, std::array{lhs, rhs});
	}
}