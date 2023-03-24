module;

#include <cassert>

export module lab:sample;

import :core;
import :estimate;

export namespace lab
{
	namespace _detail
	{
		template<typename ValueType>
		struct analysis_result_t
		{
			using value_type = ValueType;
		private:
			size_t _size;
			value_type _variance, _weight_sum, _weight2_sum;
			estimate_t<value_type> _mean;
		public:
			analysis_result_t(size_t size, estimate_t<value_type> mean, value_type variance, value_type weight_sum, value_type weight2_sum)
				: _size(size), _mean(mean), _variance(variance), _weight_sum(weight_sum), _weight2_sum(weight2_sum) {}

			analysis_result_t(size_t size, estimate_t<value_type> mean, value_type variance)
				: analysis_result_t(size, mean, variance, value_type(size), value_type(1) / size) {}

			analysis_result_t(analysis_result_t const&) = default;


			size_t size() const  { return _size; }

			estimate_t<value_type> mean() const
			{
				assert(size() != 0);
				return _mean;
			}
			value_type variance() const
			{
				assert(size() > 1);
				return _variance;
			}
			value_type stddev() const { return std::sqrt(variance()); }

			value_type weight_sum() const { return _weight_sum; }
			value_type weight2_sum() const { return _weight2_sum; }
		};

		template<typename ValueType>
		struct pair_analysis_result_t
		{
			using value_type = ValueType;
		private:
			size_t _size;
			value_type _x_mean, _x_variance, _y_mean, _y_variance, _covariance, _weight_sum, _weight2_sum;
		public:
			pair_analysis_result_t(size_t size, value_type x_mean, value_type x_variance, value_type y_mean, value_type y_variance, value_type covariance, value_type weight_sum, value_type  weight2_sum)
				: _size(size), _x_mean(x_mean), _x_variance(x_variance), _y_mean(y_mean), _y_variance(y_variance), _covariance(covariance), _weight_sum(weight_sum), _weight2_sum(weight2_sum) {}

			pair_analysis_result_t(size_t size, value_type x_mean, value_type x_variance, value_type y_mean, value_type y_variance, value_type covariance)
				: pair_analysis_result_t(size, x_mean, x_variance, y_mean, y_variance, covariance, value_type(size), value_type(1) / size) {}

			pair_analysis_result_t(pair_analysis_result_t const&) = default;


			size_t size() const { return _size; }

			value_type x_mean() const
			{
				assert(size() != 0);
				return _x_mean;
			}
			value_type x_variance() const
			{
				assert(size() > 1);
				return _x_variance;
			}
			value_type x_stddev() const { return std::sqrt(x_variance()); }
	
			value_type y_mean() const
			{
				assert(size() != 0);
				return _y_mean;
			}
			value_type y_variance() const
			{
				assert(size() > 1);
				return _y_variance;
			}
			value_type y_stddev() const { return std::sqrt(y_variance()); }
			
			value_type covariance() const { return _covariance; }
			value_type weight_sum() const { return _weight_sum; }
			value_type weight2_sum() const { return _weight2_sum; }
		};
	} // namespace _detail

	template<typename Sample>
	auto analyze_sample(Sample&& sample)
	{
		static_assert(stdr::range<Sample>);

		using range_value_t = stdr::range_value_t<Sample>;

		if constexpr (std::floating_point<range_value_t>)
		{
			using value_type = range_value_t;

			size_t size = 0;
			value_type mean = 0, variance = 0;

			for (auto x : sample)
			{
				++size;

				value_type
					delta = x - mean,
					r_delta = delta / size;

				mean += r_delta;
				variance += delta * r_delta * (size - 1);
			}
			variance /= (size - 1);
			return _detail::analysis_result_t(size, estimate_t(mean, variance / size), variance);
		}
		else if constexpr (is_estimate<range_value_t>)
		{
			using value_type = range_value_t::value_type;

			size_t size = 0;
			value_type mean = 0, variance = 0, w_sum = 0, w2_sum = 0;

			for (auto estimate : sample)
			{
				++size;
				
				value_type
					x = estimate.value(),
					w = 1 / estimate.variance(),
					delta = x - mean;
				
				w_sum += w;
				w2_sum += w * w;
				mean += (w / w_sum) * delta;
				variance += w * delta * (x - mean);
			}
			variance /= (w_sum - w2_sum / w_sum);
			return _detail::analysis_result_t(size, estimate_t(mean, 1 / w_sum), variance, w_sum, w2_sum);
		}
		else if constexpr (std::tuple_size_v<range_value_t> == 2)
		{
			using first_type = std::tuple_element_t<0, range_value_t>;
			using second_type = std::tuple_element_t<1, range_value_t>;

			if constexpr (std::floating_point<first_type> && std::floating_point<second_type>)
			{
				using value_type = std::common_type_t<first_type, second_type>;

				size_t size = 0;
				value_type x_mean = 0, x_variance = 0, y_mean = 0, y_variance = 0, covariance = 0;
				for (auto [x, y] : sample)
				{
					++size;

					value_type
						delta_x = x - x_mean,
						r_delta_x = delta_x / size,
						delta_y = y - y_mean,
						r_delta_y = delta_y / size;

					x_mean += r_delta_x;
					x_variance += delta_x * r_delta_x * (size - 1);

					y_mean += r_delta_y;
					y_variance += delta_y * r_delta_y * (size - 1);

					covariance += delta_x * (y - y_mean);
				}
				value_type factor = value_type(1) / (size - 1);
				x_variance *= factor;
				y_variance *= factor;
				covariance *= factor;

				return _detail::pair_analysis_result_t(size, estimate_t(x_mean, x_variance / size), x_variance, estimate_t(y_mean, y_variance / size), y_variance, covariance);
			}
			else if constexpr (is_estimate<first_type> && is_estimate<second_type>)
			{
				using value_type = std::common_type_t<first_type::value_type, second_type::value_type>;

				size_t size = 0;
				value_type x_mean = 0, x_variance = 0, y_mean = 0, y_variance = 0, covariance = 0, w_sum = 0, w2_sum = 0;
				for (auto [x_estimate, y_estimate] : sample)
				{
					++size;

					value_type
						x = x_estimate.value(),
						y = y_estimate.value(),
						w = 1 / y_estimate.variance(),
						delta_x = x - x_mean,
						delta_y = y - y_mean;
					
					w_sum += w;
					w2_sum += w * w;

					x_mean += (w / w_sum) * delta_x;
					x_variance += w * delta_x * (x - x_mean);

					y_mean += (w / w_sum) * delta_y;
					y_variance += w * delta_y * (y - y_mean);

					covariance += w * delta_x * (y - y_mean);
				}
				value_type factor = 1 / (w_sum - w2_sum / w_sum);
				x_variance *= factor;
				y_variance *= factor;
				covariance *= factor;

				return _detail::pair_analysis_result_t(size, x_mean, x_variance, y_mean, y_variance, covariance, w_sum, w2_sum);
			}
			else
				static_assert(false);
		}
		else
			static_assert(false);
	}
}