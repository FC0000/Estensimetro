export module lab:regression;

import :core;
import :estimate;
import :sample;

export namespace lab
{
	namespace _detail
	{
		template<typename ValueType>
		struct regression_result_t
		{
			using value_type = ValueType;
		private:
			estimate_t<value_type> _slope, _intercept;
			pair_analysis_result_t<value_type> _sample;
		public:
			regression_result_t(estimate_t<value_type> slope, estimate_t<value_type> intercept, pair_analysis_result_t<value_type> const& sample)
				: _slope(slope), _intercept(intercept), _sample(sample) {}

			auto slope() const { return _slope; }
			value_type slope_stderr() const { return _slope.stddev(); }

			auto intercept() const { return _intercept; }
			value_type intercept_stderr() const { return _intercept.stddev(); }

			auto correlation_coefficient() const { return _sample.covariance() / (_sample.x_stddev()* _sample.y_stddev()); }

			auto const& sample() const { return _sample; }
		};
	} // namespace _detail

	template<typename Sample>
	auto regression(Sample&& sample)
	{
		auto sample_data = analyze_sample(std::forward<Sample>(sample));
		auto
			w_sum = sample_data.weight_sum(),
			w2_sum = sample_data.weight2_sum(),
			covariance = sample_data.covariance(),
			x_mean = sample_data.x_mean(),
			x_variance = sample_data.x_variance(),
			y_mean = sample_data.y_mean(),
			slope = covariance / x_variance,
			intercept = y_mean - slope * x_mean,
			slope_stderr2 = 1 / ((w_sum - w2_sum / w_sum) * x_variance),
			intercept_stderr2 = 1 / w_sum + slope_stderr2 * x_mean * x_mean;

		return _detail::regression_result_t(estimate_t(slope, slope_stderr2), estimate_t(intercept, intercept_stderr2), sample_data);
	}
}