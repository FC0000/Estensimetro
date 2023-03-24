export module lab:root;

import :core;

import <TAxis.h>;
import <TCanvas.h>;
import <TGraphErrors.h>;
import <TLine.h>;
import <TF1.h>;
import <TH1F.h>;

export namespace lab
{
	template<typename Range>
	void plot_linear_regression(
		std::string_view title,
		std::string_view x_label,
		std::string_view y_label,
		Range&& data,
		auto const& regression_data,
		stdf::path const& path,
		size_t height = 4096, size_t width = 2160)
	{
		static_assert(stdr::range<Range>);

		std::vector<double> x, y, ex, ey;
		double x_min = std::numeric_limits<double>::infinity(), x_max = -x_min, y_min = x_min, y_max = x_max;
		if constexpr (stdr::sized_range<Range>)
		{
			size_t size = stdr::size(data);
			x.reserve(size);
			y.reserve(size);
			ex.reserve(size);
			ey.reserve(size);
		}
		for (auto [first, second] : data)
		{
			if (first.value() < x_min) x_min = first.value();
			if (x_max < first.value()) x_max = first.value();

			if (second.value() < y_min) y_min = second.value();
			if (y_max < second.value()) y_max = second.value();

			x.push_back(first.value());
			y.push_back(second.value());
			ex.push_back(first.stddev());
			ey.push_back(second.stddev());
		}

		TCanvas canvas;
		canvas.SetCanvasSize(UInt_t(height), UInt_t(width));
		double x_margin = (x_max - x_min) / 24, y_margin = (y_max - y_min) / 24;
		auto frame = canvas.DrawFrame(x_min - x_margin, y_min - y_margin, x_max  + x_margin, y_max + y_margin, title.data());

		frame->GetXaxis()->SetTitle(x_label.data());
		frame->GetYaxis()->SetTitle(y_label.data());

		TGraphErrors points(static_cast<Int_t>(x.size()), x.data(), y.data(), ex.data(), ey.data());
		points.Draw("P");

		double slope = regression_data.slope().value(), intercept = regression_data.intercept().value();
		TLine line(x_min, slope * x_min + intercept, x_max, slope * x_max + intercept);
		line.Draw();
			
		canvas.SaveAs(path.string().c_str());
	}
}
