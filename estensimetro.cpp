
import lab;

namespace stdv = std::views;
namespace stdr = std::ranges;
namespace stdf = std::filesystem;

struct __efn
{
	auto value_at(auto x0, auto d, auto k) const
	{
		using cnst = lab::constants<decltype(x0)>;
		return 4 * x0 / (cnst::pi * d * d * k);
	}
	auto derivative_at(auto x0, auto d, auto k) const
	{
		using cnst = lab::constants<decltype(x0)>;
		auto common = 4 / (cnst::pi * d * d * k);
		return std::array{common, -2 * x0 * common / d, -1 * x0 * common / k};
	}
} inline constexpr e_fn;

template<typename value_type = double>
auto analyze(stdf::path extension_input_path, stdf::path compression_input_path, stdf::path extension_output_path, stdf::path compression_output_path,
	lab::estimate_t<value_type> x0,
	lab::estimate_t<value_type> d,
	std::ostream& output = std::cout)
{
	using estimate_t = lab::estimate_t<value_type>;

	std::vector<estimate_t> extensions, compressions;
	for(auto [path, vec] : {std::tie(extension_input_path, extensions), std::tie(compression_input_path, compressions)})
	{
		std::ifstream input(path);
		if (!input)
			throw std::runtime_error(std::format("Cannot open {} for reading.", path.string()));
		input.exceptions(input.badbit);

		constexpr size_t chunk_size = 5;
		constexpr value_type conv_factor = 1e-6;

		vec = std::vector<estimate_t>(std::from_range,
			stdv::istream<value_type>(input) |
			stdv::chunk(chunk_size) |
			stdv::transform([conv_factor](auto x) {return conv_factor * lab::analyze_sample(x).mean(); }));
	}
	std::print(output,
		"Extension sample (m):\n{:.6f}\n\n"
		"Compression sample (m):\n{:.6f}\n\n",
		extensions,
		compressions);


	std::vector<std::pair<estimate_t, estimate_t>> regression_ext_data, regression_compr_data;

	constexpr value_type force_conversion_factor = 4 * 9.806 / 1000;

	auto forces = stdv::iota(2, 13) | stdv::transform([](int i) {return force_conversion_factor * estimate_t(i * 100, 0 /*!!!*/); });
	auto force_extension_pairs = stdv::zip(forces, extensions);
	auto [init_force, init_length] = force_extension_pairs[0];

	estimate_t extension_k, compression_k;
	std::print(output, "Extension:\n");
	{
		for (auto [force, extension] : force_extension_pairs | stdv::drop(1))
		{
			auto delta_f = force - init_force, delta_x = extension - init_length;
			regression_ext_data.push_back({ delta_f, delta_x });
			std::print(output, "{:.0f} gp\t:\t{:.6f} m\n", delta_f / force_conversion_factor, delta_x);
		}
		auto regression_result = lab::regression(regression_ext_data);
		extension_k = regression_result.slope();
		std::print(output,
			"\nAllungamento: Dx = ({:.8f}) * DF + {:.8f}\n",
			regression_result.slope(), regression_result.intercept());

		//                        vvvvv titolo
		lab::plot_linear_regression("", "\\Delta F (N)", "\\Delta x (m)", regression_ext_data | stdv::transform([](auto x) {return std::pair(x.first, estimate_t(x.second.value(), x.second.variance() * 100)); }), regression_result, extension_output_path);
	}

	std::print(output, "\nCompression:\n");
	auto force_compression_pairs = stdv::zip(forces | stdv::reverse, compressions);
	{
		for (auto [force, compression] : force_compression_pairs | stdv::take(10))
		{
			auto delta_f = force - init_force, delta_x = compression - init_length;
			regression_compr_data.push_back({ delta_f, delta_x });
			std::print(output, "{:.0f} gp\t:\t{:.6f} m\n", delta_f / force_conversion_factor, delta_x);
		}
		auto regression_result = lab::regression(regression_compr_data);
		compression_k = regression_result.slope();
		std::print(output,
			"\nAccorciamento: Dx = ({:.8f}) * DF + {:.8f}\n",
			regression_result.slope(), regression_result.intercept());
		//                        vvvvv titolo
		lab::plot_linear_regression("",  "\\Delta F (N)", "\\Delta x (m)", regression_compr_data | stdv::transform([](auto x) {return std::pair(x.first, estimate_t(x.second.value(), x.second.variance() * 100)); }), regression_result, compression_output_path);
	}
	auto k_average = lab::analyze_sample(std::array{extension_k, compression_k});
	std::print(output, "\nK = {:.8f} m/N\n", k_average.mean());



	std::print(output, "\nMetodo ISO:\n");
	std::vector<value_type> iso_ks;
	std::vector<estimate_t> es;
	for (auto [first, second] : force_extension_pairs | stdv::adjacent<2> | stdv::stride(2))
	{
		estimate_t
			deltax = std::get<1>(second) - std::get<1>(first),
			deltaf = std::get<0>(second) - std::get<0>(first),
			k = deltax / deltaf;
		iso_ks.push_back(k.value());
		es.push_back(lab::estimate(e_fn, std::array{x0, d, k}));

		std::print(output,
			"{:.0f} gp\t:\t{:.6f} m\n"
			"{:.0f} gp\t:\t{:.6f} m\n\n",
			std::get<0>(first) / force_conversion_factor, std::get<1>(first),
			std::get<0>(second) / force_conversion_factor, std::get<1>(second));
	}
	auto iso_ext = lab::analyze_sample(iso_ks);
	std::print(output,
		"Campione (m/N): {:.8f}\n"
		"K_allungamento = {:.8f} m/N\n\n",
		iso_ks,
		iso_ext.mean());

	iso_ks.clear();

	for (auto [first, second] : force_compression_pairs | stdv::adjacent<2> | stdv::stride(2))
	{
		estimate_t
			deltax = std::get<1>(second) - std::get<1>(first),
			deltaf = std::get<0>(second) - std::get<0>(first),
			k = deltax / deltaf;
		iso_ks.push_back(k.value());
		es.push_back(lab::estimate(e_fn, std::array{x0, d, k}));
		std::print(output,
			"{:.0f} gp\t:\t{:.6f} m\n"
			"{:.0f} gp\t:\t{:.6f} m\n\n",
			std::get<0>(first) / force_conversion_factor, std::get<1>(first),
			std::get<0>(second) / force_conversion_factor, std::get<1>(second));
	}
	auto iso_compr = lab::analyze_sample(iso_ks);
	std::print(output,
		"Campione (m/N): {:.8f}\n"
		"K_accorciamento = {:.8f} m/N\n",
		iso_ks,
		iso_compr.mean());

	k_average = lab::analyze_sample(std::array{iso_ext.mean(), iso_compr.mean()});
	std::print(output, "\nK = {:.8f} m/N\nE (con <K>) = {:.4} Pa\nE (con xi) = {:.4} Pa\n\n", k_average.mean(), lab::estimate(e_fn, std::array{x0, d, k_average.mean()}), lab::analyze_sample(es).mean());


	std::print(output,
		"Verifica di Dx=K*DF con K del metodo ISO\n"
	);

	value_type x2 = 0;
	size_t size = 0;
	for (auto [df, dx] : stdv::join(std::array{ regression_ext_data, regression_compr_data }))
	{
		++size;
		value_type term = dx.value() - k_average.mean().value() * df.value();
		x2 += term * term / dx.variance();
	}
	std::print(output, "N = {}\tV = 1\tGDL = {}\tX^2 = {:.2f}\tX_0^2(95%) = {}\n", size, size - 1, x2, 28.87);
	std::print(output, "Coefficiente di correlazione: {}\n", lab::regression(stdv::join(std::array{ regression_ext_data, regression_compr_data })).correlation_coefficient());

	return k_average.mean();
}

template<typename value_type = double>
void analyze_error(stdf::path input_path400, stdf::path input_path1000)
{
	using estimate_t = lab::estimate_t<value_type>;

	value_type x400_ext, x400_compr, x1000_ext, x1000_compr;
	for (auto& [path, x_ext, x_compr] : { std::tie(input_path400, x400_ext, x400_compr), std::tie(input_path1000, x1000_ext, x1000_compr) })
	{
		std::ifstream input(path);
		if (!input)
			throw std::runtime_error(std::format("Cannot open {} for reading.", path.string()));
		input.exceptions(input.badbit);

		constexpr value_type conv_factor = 1e-6;

		std::vector<value_type> ext, compr;

		for (auto r : stdv::istream<value_type>(input) | stdv::chunk(6))
		{
			ext.append_range(r | stdv::take(3));
			compr.append_range(r);
		}

		x_ext = conv_factor * lab::analyze_sample(ext).mean().value();
		x_compr = conv_factor * lab::analyze_sample(compr).mean().value();
	}
	std::print("\nErrore sistematico: delta = {} m/N\n", ((x1000_ext - x400_ext) - (x1000_compr - x400_compr)) / (600 * 4 * 9.806 / 1000));
}

int main()
{
	using value_type = double;
	using estimate_t = lab::estimate_t<value_type>;
	using cnst = lab::constants<value_type>;

	stdf::path base_path = "";

	struct obj
	{
		estimate_t x0, d;
	};

	auto id_to_index = [](int i)
	{
		if (i > 12) return i - 3;
		return i - 2;
	};

	auto objs = std::array<obj, 17>{
		obj{{1000, 2 * 2}, {0.250, 0.005 * 0.005}}, // 2
		obj{{1000, 2 * 2}, {0.500, 0.005 * 0.005}}, // 3
		obj{{950, 2 * 2}, {0.229, 0.229 * 0.229 / 10000}}, // 4
		obj{{950, 2 * 2}, {0.279, 0.279 * 0.279 / 10000}}, // 5
		obj{{950, 2 * 2}, {0.305, 0.305 * 0.305 / 10000}}, // 6
		obj{{950, 2 * 2}, {0.330, 0.330 * 0.330 / 10000}}, // 7
		obj{{950, 2 * 2}, {0.356, 0.356 * 0.356 / 10000}}, // 8
		obj{{950, 2 * 2}, {0.381, 0.381 * 0.381 / 10000}}, // 9
		obj{{950, 2 * 2}, {0.406, 0.406 * 0.406 / 10000}}, // 10
		obj{{950, 2 * 2}, {0.432, 0.432 * 0.432 / 10000}}, // 11
		obj{{900, 2 * 2}, {0.432, 0.432 * 0.432 / 10000}}, // 13
		obj{{800, 2 * 2}, {0.279, 0.279 * 0.279 / 10000}}, // 14
		obj{{700, 2 * 2}, {0.279, 0.279 * 0.279 / 10000}}, // 15
		obj{{600, 2 * 2}, {0.279, 0.279 * 0.279 / 10000}}, // 16
		obj{{500, 2 * 2}, {0.279, 0.279 * 0.279 / 10000}}, // 17
		obj{{400, 2 * 2}, {0.279, 0.279 * 0.279 / 10000}}, // 18
		obj{{300, 2 * 2}, {0.279, 0.279 * 0.279 / 10000}} // 19
	} | stdv::transform([](obj o) {return obj(o.x0 * 0.001, o.d * 0.001); });

	std::array<estimate_t, 17> ks =
	{
		estimate_t{lab::from_stddev, 0.000053970, 0.000000127}, // 2, GIO5
		estimate_t{}, // 3*
		estimate_t{}, // 4*
		estimate_t{lab::from_stddev, 0.00007391884322, 0.000000683181686057084}, // 5 GIO2
		estimate_t{lab::from_stddev, 6.3678E-05, 1.1399E-07}, // 6, GIO11
		estimate_t{lab::from_stddev, 0.00005647, 0.0000001}, // 7, GIO7
		estimate_t{lab::from_stddev, 4.470E-05, 5E-08}, // 8, GIO1
		estimate_t{lab::from_stddev, 0.00004137526288, 0.0000002255753081}, // 9 GIO2
		estimate_t{lab::from_stddev, 0.00003538422429, 0.0000002483229556}, // 10, GIO2
		estimate_t{lab::from_stddev, 3.1678E-5, 1.32E-07}, // 11, GIO9
		estimate_t{}, // 13*
		estimate_t{}, // 14*
		estimate_t{lab::from_stddev, 5.68E-05, 4E-07}, // 15, GIO1
		estimate_t{}, // 16*
		estimate_t{lab::from_stddev, 3.91E-5, 0.01E-5}, // 17, GIO8
		estimate_t{lab::from_stddev, 3.31E-05, 1E-07}, // 18, GIO1
		estimate_t{lab::from_stddev, 0.00002408, 0.0000002} // 19, GIO7
	};

	for (int i : {3, 4, 13, 14, 16})
	{
		estimate_t k = analyze(
			base_path / std::format("{}al.txt", i),
			base_path / std::format("{}ac.txt", i),
			base_path / std::format("{}al.png", i),
			base_path / std::format("{}ac.png", i),
			objs[id_to_index(i)].x0,
			objs[id_to_index(i)].d
			, std::ofstream(base_path / std::format("{}out.txt", i))
		);
		ks[id_to_index(i)] = k;
	}

	analyze_error(base_path / "4_s400.txt", base_path / "4_s1000.txt");
	
	{
		std::print("\nL = 950mm (estensimetri 4~11)\n");
		auto data = std::array{4, 5, 6, 7, 8, 9, 10, 11} | stdv::transform([&](int i) {return std::pair(4.0 / (lab::constants<value_type>::pi * objs[id_to_index(i)].d * objs[id_to_index(i)].d), ks[id_to_index(i)]); });

		auto regression_result = lab::regression(data);
		lab::plot_linear_regression("Lunghezza a riposo costante", "1/S (m^{-2})", "K (mN^{-1})", data, regression_result, base_path / "constant_L.png");

	}

	{
		std::print("D = 0.279mm (estensimetri 5, 14~19)\n");
		auto data = std::array{5, 14, 15, 16, 17, 18, 19} | stdv::transform([&](int i) {return std::pair(objs[id_to_index(i)].x0, ks[id_to_index(i)]); });
		auto regression_result = lab::regression(data);
		lab::plot_linear_regression("Sezione costante", "x_{0} (m)", "K (mN^{-1})", data, regression_result, base_path / "constant_D.png");
	}

	// do not include 3
	auto es = std::array{4, 13, 14, 16, 5, 14, 15, 16, 17, 18, 19} |
		stdv::transform([&](int i)
			{
				auto idx = id_to_index(i);
				auto o = objs[idx];
				return 4 * o.x0 / (cnst::pi * o.d * o.d * ks[idx]);
			});

	int brass = id_to_index(3);
	std::print(
		"Modulo di Young acciaio (ISO): {}\n"
		"Modulo di Young ottone: {}\n",
		lab::analyze_sample(es).mean(),
		4 * objs[brass].x0 / (cnst::pi * objs[brass].d * objs[brass].d * ks[brass])
	);
}
