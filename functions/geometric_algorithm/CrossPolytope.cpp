/*
 * CrossPolytope
 * Code implementing the "geometric algorithm" in https://arxiv.org/abs/2103.04470
 * This code uses part of Ripser, specifically loading of files and ouptut of 
 * results.
 * 
 Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

 MIT License

 Copyright (c) 2015–2019 Ulrich Bauer

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or
 upgrades to the features, functionality or performance of the source code
 ("Enhancements") to anyone; however, if you choose to make your Enhancements
 available either publicly, or directly to the author of this software, without
 imposing a separate written license agreement for such Enhancements, then you
 hereby grant the following license: a non-exclusive, royalty-free perpetual
 license to install, use, modify, prepare derivative works, incorporate into
 other computer software, distribute, and sublicense such enhancements or
 derivative works thereof, in binary and source code form.

*/

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>

typedef float value_t;
typedef int64_t index_t;

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> struct compressed_distance_matrix {
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	compressed_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	template <typename DistanceMatrix>
	compressed_distance_matrix(const DistanceMatrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();

		for (size_t i = 1; i < size(); ++i)
			for (size_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
	}

	value_t operator()(const index_t i, const index_t j) const;
	size_t size() const { return rows.size(); }
	void init_rows();
};

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;

template <> void compressed_lower_distance_matrix::init_rows() {
	value_t* pointer = &distances[0];
	for (size_t i = 1; i < size(); ++i) {
		rows[i] = pointer;
		pointer += i;
	}
}

template <> void compressed_upper_distance_matrix::init_rows() {
	value_t* pointer = &distances[0] - 1;
	for (size_t i = 0; i < size() - 1; ++i) {
		rows[i] = pointer;
		pointer += size() - i - 2;
	}
}

template <>
value_t compressed_lower_distance_matrix::operator()(const index_t i, const index_t j) const {
	return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

template <>
value_t compressed_upper_distance_matrix::operator()(const index_t i, const index_t j) const {
	return i == j ? 0 : i > j ? rows[j][i] : rows[i][j];
}

struct euclidean_distance_matrix {
	std::vector<std::vector<value_t>> points;

	euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points)
	    : points(std::move(_points)) {
		for (auto p : points) { assert(p.size() == points.front().size()); }
	}

	value_t operator()(const index_t i, const index_t j) const {
		assert(i < points.size());
		assert(j < points.size());
		return std::sqrt(std::inner_product(
		    points[i].begin(), points[i].end(), points[j].begin(), value_t(), std::plus<value_t>(),
		    [](value_t u, value_t v) { return (u - v) * (u - v); }));
	}

	size_t size() const { return points.size(); }
};

template <typename DistanceMatrix> class CrossPolytope {
	const DistanceMatrix dist;

public:
	CrossPolytope(DistanceMatrix&& _dist)
	    : dist(std::move(_dist)) {}
	
	void compute_barcodes() {
		std::cout << "persistence intervals in dim " << (dist.size()-1)/2 << ":" << std::endl;
		
		value_t tb = 0;
		value_t td = std::numeric_limits<value_t>::infinity();
		
		for (int i=0; i<dist.size(); i++) {
			value_t tb_row = 0;
			value_t td_row = dist(i,0);
			
			for (int j=1; j<dist.size(); j++) {
				value_t d = dist(i,j);
				
				if (td_row < d) {
					tb_row = td_row;
					td_row = d;
				} else if (tb_row < d) {
					tb_row = d;
				}
			}
			
			tb = std::max(tb, tb_row);
			td = std::min(td, td_row);
		}
		
		if (tb < td) {
			std::cout << " [" << tb << "," << td << ")" << std::endl;
		}
	}
};

enum file_format {
	LOWER_DISTANCE_MATRIX,
	UPPER_DISTANCE_MATRIX,
	DISTANCE_MATRIX,
	POINT_CLOUD,
	BINARY
};

static const uint16_t endian_check(0xff00);
static const bool is_big_endian = *reinterpret_cast<const uint8_t*>(&endian_check);

template <typename T> T read(std::istream& input_stream) {
	T result;
	char* p = reinterpret_cast<char*>(&result);
	if (input_stream.read(p, sizeof(T)).gcount() != sizeof(T)) return T();
	if (is_big_endian) std::reverse(p, p + sizeof(T));
	return result;
}

compressed_lower_distance_matrix read_point_cloud(std::istream& input_stream) {
	std::vector<std::vector<value_t>> points;

	std::string line;
	value_t value;
	while (std::getline(input_stream, line)) {
		std::vector<value_t> point;
		std::istringstream s(line);
		while (s >> value) {
			point.push_back(value);
			s.ignore();
		}
		if (!point.empty()) points.push_back(point);
		assert(point.size() == points.front().size());
	}

	euclidean_distance_matrix eucl_dist(std::move(points));
	index_t n = eucl_dist.size();
	std::cout << "point cloud with " << n << " points in dimension "
	          << eucl_dist.points.front().size() << std::endl;

	std::vector<value_t> distances;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_upper_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
}

compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;

	std::string line;
	value_t value;
	for (int i = 0; std::getline(input_stream, line); ++i) {
		std::istringstream s(line);
		for (int j = 0; j < i && s >> value; ++j) {
			distances.push_back(value);
			s.ignore();
		}
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_binary(std::istream& input_stream) {
	std::vector<value_t> distances;
	while (!input_stream.eof()) distances.push_back(read<value_t>(input_stream));
	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream, const file_format format) {
	switch (format) {
	case LOWER_DISTANCE_MATRIX:
		return read_lower_distance_matrix(input_stream);
	case UPPER_DISTANCE_MATRIX:
		return read_upper_distance_matrix(input_stream);
	case DISTANCE_MATRIX:
		return read_distance_matrix(input_stream);
	case POINT_CLOUD:
		return read_point_cloud(input_stream);
	default:
		return read_binary(input_stream);
	}
}

void print_usage_and_exit(int exit_code) {
	std::cerr
	    << "Usage: "
	    << "CrossPolytope "
	    << "[options] [filename]" << std::endl
	    << std::endl
	    << "Options:" << std::endl
	    << std::endl
	    << "  --help           print this screen" << std::endl
	    << "  --format         use the specified file format for the input. Options are:"
	    << std::endl
	    << "                     lower-distance (lower triangular distance matrix; default)"
	    << std::endl
	    << "                     upper-distance (upper triangular distance matrix)" << std::endl
	    << "                     distance       (full distance matrix)" << std::endl
	    << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
	    << std::endl
	    << "                     binary         (lower triangular distance matrix in binary format)"
	    << std::endl
	    << std::endl;
	exit(exit_code);
}

int main(int argc, char** argv) {
	const char* filename = nullptr;

	file_format format = DISTANCE_MATRIX;

	for (index_t i = 1; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter.rfind("lower", 0) == 0)
				format = LOWER_DISTANCE_MATRIX;
			else if (parameter.rfind("upper", 0) == 0)
				format = UPPER_DISTANCE_MATRIX;
			else if (parameter.rfind("dist", 0) == 0)
				format = DISTANCE_MATRIX;
			else if (parameter.rfind("point", 0) == 0)
				format = POINT_CLOUD;
			else if (parameter == "binary")
				format = BINARY;
			else
				print_usage_and_exit(-1);
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	} else {
		compressed_lower_distance_matrix dist =
		    read_file(filename ? file_stream : std::cin, format);
		
		std::cout << "distance matrix with " << dist.size() << " points" << std::endl;
		
		CrossPolytope<compressed_lower_distance_matrix>(std::move(dist))
			    .compute_barcodes();

		exit(0);
	}
}
