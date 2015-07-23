#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
using namespace std;

template <typename T>
inline vector<T> read(const string path)
{
	vector<T> buf;
	std::ifstream ifs(path, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cout << "score sample.bin output.bin" << endl;
		return 1;
	}
	const size_t num_usrs = 2;
	constexpr array<size_t, num_usrs> qn{{ 12, 60 }};
	constexpr array<double, num_usrs> qv{{ 1.0 / qn[0], 1.0 / qn[1] }};
	const auto qs = read<array<double, qn.back()>>(argv[1]);
	const auto ls = read<array<double, qn.back()>>(argv[2]);
	const auto q = qs.front();
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(8);
	for (const auto& l : ls)
	{
		double s = 0;
		#pragma unroll
		for (size_t i = 0, u = 0; u < num_usrs; ++u)
		{
			#pragma unroll
			for (; i < qn[u]; i += 4)
			{
				#pragma unroll
				for (size_t o = i; o < i + 4; ++o)
				{
					s += fabs(q[o] - l[o]);
				}
			}
			s = 1 / (1 + s * qv[u]);
			cout << s;
			if (u) cout << endl;
			else cout << '\t';
		}
	}
}
