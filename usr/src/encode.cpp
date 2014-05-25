#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <cmath>
#include <limits>
using namespace std;

double dist2(const array<double, 3>& p0, const array<double, 3>& p1)
{
	const auto d0 = p0[0] - p1[0];
	const auto d1 = p0[1] - p1[1];
	const auto d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

int main(int argc, char* argv[])
{
	while (true)
	{
		vector<array<double, 3>> atoms;
		atoms.reserve(41);
		for (string line, record; getline(cin, line) && (record = line.substr(0, 6)) != "TORSDO";)
		{
			if (record != "ATOM  " && record != "HETATM") continue;
			if (line[77] == 'H' && (line[78] == ' ' || line[78] == 'D')) continue;
			array<double, 3> a = { stof(line.substr(30, 8)), stof(line.substr(38, 8)), stof(line.substr(46, 8)) };
			atoms.push_back(move(a));
		}
		if (atoms.empty()) break;
		const auto n = atoms.size();
		const auto v = 1.0 / n;
		array<double, 3> ctd{};
		array<double, 3> cst{};
		array<double, 3> fct{};
		array<double, 3> ftf{};
		for (size_t i = 0; i < 3; ++i)
		{
			for (const auto& a : atoms)
			{
				ctd[i] += a[i];
			}
			ctd[i] *= v;
		}
		double cst_dist = numeric_limits<double>::max();
		double fct_dist = numeric_limits<double>::lowest();
		double ftf_dist = numeric_limits<double>::lowest();
		for (const auto& a : atoms)
		{
			const auto this_dist = dist2(a, ctd);
			if (this_dist < cst_dist)
			{
				cst = a;
				cst_dist = this_dist;
			}
			if (this_dist > fct_dist)
			{
				fct = a;
				fct_dist = this_dist;
			}
		}
		for (const auto& a : atoms)
		{
			const auto this_dist = dist2(a, fct);
			if (this_dist > ftf_dist)
			{
				ftf = a;
				ftf_dist = this_dist;
			}
		}
		for (const auto& rpt : { ctd, cst, fct, ftf })
		{
			vector<double> dists(n);
			for (size_t i = 0; i < n; ++i)
			{
				dists[i] = sqrt(dist2(atoms[i], rpt));
			}
			array<double, 3> m{};
			for (size_t i = 0; i < n; ++i)
			{
				const auto d = dists[i];
				m[0] += d;
			}
			m[0] *= v;
			for (size_t i = 0; i < n; ++i)
			{
				const auto d = dists[i] - m[0];
				m[1] += d * d;
			}
			m[1] = sqrt(m[1] * v);
			for (size_t i = 0; i < n; ++i)
			{
				const auto d = dists[i] - m[0];
				m[2] += d * d * d;
			}
			m[2] = cbrt(m[2] * v);
			cout.write(reinterpret_cast<char*>(m.data()), sizeof(m));
		}
	}
}
