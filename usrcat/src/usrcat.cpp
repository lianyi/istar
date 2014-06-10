#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <thread>
#include <immintrin.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <boost/tokenizer.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <mongo/client/dbclient.h>
#include <Poco/Net/MailMessage.h>
#include <Poco/Net/MailRecipient.h>
#include <Poco/Net/SMTPClientSession.h>
using namespace std;
using namespace std::chrono;
using namespace OpenBabel;
using namespace boost::filesystem;
using namespace boost::iostreams;
using namespace boost::gregorian;
using namespace boost::posix_time;
using namespace mongo;
using namespace bson;
using namespace Poco::Net;

inline static string now()
{
	return to_simple_string(second_clock::local_time()) + " ";
}

template <typename T>
void read(vector<T>& v, const string f)
{
	ifstream ifs(f, ios::binary);
	ifs.seekg(0, ios::end);
	const size_t num_bytes = ifs.tellg();
	v.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(v.data()), num_bytes);
}

double dist2(const array<double, 3>& p0, const array<double, 3>& p1)
{
	const auto d0 = p0[0] - p1[0];
	const auto d1 = p0[1] - p1[1];
	const auto d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

int main(int argc, char* argv[])
{
	// Check the required number of command line arguments.
	if (argc != 5)
	{
		cout << "usr host user pwd jobs_path" << endl;
		return 0;
	}

	// Fetch command line arguments.
	const auto host = argv[1];
	const auto user = argv[2];
	const auto pwd = argv[3];
	const path jobs_path = argv[4];

	DBClientConnection conn;
	{
		// Connect to host and authenticate user.
		cout << now() << "Connecting to " << host << " and authenticating " << user << endl;
		string errmsg;
		if ((!conn.connect(host, errmsg)) || (!conn.auth("istar", user, pwd, errmsg)))
		{
			cerr << now() << errmsg << endl;
			return 1;
		}
	}
	const auto collection = "istar.usr";
	const auto m256s = _mm256_set1_pd(-0. ); // -0.  = 1 << 63
	const auto qn = 60;
	const auto qv = 1.0 / qn;
	const auto epoch = date(1970, 1, 1);
	const array<string, 5> SmartsPatterns =
	{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	};

	// Read the feature bin file.
	vector<array<double, qn>> features;
	read(features, "16_usrcat.bin");
	const size_t n = features.size();

	// Read the header bin file.
	vector<size_t> headers;
	read(headers, "16_hdr.bin");
	assert(n == headers.size());

	// Search the features for records similar to the query.
	vector<double> scores(n);
	vector<size_t> scase(n);
	array<array<double, 4>, 1> aw;
	auto a = aw.front();
	string line;
	ifstream ligands("16_lig.pdbqt");
	while (true)
	{
		// Fetch jobs.
		auto cursor = conn.query(collection, QUERY("done" << BSON("$exists" << false)).sort("submitted"), 100); // Each batch processes 100 jobs.
		while (cursor->more())
		{
			const auto job = cursor->next();
			const auto _id = job["_id"].OID();
			cout << now() << "Executing job " << _id.str() << endl;

			// Obtain job properties.
			const auto ligand = job["ligand"].str(); // If .String() were used, exceptions would be thrown when the ligand property is not a string.
			const auto format = job["format"].String();

			// Split the ligand string into lines.
			vector<string> lines;
			lines.reserve(200);
			stringstream ss(ligand);
			for (string line; getline(ss, line); lines.push_back(line));

			// Parse the ligand lines.
			vector<array<double, 3>> atoms;
			atoms.reserve(80);
			bool invalid = false;
			try
			{
				if (format == "mol2")
				{
					const auto atomCount = stoul(lines[2].substr(0, 5));
					for (auto i = 0; i < atomCount; ++i)
					{
						const auto& line = lines[7 + i];
						atoms.push_back({ stod(line.substr(16, 10)), stod(line.substr(26, 10)), stod(line.substr(36, 10)) });
					}
				}
				else if (format == "sdf")
				{
					const auto atomCount = stoul(lines[3].substr(0, 3));
					for (auto i = 0; i < atomCount; ++i)
					{
						const auto& line = lines[4 + i];
						atoms.push_back({ stod(line.substr(0, 10)), stod(line.substr(10, 10)), stod(line.substr(20, 10)) });
					}
				}
				else if (format == "xyz")
				{
					boost::char_separator<char> sep(" ");
					const auto atomCount = stoul(lines[0].substr(0, 3));
					for (auto i = 0; i < atomCount; ++i)
					{
						const auto& line = lines[2+i];
						const boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
						auto ti = tokens.begin();
						atoms.push_back({ stod(*++ti), stod(*++ti), stod(*++ti) });
					}
				}
				else if (format == "pdb")
				{
					for (const auto& line : lines)
					{
						const auto record = line.substr(0, 6);
						if (record == "ATOM  " || record == "HETATM")
						{
							atoms.push_back({ stod(line.substr(30, 8)), stod(line.substr(38, 8)), stod(line.substr(46, 8)) });
						}
						else if (record == "ENDMDL") break;
					}
				}
				else if (format == "pdbqt")
				{
					for (const auto& line : lines)
					{
						const auto record = line.substr(0, 6);
						if (record == "ATOM  " || record == "HETATM")
						{
							atoms.push_back({ stod(line.substr(30, 8)), stod(line.substr(38, 8)), stod(line.substr(46, 8)) });
						}
						else if (record == "TORSDO") break;
					}
				}
			}
			catch (const exception& e)
			{
				invalid = true;
			}
			if (invalid || atoms.empty())
			{
				continue;
			}

			OBConversion obConversion;
			obConversion.SetInFormat(format.c_str());
			OBMol obMol;
			obConversion.Read(&obMol, &ss);
			array<vector<int>, 5> subsets;
			for (size_t k = 0; k < 5; ++k)
			{
				auto& subset = subsets[k];
				subset.reserve(atoms.size());
				OBSmartsPattern smarts;
				smarts.Init(SmartsPatterns[k]);
				smarts.Match(obMol);
				for (const auto& map : smarts.GetMapList())
				{
					subset.push_back(map.front() - 1);
				}
			}
			const auto& subset0 = subsets.front();
			const auto n = subset0.size();
			const auto v = 1.0 / n;
			array<double, 3> ctd{};
			array<double, 3> cst{};
			array<double, 3> fct{};
			array<double, 3> ftf{};
			for (size_t k = 0; k < 3; ++k)
			{
				for (const auto i : subset0)
				{
					const auto& a = atoms[i];
					ctd[k] += a[k];
				}
				ctd[k] *= v;
			}
			double cst_dist = numeric_limits<double>::max();
			double fct_dist = numeric_limits<double>::lowest();
			double ftf_dist = numeric_limits<double>::lowest();
			for (const auto i : subset0)
			{
				const auto& a = atoms[i];
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
			for (const auto i : subset0)
			{
				const auto& a = atoms[i];
				const auto this_dist = dist2(a, fct);
				if (this_dist > ftf_dist)
				{
					ftf = a;
					ftf_dist = this_dist;
				}
			}
			array<array<double, qn>, 1> qw;
			auto q = qw.front();
			size_t qo = 0;
			for (const auto& subset : subsets)
			{
				const auto n = subset.size();
				const auto v = 1.0 / n;
				for (const auto& rpt : { ctd, cst, fct, ftf })
				{
					vector<double> dists(n);
					for (size_t i = 0; i < n; ++i)
					{
						dists[i] = sqrt(dist2(atoms[subset[i]], rpt));
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
					for (const auto e : m)
					{
						q[qo++] = e;
					}
				}
			}

			// Compute USRCAT scores.
			for (size_t k = 0; k < n; ++k)
			{
				const auto& l = features[k];
				double s = 0;
				#pragma unroll
				for (size_t i = 0; i < qn; i += 4)
				{
					const auto m256a = _mm256_andnot_pd(m256s, _mm256_sub_pd(_mm256_load_pd(&q[i]), _mm256_load_pd(&l[i])));
					_mm256_stream_pd(a.data(), _mm256_hadd_pd(m256a, m256a));
					s += a[0] + a[2];
				}
				scores[k] = 1 / (1 + s * qv);
			}

			// Sort the scores.
			iota(scase.begin(), scase.end(), 0);
			sort(scase.begin(), scase.end(), [&scores](const size_t val1, const size_t val2)
			{
				return scores[val1] > scores[val2];
			});

			// Write results.
			const auto job_path = jobs_path / _id.str();
			create_directory(job_path);
			filtering_ostream ligands_pdbqt_gz;
			ligands_pdbqt_gz.push(gzip_compressor());
			ligands_pdbqt_gz.push(file_sink((job_path / "ligands.pdbqt.gz").string()));
			ligands_pdbqt_gz.setf(ios::fixed, ios::floatfield);
			ligands_pdbqt_gz << setprecision(8);
			for (size_t k = 0; k < 1000; ++k)
			{
				const size_t c = scase[k];
				ligands.seekg(headers[c]);
				for (size_t i = 0; i < 3 && getline(ligands, line); ++i)
				{
					ligands_pdbqt_gz << line << '\n';
				}
				ligands_pdbqt_gz << "REMARK     USRCAT SCORE: " << setw(10) << scores[c] << '\n';
				while (getline(ligands, line))
				{
					ligands_pdbqt_gz << line << '\n';
					if (line.substr(0, 6) == "TORSDO") break;
				}
			}

			// Update progress.
			const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
			conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));
			const auto err = conn.getLastError();
			if (!err.empty())
			{
				cerr << now() << err << endl;
			}

			// Send completion notification email.
			const auto email = job["email"].String();
			cout << now() << "Sending a completion notification email to " << email << endl;
			MailMessage message;
			message.setSender("istar <noreply@cse.cuhk.edu.hk>");
			message.setSubject("Your usrcat job has completed");
			message.setContent("Your usrcat job submitted on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC with description \"" + job["description"].String() + "\" was done on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC. View result at http://istar.cse.cuhk.edu.hk/usrcat/iview/?" + _id.str());
			message.addRecipient(MailRecipient(MailRecipient::PRIMARY_RECIPIENT, email));
			SMTPClientSession session("137.189.91.190");
			session.login();
			session.sendMessage(message);
			session.close();
		}

		// Sleep for a while.
		this_thread::sleep_for(std::chrono::seconds(10));
	}
}
