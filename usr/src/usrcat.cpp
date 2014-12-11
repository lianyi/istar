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
#include <boost/filesystem/fstream.hpp>
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

	// Initialize constants.
	const auto collection = "istar.usr";
	const auto m256s = _mm256_set1_pd(-0. ); // -0.  = 1 << 63
	const auto epoch = date(1970, 1, 1);
	const size_t num_usrs = 2;
	constexpr array<size_t, num_usrs> qn = { 12, 60 };
	constexpr array<double, num_usrs> qv = { 1.0 / qn[0], 1.0 / qn[1] };
	const size_t num_references = 4;
	const size_t num_subsets = 5;
	const array<string, num_subsets> SmartsPatterns =
	{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	};

	// Read the header bin file.
	vector<size_t> headers;
	{
		std::ifstream hdr_bin("16_hdr.bin", ios::binary | ios::ate);
		const size_t num_bytes = hdr_bin.tellg();
		headers.resize(num_bytes / sizeof(size_t));
		hdr_bin.seekg(0);
		hdr_bin.read(reinterpret_cast<char*>(headers.data()), num_bytes);
	}
	const size_t num_ligands = headers.size();

	// Search the features for records similar to the query.
//	vector<array<double, qn.back()>> features;
	array<array<double, qn.back()>, 2> qlw;
	array<array<double, 4>, 1> aw;
	auto q = qlw[0];
	auto l = qlw[1];
	auto a = aw.front();
	array<vector<double>, 2> scores{ vector<double>(num_ligands), vector<double>(num_ligands) };
	vector<size_t> scase(num_ligands);
	string line;
	std::ifstream usrcat_bin("16_usrcat.bin", ios::binary);
	std::ifstream lig_pdbqt("16_lig.pdbqt");
	while (true)
	{
		// Fetch jobs.
		auto cursor = conn.query(collection, QUERY("done" << BSON("$exists" << false)).sort("submitted"), 1); // Each batch processes 1 job.
		while (cursor->more())
		{
			const auto job = cursor->next();
			const auto _id = job["_id"].OID();
			cout << now() << "Executing job " << _id.str() << endl;
			const auto job_path = jobs_path / _id.str();

			// Obtain job properties.
			const auto format = job["format"].String();
			const auto email = job["email"].String();
			const auto ligand_path = job_path / ("ligand." + format);

			// Parse the user-supplied ligand.
			bool invalid = false;
			vector<array<double, 3>> atoms;
			boost::filesystem::ifstream ligand_format(ligand_path);
			try
			{
				if (format == "mol2")
				{
					// Locate @<TRIPOS>MOLECULE so as to parse the number of atoms.
					while (getline(ligand_format, line))
					{
						if (line.substr(0, 17) == "@<TRIPOS>MOLECULE") break;
					}
					getline(ligand_format, line);
					getline(ligand_format, line);
					const auto atomCount = stoul(*boost::tokenizer<boost::char_separator<char>>(line, boost::char_separator<char>(" ")).begin());
					atoms.reserve(atomCount);

					// Locate @<TRIPOS>MOLECULE so as to parse the atoms.
					while (getline(ligand_format, line))
					{
						if (line.substr(0, 13) == "@<TRIPOS>ATOM") break;
					}
					for (auto i = 0; i < atomCount; ++i)
					{
						getline(ligand_format, line);
						atoms.push_back({ stod(line.substr(16, 10)), stod(line.substr(26, 10)), stod(line.substr(36, 10)) });
					}
				}
				else if (format == "sdf") // $$$$ at the end of a molecule.
				{
					// Parse the number of atoms.
					getline(ligand_format, line);
					getline(ligand_format, line);
					getline(ligand_format, line);
					getline(ligand_format, line);
					const auto atomCount = stoul(line.substr(0, 3));
					atoms.reserve(atomCount);

					// Parse the atoms.
					for (auto i = 0; i < atomCount; ++i)
					{
						getline(ligand_format, line);
						atoms.push_back({ stod(line.substr( 0, 10)), stod(line.substr(10, 10)), stod(line.substr(20, 10)) });
					}
				}
				else if (format == "xyz")
				{
					// Parse the number of atoms.
					getline(ligand_format, line);
					const auto atomCount = stoul(line);
					atoms.reserve(atomCount);

					// Parse the atoms.
					const boost::char_separator<char> sep(" ");
					getline(ligand_format, line);
					for (auto i = 0; i < atomCount; ++i)
					{
						getline(ligand_format, line);
						const boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
						auto ti = tokens.begin();
						atoms.push_back({ stod(*++ti), stod(*++ti), stod(*++ti) });
					}
				}
				else if (format == "pdb")
				{
					atoms.reserve(80);
					while (getline(ligand_format, line))
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
					atoms.reserve(80);
					while (getline(ligand_format, line))
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
				// Update progress.
				const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
				conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));
				const auto err = conn.getLastError();
				if (!err.empty())
				{
					cerr << now() << err << endl;
				}

				// Send error notification email.
				cout << now() << "Sending an error notification email to " << email << endl;
				MailMessage message;
				message.setSender("istar <noreply@cse.cuhk.edu.hk>");
				message.setSubject("Your usr job has failed to complete");
				message.setContent("Your usr job submitted on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC with description \"" + job["description"].String() + "\" failed on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC because your provided ligand failed to be parsed.");
				message.addRecipient(MailRecipient(MailRecipient::PRIMARY_RECIPIENT, email));
				SMTPClientSession session("137.189.91.190");
				session.login();
				session.sendMessage(message);
				session.close();
			}

			OBMol obMol;
			OBConversion obConversion;
			obConversion.SetInFormat(format.c_str());
			obConversion.ReadFile(&obMol, ligand_path.string());
			array<vector<int>, num_subsets> subsets;
			for (size_t k = 0; k < num_subsets; ++k)
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
			array<array<double, 3>, num_references> references{};
			auto& ctd = references[0];
			auto& cst = references[1];
			auto& fct = references[2];
			auto& ftf = references[3];
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
			array<vector<double>, num_references> dista;
			for (size_t k = 0; k < num_references; ++k)
			{
				const auto& reference = references[k];
				auto& dists = dista[k];
				dists.resize(atoms.size());
				for (size_t i = 0; i < n; ++i)
				{
					dists[subset0[i]] = sqrt(dist2(atoms[subset0[i]], reference));
				}
			}
			size_t qo = 0;
			for (const auto& subset : subsets)
			{
				const auto n = subset.size();
				for (size_t k = 0; k < num_references; ++k)
				{
					const auto& distp = dista[k];
					vector<double> dists(n);
					for (size_t i = 0; i < n; ++i)
					{
						dists[i] = distp[subset[i]];
					}
					array<double, 3> m{};
					if (n > 2)
					{
						const auto v = 1.0 / n;
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
					}
					else if (n == 2)
					{
						m[0] = 0.5 *     (dists[0] + dists[1]);
						m[1] = 0.5 * fabs(dists[0] - dists[1]);
					}
					else if (n == 1)
					{
						m[0] = dists[0];
					}
					for (const auto e : m)
					{
						q[qo++] = e;
					}
				}
			}
			assert(qo == qn.back());

			// Read features chunk by chunk, and compute USR and USRCAT scores.
			usrcat_bin.seekg(0);
			for (size_t k = 0; k < num_ligands; ++k)
			{
//				const auto& l = features[k];
				usrcat_bin.read(reinterpret_cast<char*>(l.data()), sizeof(l));
				double s = 0;
				size_t i = 0;
				#pragma unroll
				for (size_t u = 0; u < num_usrs; ++u)
				{
					#pragma unroll
					for (; i < qn[u]; i += 4)
					{
						const auto m256a = _mm256_andnot_pd(m256s, _mm256_sub_pd(_mm256_load_pd(&q[i]), _mm256_load_pd(&l[i])));
						_mm256_stream_pd(a.data(), _mm256_hadd_pd(m256a, m256a));
						s += a[0] + a[2];
					}
					scores[u][k] = 1 / (1 + s * qv[u]);
				}
			}

			// Sort the ligands by USRCAT scores.
			const auto& uscores = scores[1];
			iota(scase.begin(), scase.end(), 0);
			sort(scase.begin(), scase.end(), [&uscores](const size_t val1, const size_t val2)
			{
				return uscores[val1] > uscores[val2];
			});

			// Write results.
			filtering_ostream ligands_pdbqt_gz;
			ligands_pdbqt_gz.push(gzip_compressor());
			ligands_pdbqt_gz.push(file_sink((job_path / "ligands.pdbqt.gz").string()));
			ligands_pdbqt_gz.setf(ios::fixed, ios::floatfield);
			ligands_pdbqt_gz << setprecision(8);
			for (size_t k = 0; k < 1000; ++k)
			{
				const size_t c = scase[k];
				lig_pdbqt.seekg(headers[c]);
				for (size_t i = 0; i < 3 && getline(lig_pdbqt, line); ++i)
				{
					ligands_pdbqt_gz << line << '\n';
				}
				ligands_pdbqt_gz << "REMARK 951    USR SCORE: " << setw(10) << scores[0][c] << '\n';
				ligands_pdbqt_gz << "REMARK 952 USRCAT SCORE: " << setw(10) << scores[1][c] << '\n';
				while (getline(lig_pdbqt, line))
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
			cout << now() << "Sending a completion notification email to " << email << endl;
			MailMessage message;
			message.setSender("istar <noreply@cse.cuhk.edu.hk>");
			message.setSubject("Your usr job has completed");
			message.setContent("Your usr job submitted on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC with description \"" + job["description"].String() + "\" was done on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC. View result at http://istar.cse.cuhk.edu.hk/usr/iview/?" + _id.str());
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
