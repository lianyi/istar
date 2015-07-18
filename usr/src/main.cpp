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

inline static string local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " ";
}

template <typename T>
inline vector<T> read(const string path)
{
	cout << local_time() << "Reading " << path << endl;
	vector<T> buf;
	std::ifstream ifs(path, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}

struct zproperty
{
	float mwt, lgp, ads, pds;
	int16_t hbd, hba, psa, chg, nrb;
};

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

	// Connect to host and authenticate user.
	DBClientConnection conn;
	{
		cout << local_time() << "Connecting to " << host << " and authenticating " << user << endl;
		string errmsg;
		if ((!conn.connect(host, errmsg)) || (!conn.auth("istar", user, pwd, errmsg)))
		{
			cerr << local_time() << errmsg << endl;
			return 1;
		}
	}

	// Initialize constants.
	cout << local_time() << "Initializing" << endl;
	const auto collection = "istar.usr";
	const auto epoch = date(1970, 1, 1);
	const size_t num_usrs = 2;
	constexpr array<size_t, num_usrs> qn{{ 12, 60 }};
	constexpr array<double, num_usrs> qv{{ 1.0 / qn[0], 1.0 / qn[1] }};
	const size_t num_references = 4;
	const size_t num_subsets = 5;
	const array<string, num_subsets> SubsetSMARTS
	{{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	}};

	// Initialize variables.
	array<array<double, qn.back()>, 1> qw;
	array<array<double, qn.back()>, 1> lw;
	auto q = qw[0];
	auto l = lw[0];
	string buf;
	buf.resize(5366); // Maximum size of a ligand out of the 23M ligands in PDBQT format.

	// Read ID file.
	cout << local_time() << "Reading 16_zincid.txt" << endl;
	string zincids;
	{
		std::ifstream ifs("16_zincid.txt", ios::binary | ios::ate);
		const size_t num_bytes = ifs.tellg();
		zincids.resize(num_bytes);
		ifs.seekg(0);
		ifs.read(const_cast<char*>(zincids.data()), num_bytes);
	}
	assert(zincids.size() % 9 == 0);
	const size_t num_ligands = zincids.size() / 9;
	array<vector<double>, 2> scores{{ vector<double>(num_ligands, 0), vector<double>(num_ligands, 0) }};
	vector<size_t> scase(num_ligands);

	// Read ZINC property file.
	cout << local_time() << "Reading ZINC property file" << endl;
	vector<zproperty> zproperties(num_ligands);
	{
		std::ifstream ifs("16_zproperty.bin", ios::binary);
		for (auto& p : zproperties)
		{
			ifs.read(reinterpret_cast<char*>(&p), 26); // sizeof(zproperty) == 28
		}
	}

	// Read SMILES file.
	cout << local_time() << "Reading SMILES file" << endl;
	vector<string> smileses(num_ligands);
	{
		std::ifstream ifs("16_smiles.txt");
		for (auto& smiles : smileses)
		{
			getline(ifs, smiles);
		}
	}

	// Read supplier file.
	cout << local_time() << "Reading supplier file" << endl;
	vector<string> suppliers(num_ligands);
	{
		std::ifstream ifs("16_supplier.txt");
		for (auto& supplier : suppliers)
		{
			getline(ifs, supplier);
		}
	}

	// Read header file.
	const auto headers = read<size_t>("16_header.bin");
	assert(headers.size() == num_ligands + 1);

	// Open files for subsequent reading.
	std::ifstream usrcat_bin("16_usrcat.bin");
	std::ifstream ligand_pdbqt("16_ligand.pdbqt");

	// Enter event loop.
	cout << local_time() << "Entering event loop" << endl;
	bool sleeping = false;
	while (true)
	{
		// Fetch an incompleted job in a first-come-first-served manner.
		if (!sleeping) cout << local_time() << "Fetching an incompleted job" << endl;
//		auto cursor = conn.query(collection, QUERY("done" << BSON("$exists" << false)).sort("submitted"), 1); // Each batch processes 1 job.
		BSONObj info;
		conn.runCommand("istar", BSON("findandmodify" << "usr" << "query" << BSON("done" << BSON("$exists" << false) << "started" << BSON("$exists" << false)) << "sort" << BSON("submitted" << 1) << "update" << BSON("$set" << BSON("started" << Date_t(duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count()))) << "fields" << BSON("_id" << 1 << "format" << 1 << "email" << 1 << "submitted" << 1 << "description" << 1)), info); // conn.findAndModify() is available since MongoDB C++ Driver legacy-1.0.0
		const auto value = info["value"];
		if (value.isNull())
		{
			// No incompleted jobs. Sleep for a while.
			if (!sleeping) cout << local_time() << "Sleeping for 10 seconds" << endl;
			sleeping = true;
			this_thread::sleep_for(chrono::seconds(10));
			continue;
		}
		sleeping = false;
		const auto job = value.Obj();

		// Obtain job properties.
		const auto _id = job["_id"].OID();
		cout << local_time() << "Executing job " << _id.str() << endl;
		const auto job_path = jobs_path / _id.str();
		const auto format = job["format"].String();
		const auto email = job["email"].String();

		// Parse the user-supplied ligand.
		OBMol obMol;
		OBConversion obConversion;
		obConversion.SetInFormat(format.c_str());
		obConversion.ReadFile(&obMol, (job_path / ("ligand." + format)).string());
		const auto num_atoms = obMol.NumAtoms();
//		obMol.AddHydrogens(); // Adding hydrogens does not seem to affect SMARTS matching.

		// Classify subset atoms.
		array<vector<int>, num_subsets> subsets;
		for (size_t k = 0; k < num_subsets; ++k)
		{
			auto& subset = subsets[k];
			subset.reserve(num_atoms);
			OBSmartsPattern smarts;
			smarts.Init(SubsetSMARTS[k]);
			smarts.Match(obMol);
			for (const auto& map : smarts.GetMapList())
			{
				subset.push_back(map.front());
			}
		}
		const auto& subset0 = subsets.front();

		// Check user-provided ligand validity.
		if (subset0.empty())
		{
			// Record job completion time stamp.
			const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
			conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));

			// Send error notification email.
			cout << local_time() << "Sending an error notification email to " << email << endl;
			MailMessage message;
			message.setSender("usr <noreply@cse.cuhk.edu.hk>");
			message.setSubject("Your usr job has failed");
			message.setContent("Description: " + job["description"].String() + "\nSubmitted: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC\nFailed: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC\nReason: failed to parse the provided ligand.");
			message.addRecipient(MailRecipient(MailRecipient::PRIMARY_RECIPIENT, email));
			SMTPClientSession session("137.189.91.190");
			session.login();
			session.sendMessage(message);
			session.close();
			continue;
		}

		// Calculate the four reference points.
		const auto n = subset0.size();
		const auto v = 1.0 / n;
		array<vector3, num_references> references{};
		auto& ctd = references[0];
		auto& cst = references[1];
		auto& fct = references[2];
		auto& ftf = references[3];
		for (const auto i : subset0)
		{
			ctd += obMol.GetAtom(i)->GetVector();
		}
		ctd *= v;
		double cst_dist = numeric_limits<double>::max();
		double fct_dist = numeric_limits<double>::lowest();
		double ftf_dist = numeric_limits<double>::lowest();
		for (const auto i : subset0)
		{
			const auto& a = obMol.GetAtom(i)->GetVector();
			const auto this_dist = a.distSq(ctd);
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
			const auto& a = obMol.GetAtom(i)->GetVector();
			const auto this_dist = a.distSq(fct);
			if (this_dist > ftf_dist)
			{
				ftf = a;
				ftf_dist = this_dist;
			}
		}

		// Precalculate the distances between each atom and each reference point.
		array<vector<double>, num_references> dista;
		for (size_t k = 0; k < num_references; ++k)
		{
			const auto& reference = references[k];
			auto& dists = dista[k];
			dists.resize(1 + num_atoms); // OpenBabel atom index starts from 1. dists[0] is dummy.
			for (size_t i = 0; i < n; ++i)
			{
				dists[subset0[i]] = sqrt(obMol.GetAtom(subset0[i])->GetVector().distSq(reference));
			}
		}

		// Calculate USR and USRCAT features of the input ligand.
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
				#pragma unroll
				for (const auto e : m)
				{
					q[qo++] = e;
				}
			}
		}
		assert(qo == qn.back());

		// Compute USR and USRCAT scores.
		usrcat_bin.seekg(0);
		for (size_t k = 0; k < num_ligands; ++k)
		{
			usrcat_bin.read(reinterpret_cast<char*>(l.data()), sizeof(l));
			double s = 0;
			#pragma unroll
			for (size_t i = 0, u = 0; u < num_usrs; ++u)
			{
				const auto qnu = qn[u];
				#pragma unroll
				for (; i < qnu; ++i)
				{
					s += fabs(q[i] - l[i]);
				}
				scores[u][k] = 1 / (1 + s * qv[u]); // TODO: use s for sorting.
			}
		}
		assert(usrcat_bin.tellg() == sizeof(l) * num_ligands);

		// Sort ligands by USRCAT score.
		const size_t u = 1;
		const auto& uscores = scores[u];
		iota(scase.begin(), scase.end(), 0);
		sort(scase.begin(), scase.end(), [&uscores](const size_t val1, const size_t val2)
		{
			return uscores[val1] > uscores[val2];
		});

		// Write results.
		filtering_ostream log_csv_gz;
		log_csv_gz.push(gzip_compressor());
		log_csv_gz.push(file_sink((job_path / "log.csv.gz").string()));
		log_csv_gz.setf(ios::fixed, ios::floatfield);
		log_csv_gz << "ZINC ID,USR score,USRCAT score\n" << setprecision(8);
		filtering_ostream ligands_pdbqt_gz;
		ligands_pdbqt_gz.push(gzip_compressor());
		ligands_pdbqt_gz.push(file_sink((job_path / "ligands.pdbqt.gz").string()));
		ligands_pdbqt_gz.setf(ios::fixed, ios::floatfield);
		for (size_t t = 0; t < 10000; ++t)
		{
			const size_t k = scase[t];
			const auto zincid = zincids.substr(9 * k, 8);
			const auto u0score = scores[0][k];
			const auto u1score = scores[1][k];
			log_csv_gz << zincid << ',' << u0score << ',' << u1score << '\n';

			// Only write conformations of the top ligands to ligands.pdbqt.gz.
			if (t >= 1000) continue;

			const auto zp = zproperties[k];
			ligands_pdbqt_gz
				<< "MODEL " << '\n'
				<< "REMARK 911 " << zincid
				<< setprecision(3)
				<< ' ' << setw(8) << zp.mwt
				<< ' ' << setw(8) << zp.lgp
				<< ' ' << setw(8) << zp.ads
				<< ' ' << setw(8) << zp.pds
				<< ' ' << setw(3) << zp.hbd
				<< ' ' << setw(3) << zp.hba
				<< ' ' << setw(3) << zp.psa
				<< ' ' << setw(3) << zp.chg
				<< ' ' << setw(3) << zp.nrb
				<< '\n'
				<< "REMARK 912 " << smileses[k] << '\n'
				<< "REMARK 913 " << suppliers[k] << '\n'
				<< setprecision(8)
				<< "REMARK 951    USR SCORE: " << setw(10) << u0score << '\n'
				<< "REMARK 952 USRCAT SCORE: " << setw(10) << u1score << '\n'
			;
			ligand_pdbqt.seekg(headers[k]);
			const size_t length = headers[k + 1] - headers[k];
			if (length > buf.size()) buf.resize(length);
			ligand_pdbqt.read(const_cast<char*>(buf.data()), length);
			ligands_pdbqt_gz.write(buf.data(), length);
			ligands_pdbqt_gz << "ENDMDL\n";
		}

		// Update progress.
		cout << local_time() << "Setting done time" << endl;
		const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
		conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));

		// Send completion notification email.
		cout << local_time() << "Sending a completion notification email to " << email << endl;
		MailMessage message;
		message.setSender("istar <noreply@cse.cuhk.edu.hk>");
		message.setSubject("Your usr job has completed");
		message.setContent("Description: " + job["description"].String() + "\nSubmitted: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC\nCompleted: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC\nResult: http://istar.cse.cuhk.edu.hk/usr/iview/?" + _id.str());
		message.addRecipient(MailRecipient(MailRecipient::PRIMARY_RECIPIENT, email));
		SMTPClientSession session("137.189.91.190");
		session.login();
		session.sendMessage(message);
		session.close();
	}
}
