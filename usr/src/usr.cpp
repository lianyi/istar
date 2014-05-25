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
	std::ifstream ifs(f, std::ifstream::binary);
	ifs.seekg(0, std::ifstream::end);
	const size_t num_bytes = ifs.tellg();
	v.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(v.data()), num_bytes);
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
	const double qv = 1.0 / 12;
	const auto epoch = date(1970, 1, 1);

	// Read the feature bin file.
	vector<array<double, 12>> features;
	read(features, "16_usr.bin");
	const size_t n = features.size();

	// Read the header bin file.
	vector<size_t> headers;
	read(headers, "16_hdr.bin");
	assert(n == headers.size());

	// Search the features for records similar to the query.
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(4);
	vector<double> scores(n);
	vector<size_t> scase(n);
	array<double, 4> a;
	string line;
	std::ifstream ligands("16_lig.pdbqt");
	while (true)
	{
		// Fetch jobs.
		auto cursor = conn.query(collection, QUERY("done" << BSON("$exists" << false)).sort("submitted"), 100); // Each batch processes 100 jobs.
		while (cursor->more())
		{
			const auto job = cursor->next();
			const auto _id = job["_id"].OID();
			const auto job_path = jobs_path / _id.str();
			cout << now() << "Executing job " << _id.str() << endl;

			// Obtain job properties.
			const auto ligand = job["ligand"].Array();
			array<double, 12> q;
			for (size_t i = 0; i < 12; ++i)
			{
				q[i] = stod(ligand[i].String());
			}

			for (size_t k = 0; k < n; ++k)
			{
				const auto& l = features[k];
				double s = 0;
				for (size_t i = 0; i < 12; i += 4)
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
			boost::filesystem::ofstream log_csv(job_path / "log.csv.gz");
			boost::filesystem::ofstream ligands_pdbqt(job_path / "ligands.pdbqt.gz");
			filtering_ostream log_csv_gz;
			filtering_ostream ligands_pdbqt_gz;
			log_csv_gz.push(gzip_compressor());
			ligands_pdbqt_gz.push(gzip_compressor());
			log_csv_gz.push(log_csv);
//			log_csv_gz.push(file_sink((job_path / "log.csv.gz").string()));
			ligands_pdbqt_gz.push(ligands_pdbqt);
			log_csv_gz.setf(ios::fixed, ios::floatfield);
			log_csv_gz << "ZINC ID,USR score,Molecular weight (g/mol),Partition coefficient xlogP,Apolar desolvation (kcal/mol),Polar desolvation (kcal/mol),Hydrogen bond donors,Hydrogen bond acceptors,Polar surface area tPSA (A^2),Net charge,Rotatable bonds,SMILES,Substance information,Suppliers\n" << setprecision(3);
			for (const size_t c : scase)
			{
//				ligands.seekg(headers[c]);
				getline(ligands, line);
//				cout << c << '\t' << scores[c] << endl;
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
