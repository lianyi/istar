#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <immintrin.h>
/*#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <mongo/client/dbclient.h>
#include <Poco/Net/MailMessage.h>
#include <Poco/Net/MailRecipient.h>
#include <Poco/Net/SMTPClientSession.h>*/
using namespace std;
/*using namespace std::chrono;
using namespace boost::filesystem;
using namespace boost::iostreams;
using namespace boost::gregorian;
using namespace boost::posix_time;
using namespace mongo;
using namespace bson;
using namespace Poco::Net;*/

/*inline static string now()
{
	return to_simple_string(second_clock::local_time()) + " ";
}*/

vector<double> parse(const string& line)
{
	vector<double> r;
	if (line.size())
	{
		r.reserve(12);
		for (size_t b = 0, e; true; b = e + 1)
		{
			if ((e = line.find(',', b + 6)) == string::npos)
			{
				r.push_back(stof(line.substr(b)));
				break;
			}
			r.push_back(stof(line.substr(b, e - b)));
		}
	}
	return r;
}

int main(int argc, char* argv[])
{
/*	// Check the required number of command line arguments.
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
	const auto collection = "istar.usr";*/

	const auto m256s = _mm256_set1_pd(-0. ); // -0.  = 1 << 63

	// Read the feature file.
	string line;
	vector<vector<double>> features;
	features.reserve(23129083);
	for (ifstream ifs(argv[1]); getline(ifs, line); features.push_back(parse(line)));
	const size_t n = features.size();

	// Read the header file.
	vector<string> headers;
	headers.reserve(n);
	for (ifstream ifs(argv[2]); getline(ifs, line); headers.push_back(move(line)));

	// Search the features for records similar to each of the queries.
	vector<double> scores(n);
	vector<size_t> scase(n);
	array<double, 4> a;
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(4);
	const double qv = 1.0 / 12;
	while (getline(cin, line))
	{
		const auto& q = parse(line);
		for (size_t k = 0; k < n; ++k)
		{
			const auto& l = features[k];
			double s = 0;
			for (size_t i = 0; i < 12; i += 4)
			{
				_mm256_stream_pd(a.data(), _mm256_andnot_pd(m256s, _mm256_sub_pd(_mm256_load_pd(&q[i]), _mm256_load_pd(&l[i]))));
				s += a[0] + a[1] + a[2] + a[3];
			}
			scores[k] = 1 / (1 + s * qv);
		}
		iota(scase.begin(), scase.end(), 0);
		sort(scase.begin(), scase.end(), [&scores](const size_t val1, const size_t val2)
		{
			return scores[val1] > scores[val2];
		});
		for (const size_t c : scase)
		{
			cout << c << '\t' << headers[c] << '\t' << scores[c] << endl;
		}
	}

/*	// Initialize epoch
	const auto epoch = date(1970, 1, 1);

	while (true)
	{
		// Fetch jobs.
		auto cursor = conn.query(collection, QUERY("done" << BSON("$exists" << false)).sort("submitted"), 100); // Each batch processes 100 jobs.
		while (cursor->more())
		{
			const auto job = cursor->next();
			const auto _id = job["_id"].OID();
			cout << now() << "Executing job " << _id.str() << endl;

			// Obtain the target genome via taxid.
			const auto ligand = job["ligand"].Int();

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
			message.setSender("igrep <noreply@cse.cuhk.edu.hk>");
			message.setSubject("Your igrep job has completed");
			message.setContent("Your igrep job submitted on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC searching the genome of " + g.name + " for " + to_string(qi) + " patterns was done on " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC. View result at http://istar.cse.cuhk.edu.hk/igrep");
			message.addRecipient(MailRecipient(MailRecipient::PRIMARY_RECIPIENT, email));
			SMTPClientSession session("137.189.91.190");
			session.login();
			session.sendMessage(message);
			session.close();
		}

		// Sleep for a while.
		this_thread::sleep_for(std::chrono::seconds(10));
	}*/
}
